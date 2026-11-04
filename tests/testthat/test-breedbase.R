# A small fake Breedbase server for testing the BrAPI client without a
# network. It keeps state (folders, studies, plots) and checks the same
# things the real server checks (see CXGN::BrAPI::v2 in solgenomics/sgn).
reset_found_cache <- function() rm(list = ls(breedtools:::bb_found_cache), envir = breedtools:::bb_found_cache)

fake_breedbase <- function(accessions = c(paste0("G", 1:200), paste0("C", 1:6), LETTERS[1:4]),
                           existing_studies = character(0), password = "secret", max_search = Inf,
                           slow_writes = character(0), lost_writes = character(0)) {
  reset_found_cache()
  state <- new.env()
  state$folders <- list(list(trialDbId = "500", trialName = "Old folder", programDbId = "10"))
  state$studies <- lapply(existing_studies, function(n) list(studyDbId = "900", studyName = n))
  state$plots <- list()
  state$observations <- list()
  state$variables <- list(list(observationVariableDbId = "70001", observationVariableName = "fresh root yield|CO_334:0000013"),
                          list(observationVariableDbId = "70002", observationVariableName = "dry matter content percentage|CO_334:0000092"),
                          list(observationVariableDbId = "70003", observationVariableName = "fresh root yield spatially adjusted|COMP:0000123"))
  state$analyses <- list()
  # SGNSTAT terms and composed traits (by component ids), like the chado cvterm tables.
  state$statistics <- list(list(cvterm_id = 80054, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted Means from Spatial Correction using SpATS R",
                                accession = "SGNSTAT:0000054"),
                           list(cvterm_id = 80036, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted means from BLUEs using LMER R", accession = "SGNSTAT:0000036"))
  state$composed <- list()
  state$requests <- list()
  state$next_id <- 1000
  # slow_writes: the write is done but the gateway answers 504 (once).
  # lost_writes: the gateway answers 504 and the write is not done (once).
  state$slow <- slow_writes
  state$lost <- lost_writes
  ok <- function(data = NULL, extra = list()) {
    httr2::response_json(200, body = c(list(metadata = list(pagination = list(totalPages = 1), status = list()),
                                     result = list(data = data)), extra))
  }
  err <- function(code, msg) {
    httr2::response_json(code, body = list(metadata = list(status = list(list(message = msg, messageType = "ERROR"))),
                                    result = list()))
  }
  mock <- function(req) {
    url <- httr2::url_parse(req$url)
    path <- sub("^/brapi/v2/", "", url$path)
    method <- req$method %||% if (!is.null(req$body)) "POST" else "GET"
    body <- req$body$data
    auth <- httr2::req_get_headers(req, redacted = "reveal")$Authorization
    state$requests[[length(state$requests) + 1]] <- list(method = method, path = path, body = body, auth = auth)
    if (path == "token") {
      if (identical(body$password, password)) {
        return(httr2::response_json(200, body = list(access_token = "TOKEN123", userDisplayName = "Jane Doe")))
      }
      return(httr2::response_json(200, body = list(access_token = NULL,
                                            metadata = list(status = list(list(message = "Incorrect password"))))))
    }
    # Breedbase AJAX (not BrAPI): login token as sgn_session_id; errors as
    # {"error": ...} with HTTP 200, like CXGN's controllers.
    if (grepl("^/?ajax/", path)) {
      path <- sub("^/", "", path)
      form <- if (!is.null(req$body)) lapply(req$body$data, function(x) utils::URLdecode(as.character(x))) else list()
      session <- form$sgn_session_id %||% url$query$sgn_session_id
      if (!identical(session, "TOKEN123")) return(httr2::response_json(200, body = list(error = "You must be logged in to do this!")))
      if (path == "ajax/onto/match") {
        q <- url$query
        hits <- Filter(function(s) grepl(tolower(q$term_name), tolower(s$cvterm_name), fixed = TRUE) &&
                         startsWith(s$accession, paste0(q$db_name, ":")), state$statistics)
        return(httr2::response_json(200, body = hits))
      }
      if (path == "ajax/onto/get_trait_from_exact_components") {
        ids <- regmatches(req$url, gregexpr("ids(%5B%5D|\\[\\])=[0-9]+", req$url))[[1]]
        ids <- sort(sub("^.*=", "", ids))
        hit <- Filter(function(cp) identical(sort(cp$components), ids), state$composed)
        if (!length(hit)) return(httr2::response_json(200, body = list(error = "No exact matches found.")))
        return(httr2::response_json(200, body = list(trait_id = as.integer(hit[[1]]$id))))
      }
      if (path == "ajax/onto/store_composed_term") {
        new <- jsonlite::fromJSON(form$new_trait_names, simplifyVector = FALSE)
        for (nm in names(new)) {
          comps <- strsplit(new[[nm]], ",")[[1]]
          if (length(comps) < 2) return(httr2::response_json(200, body = list(error = "less than 2 components")))
          id <- as.character(state$next_id <- state$next_id + 1)
          state$composed[[length(state$composed) + 1]] <- list(id = id, name = nm, components = comps)
        }
        return(httr2::response_json(200, body = list(success = "Saved new trait", names = as.list(names(new)))))
      }
      if (path == "ajax/analysis/store/json") {
        if (form$analysis_name %in% c(vapply(state$analyses, `[[`, "", "name"), vapply(state$studies, `[[`, "", "studyName"))) {
          return(httr2::response_json(200, body = list(error = paste("An analysis with name", form$analysis_name, "already exists"))))
        }
        values <- jsonlite::fromJSON(form$analysis_result_values, simplifyVector = FALSE)
        traits <- unlist(jsonlite::fromJSON(form$analysis_trait_names))
        known <- vapply(state$variables, `[[`, "", "observationVariableName")
        if (!all(traits %in% known)) return(httr2::response_json(200, body = list(error = "trait not found")))
        id <- as.character(state$next_id <- state$next_id + 1)
        state$analyses[[length(state$analyses) + 1]] <- list(id = id, name = form$analysis_name, form = form,
                                                             values = values)
        return(httr2::response_json(200, body = list(success = 1, analysis_id = as.integer(id), model_id = 5)))
      }
      if (path == "ajax/analyses/by_user") {
        rows <- lapply(state$analyses, function(a) list(paste0('<a href="/analyses/', a$id, '">', a$name, "</a>"), "d"))
        return(httr2::response_json(200, body = list(data = rows)))
      }
      if (grepl("^ajax/analysis/[0-9]+/retrieve$", path)) {
        id <- sub("^ajax/analysis/([0-9]+)/retrieve$", "\\1", path)
        a <- Filter(function(a) a$id == id, state$analyses)[[1]]
        rows <- c(list(list("header")), lapply(names(a$values), function(n) list(n)))
        return(httr2::response_json(200, body = list(analysis_name = a$name, data = rows)))
      }
      return(httr2::response_json(404, body = list(error = "unknown ajax call")))
    }
    if (!identical(auth, "Bearer TOKEN123")) return(err(401, "You must be logged in"))
    if (path == "programs") return(ok(list(list(programDbId = "10", programName = "IITA"),
                                         list(programDbId = "11", programName = "NRCRI"))))
    if (path == "locations") return(ok(list(list(locationDbId = "1", locationName = "Ibadan"),
                                          list(locationDbId = "2", locationName = "Ubiaja"),
                                          list(locationDbId = "3", locationName = "Ikenne"))))
    if (path == "trials" && method == "GET") {
      return(ok(Filter(function(f) f$programDbId == url$query$programDbId, state$folders)))
    }
    if (path == "trials" && method == "POST") {
      taken <- c(vapply(state$folders, `[[`, "", "trialName"), vapply(state$studies, `[[`, "", "studyName"))
      if (body[[1]]$trialName %in% taken) {
        return(err(409, sprintf("A trial or study with the name %s already exists in the database.", body[[1]]$trialName)))
      }
      id <- as.character(state$next_id <- state$next_id + 1)
      state$folders[[length(state$folders) + 1]] <- list(trialDbId = id, trialName = body[[1]]$trialName,
                                                         programDbId = body[[1]]$programDbId)
      return(ok(list(list(trialDbId = id))))
    }
    if (path == "search/germplasm") {
      # A slow server: too many names in one search ends in a gateway timeout.
      if (length(unlist(body$germplasmNames)) > max_search) {
        return(httr2::response(504, headers = list(`Content-Type` = "text/html"), body = charToRaw("<html>Gateway Timeout</html>")))
      }
      hits <- intersect(unlist(body$germplasmNames), accessions)
      return(ok(lapply(hits, function(n) list(germplasmDbId = as.character(5000 + match(n, accessions)),
                                              germplasmName = n))))
    }
    if (path == "search/observationunits") {
      hits <- Filter(function(u) u$observationUnitName %in% unlist(body$observationUnitNames), state$plots)
      return(ok(hits))
    }
    if (path == "search/variables") {
      wanted <- unlist(body$observationVariableNames)
      hits <- Filter(function(v) sub("\\|.*$", "", v$observationVariableName) %in% wanted, state$variables)
      return(ok(hits))
    }
    if (path == "search/observations") {
      hits <- Filter(function(o) o$observationUnitDbId %in% unlist(body$observationUnitDbIds) &&
                       o$observationVariableDbId %in% unlist(body$observationVariableDbIds), state$observations)
      return(ok(hits))
    }
    if (path == "observations" && method == "POST") {
      # Like Breedbase: the three fields are required, the plot and trait must
      # exist, and a plot-trait-time combination is stored once.
      unit_ids <- vapply(state$plots, `[[`, "", "observationUnitDbId")
      var_ids <- c(vapply(state$variables, `[[`, "", "observationVariableDbId"),
                   vapply(state$composed, `[[`, "", "id"))
      seen <- vapply(state$observations, function(o) paste(o$observationUnitDbId, o$observationVariableDbId,
                                                            o$observationTimeStamp %||% ""), "")
      for (o in body) {
        if (is.null(o$observationUnitDbId) || is.null(o$observationVariableDbId) || is.null(o$value)) {
          return(err(400, "An observation object is missing a required field"))
        }
        if (!o$observationUnitDbId %in% unit_ids) return(err(400, "observationUnitDbId not found"))
        if (!o$observationVariableDbId %in% var_ids) return(err(400, "observationVariableDbId not found"))
        key <- paste(o$observationUnitDbId, o$observationVariableDbId, o$observationTimeStamp %||% "")
        if (key %in% seen) return(err(400, "The combination already exists in the database"))
        seen <- c(seen, key)
      }
      state$observations <- c(state$observations, body)
      return(ok(list()))
    }
    if (path == "search/studies") {
      hits <- Filter(function(s) s$studyName %in% unlist(body$studyNames), state$studies)
      return(ok(hits))
    }
    if (path == "studies" && method == "POST") {
      s <- body[[1]]
      allowed <- c("CRD", "Alpha", "MAD", "Lattice", "Augmented", "RCBD", "p-rep", "splitplot", "stripplot",
                   "greenhouse", "Westcott", "Analysis")
      if (is.null(s$trialDbId) || is.null(s$locationDbId) || is.null(s$studyType) || is.null(s$experimentalDesign$PUI)) {
        return(err(400, "missing required field"))
      }
      if (!s$experimentalDesign$PUI %in% allowed) return(err(400, "Experimental Design not supported"))
      if (s$studyName %in% vapply(state$studies, `[[`, "", "studyName")) return(err(409, "Study exists"))
      id <- as.character(state$next_id <- state$next_id + 1)
      state$studies[[length(state$studies) + 1]] <- c(s, list(studyDbId = id))
      return(ok(list(list(studyDbId = id))))
    }
    if (path == "observationunits" && method == "POST") {
      nums <- vapply(body, function(u) u$observationUnitPosition$observationLevel$levelCode, "")
      if (anyDuplicated(paste(vapply(body, `[[`, "", "studyDbId"), nums))) return(err(422, "Plot number duplicated"))
      for (u in body) {
        # Like Breedbase: the germplasmDbId is required and must match.
        id <- u$germplasmDbId %||% ""
        if (!grepl("^[0-9]+$", id)) return(err(404, "Verify Germplasm Id."))
        if (!identical(accessions[as.integer(id) - 5000], u$germplasmName)) return(err(404, "germplasm not found"))
        state$next_id <- state$next_id + 1
        u$observationUnitDbId <- as.character(state$next_id)
        state$plots[[length(state$plots) + 1]] <- u
      }
      return(ok(list()))
    }
    if (path == "observationunits" && method == "GET") {
      return(ok(Filter(function(u) u$studyDbId == url$query$studyDbId, state$plots)))
    }
    err(404, paste("unknown call", path))
  }
  timeout <- httr2::response(504, headers = list(`Content-Type` = "text/html"),
                             body = charToRaw("<html>Gateway Timeout</html>"))
  mock_with_timeouts <- function(req) {
    path <- sub("^/brapi/v2/", "", httr2::url_parse(req$url)$path)
    is_write <- !is.null(req$body) && !grepl("^search|^token", path)
    if (is_write && path %in% state$lost) {
      state$lost <- setdiff(state$lost, path)
      return(timeout)
    }
    res <- mock(req)
    if (is_write && path %in% state$slow) {
      state$slow <- setdiff(state$slow, path)
      return(timeout)
    }
    res
  }
  list(mock = mock_with_timeouts, state = state)
}

test_that("dry run without a connection returns the BrAPI records", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  payload <- save_breedbase(rc, trial_name = "AYT_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial")
  expect_length(payload$studies, 1)
  study <- payload$studies[[1]]
  expect_equal(study$studyName, "AYT_2026")
  expect_equal(study$experimentalDesign$PUI, "Augmented")
  expect_equal(study$seasons[[1]], "2026")
  plots <- payload$plots[[1]]
  expect_length(plots, 144)
  expect_equal(sum(vapply(plots, function(p) p$observationUnitPosition$entryType == "CHECK", logical(1))), 36)
  p1 <- plots[[1]]
  first <- rc$design[rc$design$plots == 1, ]
  expect_equal(p1$observationUnitName,
               paste("AYT_2026", first$all_entries, first$block, first$rep, sep = "_"))
  names_all <- vapply(plots, `[[`, "", "observationUnitName")
  expect_false(anyDuplicated(names_all) > 0)
  expect_equal(p1$observationUnitPosition$observationLevel$levelName, "plot")
  expect_equal(p1$observationUnitPosition$positionCoordinateY, as.character(rc$design$row[rc$design$plots == 1]))
  expect_equal(vapply(p1$observationUnitPosition$observationLevelRelationships, `[[`, "", "levelName"), c("rep", "block"))
  expect_error(save_breedbase(rc, trial_name = "x", year = 2026, location = "Ibadan",
  description = "Test trial"), "breedbase_login")
  expect_error(save_breedbase(rc, trial_name = "x", year = 2026, location = "Ibadan", design_type = "RRC",
                              dry_run = TRUE,
                              description = "Test trial"), "Breedbase accepts")
})

test_that("p-rep designs become one trial per location; fillers are skipped", {
  prep <- augmented_prep_design(paste0("G", 1:42), c("North", "South"), block_size = 6, field_cols = 7, seed = 1)
  payload <- save_breedbase(prep, trial_name = "PREP_2026", year = 2026,
                            location = c(North = "Ibadan", South = "Ubiaja"), dry_run = TRUE,
                            description = "Test trial")
  expect_equal(vapply(payload$studies, `[[`, "", "studyName"), c("PREP_2026_North", "PREP_2026_South"))
  expect_equal(vapply(payload$plots, length, integer(1)), as.vector(table(prep$design$location)[c("North", "South")]))
  expect_equal(payload$studies[[1]]$experimentalDesign$PUI, "p-rep")
  nums <- vapply(payload$plots[[2]], function(p) p$observationUnitPosition$observationLevel$levelCode, "")
  expect_false(anyDuplicated(nums) > 0)

  rrc <- resolvable_row_column_design(paste0("P", 1:19), 2, 4, 5, n_starts = 1, n_iter = 200, seed = 1)
  rp <- save_breedbase(rrc, trial_name = "RRC_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial")
  expect_length(rp$plots[[1]], 38)
  expect_false("FILLER" %in% vapply(rp$plots[[1]], `[[`, "", "germplasmName"))
  expect_equal(rp$studies[[1]]$experimentalDesign$PUI, "Alpha")
})

test_that("login, checks and saving work against a fake Breedbase", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase()
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    expect_error(breedbase_login("https://bb.example.org", "jdoe", password = "wrong"), "Login failed")
    con <- suppressMessages(breedbase_login("bb.example.org/", "jdoe", password = "secret"))
    expect_s3_class(con, "breedbase_connection")
    expect_equal(con$url, "https://bb.example.org")
    expect_false(any(grepl("secret", capture.output(print(con)))))
    expect_null(con$password)

    expect_equal(breedbase_programs(con)$program, c("IITA", "NRCRI"))
    expect_equal(breedbase_folders(con, "IITA")$folder, "Old folder")
    expect_error(breedbase_folders(con, "Nope"), "Programs available: IITA, NRCRI")

    # Dry run with a connection checks everything and writes nothing.
    n_before <- length(fake$state$requests)
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, dry_run = TRUE,
                                                description = "Test trial"))
    expect_true(any(grepl("Nothing was written", out)))
    posts <- Filter(function(r) r$method == "POST" && !grepl("^search|token", r$path), fake$state$requests)
    expect_length(posts, 0)

    # Real save.
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, confirm = FALSE,
                                                description = "Test trial"))
    expect_equal(res$plots_sent, 144)
    expect_equal(res$plots_found, 144)
    expect_match(res$trial_page, "^https://bb.example.org/breeders/trial/")
    expect_true("AYT 2026" %in% vapply(fake$state$folders, `[[`, "", "trialName"))
    study <- fake$state$studies[[1]]
    expect_equal(study$locationDbId, "1")
    expect_equal(study$trialDbId, fake$state$folders[[2]]$trialDbId)
    expect_true(all(vapply(Filter(function(r) r$path != "token", fake$state$requests),
                           function(r) identical(r$auth, "Bearer TOKEN123"), logical(1))))

    # Saving again with the same name is refused before writing.
    expect_error(capture.output(save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                               trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                               confirm = FALSE,
                                               description = "Test trial")), "already exists")
  })
})

test_that("clear errors for missing accessions, folders, locations and logins", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = paste0("G", 1:100))
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial")),
                 "12 of 112 entries/checks are not accessions")
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Missing", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial")),
                 "create_folder = TRUE")
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Mars", confirm = FALSE,
    description = "Test trial")),
                 "Location 'Mars' was not found")
    bad <- con
    bad$token <- "expired"
    expect_error(breedbase_programs(bad), "HTTP 401.*log in again")
  })
  expect_false(any(grepl("secret", unlist(lapply(fake$state$requests, function(r) r$auth)))))
})

test_that("plot names: automatic pattern, own pattern, and user-given names", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                      plot_name_pattern = "{trial}-{plot}-r{row}c{col}",
                      description = "Test trial")
  expect_equal(p$plots[[1]][[1]]$observationUnitName, sprintf("T-1-r%dc%d", d$row[d$plots == 1], d$col[d$plots == 1]))
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_name_pattern = "{trial}_{plant}",
                              description = "Test trial"), "Unknown field")

  own <- paste0("MYPLOT-", seq_len(nrow(d)))
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE, plot_names = own,
  description = "Test trial")
  sent <- vapply(p$plots[[1]], `[[`, "", "observationUnitName")
  expect_setequal(sent, own)
  d$my_names <- paste0("X", d$plots)
  rc2 <- rc
  rc2$design <- d
  p <- save_breedbase(rc2, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE, plot_names = "my_names",
  description = "Test trial")
  expect_equal(p$plots[[1]][[1]]$observationUnitName, "X1")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_names = rep("A", nrow(d)),
                              description = "Test trial"), "unique")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_names = "A",
                              description = "Test trial"), "one name per plot")

  # Extra check replication: the same check twice in a block gets the plot number added.
  extra <- augmented_row_column_design(paste0("E", 1:90), LETTERS[1:4], 12, 12, 4, 4,
                                       checks_per_block = 6, n_candidates = 2, seed = 1)
  expect_message(p <- save_breedbase(extra, trial_name = "X", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial"),
                 "plot number was added")
  sent <- vapply(p$plots[[1]], `[[`, "", "observationUnitName")
  expect_false(anyDuplicated(sent) > 0)
})

test_that("plot names already in the database are refused before writing", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase()
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial"))
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T2", 2026, "Ibadan",
                                               plot_names = vapply(fake$state$plots, `[[`, "", "observationUnitName")[
                                                 order(vapply(fake$state$plots, function(u) as.integer(u$observationUnitPosition$observationLevel$levelCode), 1L))],
                                               confirm = FALSE,
                                               description = "Test trial")),
                 "plot names already exist")
  })
})

test_that("a trial description is required and is sent to Breedbase", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2, seed = 1)
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE),
               "Give a `description` of the trial")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", description = "  ",
                              dry_run = TRUE), "Give a `description`")
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan",
                      description = "Yield trial of 108 lines, Ibadan 2026", dry_run = TRUE)
  expect_equal(p$studies[[1]]$studyDescription, "Yield trial of 108 lines, Ibadan 2026")
})

test_that("breedbase_check reports every problem at once, with suggestions", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = c(paste0("G", 1:108), paste0("C", 1:4), "TMEB419"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))

    out <- capture.output(res <- breedbase_check(con, accessions = c(paste0("G", 1:100), "tmeb419", "X9", "G1"),
                                                 location = "ibadan", program = "IITA", folder = "New one",
                                                 trial_name = "T1"))
    expect_false(attr(res, "ok"))
    expect_equal(res$status[res$what == "location"], "MISSING")
    expect_match(res$suggestion[res$what == "location"], "did you mean 'Ibadan'")
    expect_equal(res$status[res$what == "folder"], "NEW")
    expect_equal(res$status[res$what == "trial"], "OK")
    acc <- res[res$what == "accession", ]
    expect_match(acc$suggestion[acc$name == "tmeb419"], "did you mean 'TMEB419'")
    expect_match(acc$suggestion[acc$name == "X9"], "add it to Breedbase first")
    expect_equal(acc$status[acc$name == "G1"], c("REPEATED", "OK"))
    expect_true(any(grepl("100 of 102 accessions are in Breedbase", out)))
    expect_true(any(grepl("needs? fixing", out)))

    # A design can be given directly; all good -> ok.
    rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                      n_candidates = 5, seed = 1)
    out <- capture.output(res <- breedbase_check(con, rc, location = "Ibadan"))
    expect_true(attr(res, "ok"))
    expect_true(any(grepl("All names are ready", out)))
    # No write calls were made by the checks.
    posts <- Filter(function(r) r$method == "POST" && !grepl("^search|token", r$path), fake$state$requests)
    expect_length(posts, 0)
  })
})

test_that("save_breedbase lists location and accession problems together", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = paste0("G", 1:100))
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    err <- tryCatch(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "ibadan",
                                                  confirm = FALSE, description = "Test trial")),
                    error = conditionMessage)
    expect_match(err, "Location 'ibadan' was not found.*did you mean 'Ibadan'")
    expect_match(err, "12 of 112 entries/checks are not accessions")
  })
})

test_that("slow servers: searches that time out are split and retried", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.retry_pause = 0)
  on.exit(options(old_options), add = TRUE)
  fake <- fake_breedbase(max_search = 12)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, dry_run = TRUE, description = "Test trial"))
    expect_true(any(grepl("all 112 entries and checks exist as accessions", out)))
    expect_true(any(grepl("searched 112 accessions in [0-9]+ batches .*split and retried", out)))
    # Found accessions are remembered: a second check searches nothing.
    n_search <- length(Filter(function(r) r$path == "search/germplasm", fake$state$requests))
    out2 <- capture.output(breedbase_check(con, accessions = rc))
    expect_true(any(grepl("all 112 accessions were already found in this session", out2)))
    expect_equal(length(Filter(function(r) r$path == "search/germplasm", fake$state$requests)), n_search)
    sizes <- vapply(Filter(function(r) r$path == "search/germplasm", fake$state$requests),
                    function(r) length(unlist(r$body$germplasmNames)), integer(1))
    expect_true(all(sizes <= 25))
    expect_true(any(sizes > 12))

    # A server that cannot answer even one name gives a clear message.
    dead <- fake_breedbase(max_search = 0)
    httr2::with_mocked_responses(dead$mock, {
      con2 <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
      expect_error(capture.output(breedbase_check(con2, accessions = c("G1", "G2"))),
                   "HTTP 504.*too slow or busy")
    })
  })
})

test_that("timeouts while writing: the result is checked, retried, and a re-run continues", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.write_wait = 0, breedtools.write_checks = 2, breedtools.retry_pause = 0)
  on.exit(options(old_options), add = TRUE)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  save <- function(con) save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                       trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                       create_folder = TRUE, confirm = FALSE, description = "Test trial")
  n_studies <- function(f) length(f$state$studies)

  # The server finishes every write but the gateway times out: no duplicates.
  fake <- fake_breedbase(slow_writes = c("trials", "studies", "observationunits"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save(con))
  })
  expect_equal(res$plots_found, 144)
  expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "AYT 2026"), 1)
  expect_equal(n_studies(fake), 1)
  expect_length(fake$state$plots, 144)
  expect_true(any(grepl("slow to create folder 'AYT 2026'", out)))
  expect_true(any(grepl("confirmed: done", out)))

  # Writes that are lost are sent once more.
  fake <- fake_breedbase(lost_writes = c("trials", "studies", "observationunits"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save(con))
  })
  expect_equal(res$plots_found, 144)
  expect_equal(n_studies(fake), 1)
  expect_true(any(grepl("trying again", out)))

  # A save cut off after the trial was made: the re-run adds the plots.
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    # Lose the plots on every try.
    lose_all <- function(req) {
      if (grepl("/brapi/v2/observationunits$", req$url) && !is.null(req$body)) return(httr2::response(504))
      fake$mock(req)
    }
    httr2::with_mocked_responses(lose_all, {
      expect_error(capture.output(save(con)), "kept timing out.*run save_breedbase\\(\\) again")
    })
    expect_equal(n_studies(fake), 1)
    expect_length(fake$state$plots, 0)
    out <- capture.output(res <- save(con))
  })
  expect_true(any(grepl("exists without plots", out)))
  expect_equal(res$plots_found, 144)
  expect_equal(n_studies(fake), 1)
})

test_that("an existing trial is reported; an existing folder is used", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.write_wait = 0, breedtools.write_checks = 1)
  on.exit(options(old_options), add = TRUE)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", folder = "Old folder", trial_name = "T_2026",
                                  year = 2026, location = "Ibadan", confirm = FALSE, description = "Test"))
    # The folder existed: no new folder, the trial is inside it.
    expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "Old folder"), 1)
    expect_equal(fake$state$studies[[1]]$trialDbId, "500")

    out <- capture.output(res <- breedbase_check(con, accessions = rc, location = "Ibadan",
                                                 program = "IITA", folder = "Old folder", trial_name = "T_2026"))
    expect_false(attr(res, "ok"))
    expect_equal(res$status[res$what == "trial"], "TAKEN")
    expect_true(any(grepl("\\[FIX\\]  trial 'T_2026' already exists.*choose another trial_name", out)))
    expect_true(any(grepl("folder 'Old folder' exists: the trial will be saved inside it", out)))

    # A folder made by an earlier save that timed out (not listed yet when
    # checked): Breedbase answers 409 and the existing folder is used.
    fake$state$folders[[length(fake$state$folders) + 1]] <- list(trialDbId = "777", trialName = "AYT 2026",
                                                                 programDbId = "10")
    hide_once <- TRUE
    listing_hides_new <- function(req) {
      if (hide_once && grepl("/brapi/v2/trials\\?", req$url) && is.null(req$body)) {
        hide_once <<- FALSE
        return(httr2::response_json(200, body = list(metadata = list(pagination = list(totalPages = 1)),
                                                     result = list(data = list(fake$state$folders[[1]])))))
      }
      fake$mock(req)
    }
    httr2::with_mocked_responses(listing_hides_new, {
      out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                  trial_name = "T2_2026", year = 2026, location = "Ibadan",
                                                  create_folder = TRUE, confirm = FALSE, description = "Test"))
    })
    expect_true(any(grepl("folder 'AYT 2026' already exists: the trial is saved inside it", out)))
    expect_equal(fake$state$studies[[2]]$trialDbId, "777")
    expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "AYT 2026"), 1)

    # A trial left without plots is a note, not a problem.
    fake$state$studies[[3]] <- list(studyDbId = "950", studyName = "EMPTY_2026")
    out <- capture.output(res <- breedbase_check(con, trial_name = "EMPTY_2026"))
    expect_true(attr(res, "ok"))
    expect_true(any(grepl("\\[NOTE\\] trial 'EMPTY_2026' exists without plots", out)))
  })
})

test_that("the folder is optional: the trial goes directly in the program", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", trial_name = "NoFolder_2026",
                                                year = 2026, location = "Ibadan", confirm = FALSE,
                                                description = "Test"))
  })
  expect_true(any(grepl("no folder given: the trial will be saved directly in program 'IITA'", out)))
  expect_false(any(grepl("folder ''", out, fixed = TRUE)))
  expect_equal(fake$state$studies[[1]]$trialDbId, "10")   # the program's id
  expect_length(fake$state$folders, 1)                     # no folder created
  expect_equal(res$plots_found, 144)
  expect_error(save_breedbase(rc, con, program = "IITA", folder = "", trial_name = "X", year = 2026,
                              location = "Ibadan", description = "Test"), "or NULL")
})

test_that("phenotype_template gives the plot names that save_breedbase uses", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  sheet <- phenotype_template(rc, "AYT_2026", traits = c("fresh root yield", "dry matter content percentage"))
  sent <- save_breedbase(rc, trial_name = "AYT_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
                         description = "x")
  expect_setequal(sheet$plot_name, vapply(sent$plots[[1]], `[[`, "", "observationUnitName"))
  expect_equal(nrow(sheet), 144)
  expect_true(all(c("fresh root yield", "dry matter content percentage") %in% names(sheet)))
  expect_equal(sum(sheet$is_check), sum(rc$design$is_control))
  expect_false(anyNA(sheet$row))

  prep <- augmented_prep_design(paste0("E", 1:30), locations = c("North", "South"), seed = 1)
  ps <- phenotype_template(prep, "PREP_2026", traits = "fresh root yield")
  expect_setequal(unique(ps$trial_name), c("PREP_2026_North", "PREP_2026_South"))
})

test_that("phenotypes: dry run, wide and long data, checks and upload", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  sheet <- phenotype_template(rc, "AYT_2026", traits = c("fresh root yield", "dry matter content percentage|CO_334:0000092"))
  set.seed(2)
  sheet$`fresh root yield` <- round(rnorm(144, 25, 4), 1)
  sheet$`dry matter content percentage|CO_334:0000092` <- round(rnorm(144, 32, 3), 1)
  sheet$`fresh root yield`[1:4] <- NA        # not measured: skipped

  obs <- suppressMessages(save_breedbase_phenotypes(sheet, dry_run = TRUE))
  expect_equal(nrow(obs), 2 * 144 - 4)
  expect_setequal(unique(obs$trait), c("fresh root yield", "dry matter content percentage|CO_334:0000092"))
  expect_error(save_breedbase_phenotypes(sheet, plot = "plot"), "Columns are")
  dup <- rbind(sheet[1:2, ], sheet[2, ])
  expect_error(save_breedbase_phenotypes(dup, dry_run = TRUE), "more than one value")

  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    # Plots not saved yet: listed together with a wrong trait.
    bad <- sheet
    names(bad)[names(bad) == "fresh root yield"] <- "root yield"
    err <- tryCatch(capture.output(save_breedbase_phenotypes(bad, con, confirm = FALSE)), error = conditionMessage)
    expect_match(err, "144 of 144 plots are not in Breedbase")
    expect_match(err, "not found in Breedbase: 'root yield'")

    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    n_before <- length(fake$state$requests)
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, dry_run = TRUE))
    expect_true(any(grepl("Nothing was written", out)))
    expect_length(fake$state$observations, 0)

    old_options <- options(breedtools.observation_batch = 100)
    on.exit(options(old_options), add = TRUE)
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    expect_length(fake$state$observations, 2 * 144 - 4)
    expect_equal(res$values_found, res$values_sent)
    expect_equal(res$observation_variable_db_id, c("70001", "70002"))
    expect_true(any(grepl("batch 3 of 3 saved", out)))
    one <- fake$state$observations[[1]]
    expect_true(all(c("observationUnitDbId", "observationVariableDbId", "value") %in% names(one)))
    expect_type(one$value, "character")

    # Saving the same values again is refused by Breedbase (no duplicates).
    expect_error(capture.output(save_breedbase_phenotypes(sheet, con, confirm = FALSE)), "already exists")

    # Long data with time stamps: repeated measurements on other dates.
    long <- data.frame(plot = sheet$plot_name[1:10], trait = "fresh root yield", value = 20:29,
                       when = as.Date("2026-09-01"))
    out <- capture.output(res <- save_breedbase_phenotypes(long, con, plot = "plot", trait = "trait",
                                                           value = "value", timestamp = "when", confirm = FALSE))
    expect_equal(res$values_found, 10)
    expect_equal(fake$state$observations[[length(fake$state$observations)]]$observationTimeStamp,
                 "2026-09-01 00:00:00")
  })
})

test_that("first-stage plot values: one per field plot, design effects removed", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  d$dm <- rnorm(nrow(d), 32, 2)
  fit <- analyze_augmented_row_column_design(d, "yield")
  v <- first_stage_values(list("fresh root yield" = fit,
                               "dry matter content percentage" = analyze_augmented_row_column_design(d, "dm")),
                          rc, "AYT_2026")
  expect_equal(nrow(v), 144)
  expect_setequal(v$plot_name, phenotype_template(rc, "AYT_2026")$plot_name)
  # value = adjusted mean of the genotype + residual of the plot
  pr <- fit$plot_results
  am <- setNames(fit$adjusted_means$adjusted_mean, fit$adjusted_means$genotype)
  i <- match(v$plot_number, pr$plots)
  expect_equal(v$`fresh root yield`, unname(am[pr$all_entries[i]] + pr$residual[i]))
  # Row groups are removed (the trend within groups needs a spatial model).
  expect_lt(abs(coef(lm(v$`fresh root yield` ~ v$row))[2]), 0.2)
  expect_equal(attr(v, "method"), "mixed")
  if (requireNamespace("SpATS", quietly = TRUE)) {
    vs <- first_stage_values(analyze_augmented_row_column_design(d, "yield", spatial = TRUE), rc, "AYT_2026",
                             trait = "fresh root yield")
    expect_lt(abs(coef(lm(vs$`fresh root yield` ~ vs$row))[2]), 0.15)   # most of the row trend (0.4) is removed
    expect_equal(attr(vs, "method"), "spatial")
  }
  expect_error(first_stage_values(fit, rc, "AYT_2026"), "Give the Breedbase `trait`")

  # p-rep two-stage: stage-1 fits per location, one trial per location
  prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
  p <- prep$design
  p$y <- rnorm(nrow(p), 25, 3)
  v2 <- first_stage_values(analyze_augmented_prep_two_stage(p, "y"), prep, "PREP_2026", trait = "fresh root yield")
  expect_equal(nrow(v2), nrow(p))
  expect_false(anyNA(v2$`fresh root yield`))
  expect_setequal(unique(v2$trial_name), c("PREP_2026_North", "PREP_2026_South"))
})

test_that("resolvable row-column analyses give plot residuals", {
  r <- resolvable_row_column_design(paste0("G", 1:24), n_rep = 2, rows = 4, cols = 6, seed = 1)
  dd <- r$design
  set.seed(3)
  dd$y <- 20 + rnorm(nrow(dd)) + 0.5 * dd$row
  f <- analyze_resolvable_row_column(dd, "y", compare = FALSE)
  pr <- f$plot_results
  expect_equal(nrow(pr), nrow(dd))
  expect_equal(pr$fitted + pr$residual, pr$y)
  # X'e = 0: residuals sum to zero within each replicate and each treatment
  expect_true(all(abs(tapply(pr$residual, pr$replicate, sum)) < 1e-6))
  expect_true(all(abs(tapply(pr$residual, pr$all_entries, sum)) < 1e-6))
  v <- first_stage_values(f, r, "RRC_2026", trait = "fresh root yield")
  expect_equal(nrow(v), nrow(dd))
  expect_false(anyNA(v$`fresh root yield`))
})

test_that("first stage saved as a trait of the field plots; second stage as an analysis", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  fit <- analyze_augmented_row_column_design(d, "yield")
  v <- first_stage_values(fit, rc, "AYT_2026", trait = "fresh root yield")
  adjusted <- c("fresh root yield" = "fresh root yield spatially adjusted|COMP:0000123")

  # Without a connection the adjusted trait is "trait|statistic".
  obs0 <- suppressMessages(save_breedbase_first_stage(v, dry_run = TRUE))
  expect_equal(unique(obs0$trait), "fresh root yield|Adjusted Means from Spatial Correction using SpATS R")
  expect_error(save_breedbase_first_stage(v, as_trait = c("fresh root yield" = "fresh root yield|CO_334:0000013"),
                                          dry_run = TRUE), "mixing raw and adjusted")
  obs <- suppressMessages(save_breedbase_first_stage(v, as_trait = adjusted, dry_run = TRUE))
  expect_equal(nrow(obs), 144)
  expect_equal(unique(obs$trait), unname(adjusted))

  means <- second_stage_values(fit, trait = "fresh root yield")
  expect_equal(nrow(means), 112)
  sent <- suppressMessages(save_breedbase_analysis(means, analysis_name = "AYT_2026_means", program = "IITA",
                                                   description = "Row-column model", dry_run = TRUE))
  expect_equal(sent$analysis_result_values_type, "analysis_result_values_match_accession_names")
  expect_equal(sent$analysis_statistical_ontology_term, "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036")
  expect_length(sent$analysis_result_values, 112)
  expect_error(save_breedbase_analysis(means, analysis_name = "A", program = "IITA", dry_run = TRUE),
               "Give a `description`")

  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    out <- capture.output(res <- save_breedbase_first_stage(v, con, as_trait = adjusted, confirm = FALSE))
    expect_equal(res$values_found, 144)
    expect_equal(unique(vapply(fake$state$observations, `[[`, "", "observationVariableDbId")), "70003")

    out <- capture.output(res <- save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means",
                                                         program = "IITA", description = "Row-column model",
                                                         confirm = FALSE))
    expect_equal(res$accessions_found, 112)
    expect_match(res$analysis_page, "/analyses/[0-9]+$")
    a <- fake$state$analyses[[1]]
    expect_equal(a$form$analysis_trait_names, '["fresh root yield|CO_334:0000013"]')   # a JSON array
    expect_equal(a$form$analysis_breeding_program_id, "10")
    expect_equal(a$form$sgn_session_id, "TOKEN123")
    v1 <- a$values[["G1"]][["fresh root yield|CO_334:0000013"]]
    expect_equal(as.numeric(v1[[1]]), means$`fresh root yield`[means$accession == "G1"], tolerance = 1e-8)
    expect_equal(v1[[3]], "jdoe")

    # The same name again is refused before writing; unknown traits are listed.
    expect_error(capture.output(save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means",
                                                        program = "IITA", description = "x", confirm = FALSE)),
                 "already exists")
    bad <- means
    names(bad)[2] <- "root weight"
    expect_error(capture.output(save_breedbase_analysis(bad, con, analysis_name = "Other", program = "IITA",
                                                        description = "x", confirm = FALSE)),
                 "not found in Breedbase: 'root weight'")
  })
})

test_that("first stage: the adjusted trait is found or created (trait + statistic)", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  v <- first_stage_values(analyze_augmented_row_column_design(d, "yield"), rc, "AYT_2026", trait = "fresh root yield")
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    a <- breedbase_adjusted_trait(con, "fresh root yield")
    expect_equal(a$status, "missing")
    expect_equal(a$adjusted_trait, "fresh root yield|Adjusted Means from Spatial Correction using SpATS R")

    # A dry run reports the new trait and creates nothing.
    out <- capture.output(save_breedbase_first_stage(v, con, dry_run = TRUE))
    expect_true(any(grepl("\\[NEW\\]  trait 'fresh root yield\\|Adjusted Means.*will be created", out)))
    expect_length(fake$state$composed, 0)
    expect_error(capture.output(save_breedbase_first_stage(v, con, create_trait = FALSE, confirm = FALSE)),
                 "do not exist")

    # Saving creates the composed trait (trait id + statistic id) and stores the values under it.
    out <- capture.output(res <- save_breedbase_first_stage(v, con, confirm = FALSE))
    expect_length(fake$state$composed, 1)
    expect_setequal(fake$state$composed[[1]]$components, c("70001", "80054"))
    expect_equal(res$values_found, 144)
    expect_equal(unique(vapply(fake$state$observations, `[[`, "", "observationVariableDbId")),
                 fake$state$composed[[1]]$id)

    # The next time it is found, not created again.
    a2 <- breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)
    expect_equal(a2$status, "exists")
    expect_length(fake$state$composed, 1)
    # Another statistic gives another trait; an unknown statistic is a clear error.
    a3 <- breedbase_adjusted_trait(con, "fresh root yield", statistic = "Adjusted means from BLUEs using LMER R",
                                   create = TRUE)
    expect_equal(a3$status, "created")
    expect_error(breedbase_adjusted_trait(con, "fresh root yield", statistic = "Magic means"), "was not found")
  })
})
