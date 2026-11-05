# Extracted from test-breedbase.R:1220

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "breedtools", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
reset_found_cache <- function() rm(list = ls(breedtools:::bb_found_cache), envir = breedtools:::bb_found_cache)
fake_breedbase <- function(accessions = c(paste0("G", 1:200), paste0("C", 1:6), LETTERS[1:4]),
                           existing_studies = character(0), password = "secret", max_search = Inf,
                           slow_writes = character(0), lost_writes = character(0), site_auth = NULL,
                           statistics = NULL) {
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
                                cvterm_name = "Adjusted means from BLUEs using LMER R", accession = "SGNSTAT:0000036"),
                           list(cvterm_id = 80034, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted means from BLUPs using LMER R", accession = "SGNSTAT:0000034"),
                           list(cvterm_id = 80037, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Phenotypic BLUEs using LMER R", accession = "SGNSTAT:0000037"))
  if (!is.null(statistics)) state$statistics <- statistics
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
      # A test site with a web password (nginx basic auth) in front of the website.
      if (!is.null(site_auth)) {
        want <- paste("Basic", jsonlite::base64_enc(paste0(site_auth[1], ":", site_auth[2])))
        got <- httr2::req_get_headers(req, redacted = "reveal")$Authorization %||% ""
        if (!identical(got, want)) {
          return(httr2::response(401, headers = list(`WWW-Authenticate` = 'Basic realm="Restricted Test Site"',
                                                     `Content-Type` = "text/html"),
                                 body = charToRaw("<html>401 Authorization Required</html>")))
        }
      }
      form <- if (is.null(req$body)) {
        list()
      } else if (identical(req$body$type, "multipart")) {
        lapply(req$body$data, function(x) if (inherits(x, "form_file")) x else as.character(x))
      } else {
        lapply(req$body$data, function(x) utils::URLdecode(as.character(x)))
      }
      cookie <- httr2::req_get_headers(req, redacted = "reveal")$Cookie %||% ""
      session <- form$sgn_session_id %||% url$query$sgn_session_id %||% sub("^.*sgn_session_id=([^;]*).*$", "\\1", cookie)
      if (!identical(session, "TOKEN123")) return(httr2::response_json(200, body = list(error = "You must be logged in to do this!")))
      if (path == "ajax/breeder/check_status") {
        return(httr2::response_json(200, body = list(timestamp = "2026-11-04")))
      }
      if (path == "ajax/breeder/refresh") {
        state$refreshed <- c(state$refreshed, url$query$matviews)
        return(httr2::response_json(200, body = list(message = "refreshing")))
      }
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
      if (path %in% c("ajax/phenotype/upload_verify/spreadsheet", "ajax/phenotype/upload_store/spreadsheet")) {
        # Like Breedbase: units may be analysis plots ("<analysis>_<stock>"), traits must exist,
        # existing values are warned about unless overwritten.
        f <- utils::read.csv(form$upload_spreadsheet_phenotype_file_input$path, check.names = FALSE, colClasses = "character")
        units <- unlist(lapply(state$analyses, function(a) paste0(a$name, "_", names(a$values))))
        bad <- setdiff(f$observationunit_name, units)
        if (length(bad)) return(httr2::response_json(200, body = list(error = list(paste("Plots not in database:", paste(bad, collapse = ","))), success = list())))
        known <- c(vapply(state$variables, `[[`, "", "observationVariableName"),
                   vapply(state$composed, function(cp) cp$name, ""),
                   vapply(state$composed, function(cp) paste0(cp$name, "|COMP:000", cp$id), ""))
        badt <- setdiff(names(f)[-1], known)
        if (length(badt)) return(httr2::response_json(200, body = list(error = list(paste("Traits not in database:", paste(badt, collapse = ","))), success = list())))
        key <- unlist(lapply(names(f)[-1], function(tr) paste(f$observationunit_name[nzchar(f[[tr]])], tr)))
        existing <- intersect(key, names(state$analysis_values))
        warn <- if (length(existing) && is.null(form$phenotype_upload_overwrite_values)) as.list(paste(existing, "already has a <strong>different value</strong>")) else list()
        if (grepl("store", path)) {
          if (length(warn)) return(httr2::response_json(200, body = list(error = list("values exist"), success = list())))
          for (tr in names(f)[-1]) for (i in which(nzchar(f[[tr]]))) state$analysis_values[[paste(f$observationunit_name[i], tr)]] <- f[[tr]][i]
        }
        return(httr2::response_json(200, body = list(success = list("File data verified."), warning = warn, error = list())))
      }
      if (path == "ajax/analysis/store/spreadsheet") {
        # Like Breedbase: name unique; "stock_names" then traits; stocks are
        # plots or accessions depending on "Analysis of"; traits must exist.
        nm <- form$upload_new_analysis_name
        if (nm %in% c(vapply(state$analyses, `[[`, "", "name"), vapply(state$studies, `[[`, "", "studyName"))) {
          return(httr2::response_json(200, body = list(error = paste("An analysis with name", nm, "already exists"))))
        }
        f <- utils::read.csv(form$upload_new_analysis_file$path, check.names = FALSE, colClasses = "character")
        if (names(f)[1] != "stock_names") return(httr2::response_json(200, body = list(error_messages = list("Header row must contain stock_names"))))
        known <- c(vapply(state$variables, `[[`, "", "observationVariableName"),
                   vapply(state$composed, function(cp) paste0(cp$name, "|COMP:000", cp$id), ""))
        bad_traits <- setdiff(names(f)[-1], known)
        if (length(bad_traits)) return(httr2::response_json(200, body = list(error_messages = list(paste("The following traits are not in the database:", paste(bad_traits, collapse = ","))))))
        plots_level <- form$upload_new_analysis_result_values_type == "analysis_result_values_match_precomputed_design"
        valid <- if (plots_level) vapply(state$plots, `[[`, "", "observationUnitName") else accessions
        bad_stocks <- setdiff(f$stock_names, valid)
        if (length(bad_stocks)) return(httr2::response_json(200, body = list(error_messages = list(paste("The following stocks are not in the database:", paste(bad_stocks, collapse = ","))))))
        id <- as.character(state$next_id <- state$next_id + 1)
        state$analyses[[length(state$analyses) + 1]] <- list(id = id, name = nm, form = form, file = f,
                                                             values = stats::setNames(as.list(f$stock_names), f$stock_names))
        return(httr2::response_json(200, body = list(success = 1, analysis_id = as.integer(id), model_id = 5)))
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
      composed <- Filter(function(cp) cp$name %in% wanted, state$composed)
      hits <- c(hits, lapply(composed, function(cp) list(observationVariableDbId = cp$id,
                                                          observationVariableName = paste0(cp$name, "|COMP:000", cp$id))))
      return(ok(hits))
    }
    if (path == "search/observations") {
      unit_study <- vapply(state$plots, function(u) as.character(u$studyDbId), "")
      unit_name <- vapply(state$plots, `[[`, "", "observationUnitName")
      names(unit_study) <- names(unit_name) <- vapply(state$plots, `[[`, "", "observationUnitDbId")
      hits <- Filter(function(o) {
        (is.null(body$observationUnitDbIds) || o$observationUnitDbId %in% unlist(body$observationUnitDbIds)) &&
          (is.null(body$studyDbIds) || unit_study[[o$observationUnitDbId]] %in% unlist(body$studyDbIds)) &&
          o$observationVariableDbId %in% unlist(body$observationVariableDbIds)
      }, state$observations)
      hits <- lapply(hits, function(o) c(o, list(observationUnitName = unit_name[[o$observationUnitDbId]])))
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

# test -------------------------------------------------------------------------
skip_if_not_installed("httr2")
skip_if_not_installed("lme4")
sites <- c("Ibadan", "Ikenne")
designs <- lapply(1:2, function(i) augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                                                 n_candidates = 5, seed = i))
set.seed(1)
g <- stats::setNames(rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
fake <- fake_breedbase(accessions = c(paste0("G", 1:108), paste0("C", 1:4)))
httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    for (i in 1:2) {
      capture.output(save_breedbase(designs[[i]], con, program = "IITA", trial_name = paste0("T_", sites[i]),
                                    year = 2026, location = sites[i], confirm = FALSE, description = "x"))
      sheet <- phenotype_template(designs[[i]], paste0("T_", sites[i]), traits = "fresh root yield")
      sheet$`fresh root yield` <- round(50 + g[sheet$accession] + 0.5 * sheet$row + rnorm(nrow(sheet)), 2)
      sheet$`fresh root yield`[5] <- 500                         # an outlier
      capture.output(save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    }
    out <- capture.output(pheno <- breedbase_phenotypes(con, trials = paste0("T_", sites), traits = "fresh root yield"))
    expect_equal(nrow(pheno), 288)
    expect_equal(sum(!is.na(pheno$`fresh root yield`)), 288)
    expect_equal(attr(pheno, "trait_names")[["fresh root yield"]], "fresh root yield|CO_334:0000013")
    expect_true(all(c("plot_name", "accession", "rep", "block", "row", "col") %in% names(pheno)))
    expect_error(breedbase_phenotypes(con, "NOPE", "fresh root yield"), "Trial\\(s\\) not found")

    clean <- remove_outliers(pheno, "fresh root yield", by = "trial_name", quiet = TRUE)
    expect_true(500 %in% attr(clean, "outliers")$value)
    two <- analyze_two_stage(clean, "fresh root yield", stage1 = "blocks", genotype_col = "accession",
                             location_col = "trial_name")
    one <- analyze_single_stage(clean, "fresh root yield")
    expect_equal(nrow(two$adjusted_means), 112)
    expect_equal(nrow(one$adjusted_means), 112)
    expect_gt(cor(two$adjusted_means$adjusted_mean, g[two$adjusted_means$genotype]), 0.8)
    m <- merge(two$adjusted_means, one$adjusted_means, by = "genotype")
    expect_gt(cor(m$adjusted_mean.x, m$adjusted_mean.y), 0.95)      # two-stage and single-stage agree

    s1 <- first_stage_values(two, trait = "fresh root yield")
    expect_setequal(s1$plot_name, clean$plot_name[!is.na(clean$`fresh root yield`)])
    s2 <- second_stage_values(two, trait = "fresh root yield")
    ss <- second_stage_values(one, trait = "fresh root yield")

    capture.output(save_breedbase_analysis(s1, con, analysis_name = "MET_stage1", program = "IITA",
                                           description = "first stage", level = "plot", estimate = "blue",
                                           confirm = FALSE))
    capture.output(save_breedbase_analysis(s2, con, analysis_name = "MET_means", program = "IITA",
                                           description = "two-stage means", confirm = FALSE))
    # A second analysis in the same project: single-stage means under "trait|statistic".
    out <- capture.output(res <- add_breedbase_analysis_values(ss, con, analysis_name = "MET_means", level = "accession",
                                                               statistic = "Phenotypic BLUEs using LMER R",
                                                               confirm = FALSE))
    expect_equal(res$plots, 112)
    expect_true(all(grepl("^MET_means_.* fresh root yield\\|Phenotypic BLUEs using LMER R\\|COMP:", names(fake$state$analysis_values))))
  })
