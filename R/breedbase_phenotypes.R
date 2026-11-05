# Phenotypic data to Breedbase: a template from the design and the upload of
# observations through BrAPI v2 (POST /observations).
#
# Breedbase stores an observation from three ids: the plot
# (observationUnitDbId), the trait (observationVariableDbId) and the value.
# The ids are looked up from the plot and trait names given by the user.

#' Phenotyping template with the Breedbase plot names of a design
#'
#' Makes a table with one row per plot, the plot names that
#' [save_breedbase()] gives (or gave) to the plots, the design columns
#' (accession, block, rep, row, column) and one empty column per trait. Fill
#' it in (in R or in a spreadsheet) and upload it with
#' [save_breedbase_phenotypes()].
#'
#' @param design A design from [augmented_row_column_design()],
#'   [augmented_prep_design()], [resolvable_row_column_design()],
#'   [as_breedbase_design()] or [prepare_breedbase_trial()].
#' @param trial_name The trial name used when saving the design. For p-rep
#'   designs saved as one trial per location, the location is added as in
#'   [save_breedbase()]. Can be omitted for a prepared Breedbase trial.
#' @param traits Trait names as in Breedbase, e.g. `"fresh root yield"` or
#'   `"fresh root yield|CO_334:0000013"`. Can be empty (no trait columns).
#' @param plot_name_pattern,plot_names The same plot naming used in
#'   [save_breedbase()].
#'
#' @return A tibble: `plot_name`, `trial_name`, `accession`, `plot_number`,
#'   `block`, `rep`, `row`, `col`, `is_check` and one column per trait.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' sheet <- phenotype_template(rc, "AYT_2026_Ibadan",
#'                             traits = c("fresh root yield", "dry matter content percentage"))
#' head(sheet)
#'
#' @export
phenotype_template <- function(design, trial_name, traits = character(0),
                               plot_name_pattern = "{trial}_{accession}_{block}_{rep}", plot_names = NULL) {
  if (missing(trial_name) && inherits(design, "breedbase_trial")) trial_name <- design$trial$trial_name
  out <- bb_design_plot_table(design, trial_name, plot_name_pattern, plot_names)
  out$.design_location <- NULL
  out$.design_plots <- NULL
  for (t in unique(as.character(traits))) out[[t]] <- NA_real_
  out
}

# Breedbase plot names of a design with the design's own plot number (and
# location for p-rep designs), to match analysis results.
bb_design_plot_table <- function(design, trial_name, plot_name_pattern = "{trial}_{accession}_{block}_{rep}",
                                 plot_names = NULL) {
  if (is.null(trial_name) && inherits(design, "breedbase_trial")) trial_name <- design$trial$trial_name
  studies <- bb_design_to_studies(design, trial_name, year = 2000, location = bb_template_location(design),
                                  description = "template", trial_type = "phenotyping_trial", design_type = NULL,
                                  plot_names = plot_names, plot_name_pattern = plot_name_pattern)
  rows <- lapply(studies, function(s) {
    level <- function(p, name) {
      rel <- Filter(function(r) identical(r$levelName, name), p$observationUnitPosition$observationLevelRelationships)
      if (length(rel)) as.integer(rel[[1]]$levelCode) else NA_integer_
    }
    tibble::tibble(
      plot_name = vapply(s$plots, `[[`, character(1), "observationUnitName"),
      trial_name = s$study$studyName,
      accession = vapply(s$plots, `[[`, character(1), "germplasmName"),
      plot_number = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$observationLevel$levelCode), integer(1)),
      block = vapply(s$plots, level, integer(1), name = "block"),
      rep = vapply(s$plots, level, integer(1), name = "rep"),
      row = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$positionCoordinateY %||% NA), integer(1)),
      col = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$positionCoordinateX %||% NA), integer(1)),
      is_check = vapply(s$plots, function(p) identical(p$observationUnitPosition$entryType, "CHECK"), logical(1))
    )
  })
  out <- dplyr::bind_rows(rows)
  # The design's own plot number (and location) of each row, in the same
  # order save_breedbase() uses: per location, fillers left out, by plot.
  d <- as.data.frame(design$design, stringsAsFactors = FALSE)
  if ("is_filler" %in% names(d)) d <- d[d$is_filler != 1, , drop = FALSE]
  groups <- if (inherits(design, "augmented_prep_design")) {
    lapply(design$input$locations, function(l) d[d$location == l, , drop = FALSE])
  } else {
    list(d)
  }
  ordered <- dplyr::bind_rows(lapply(groups, function(x) x[order(x$plots), , drop = FALSE]))
  out$.design_location <- if ("location" %in% names(ordered)) as.character(ordered$location) else NA_character_
  out$.design_plots <- as.integer(ordered$plots)
  out
}

# The template only needs the location names of p-rep designs (to name the
# trials per location); other designs need a placeholder.
bb_template_location <- function(design) {
  if (inherits(design, "augmented_prep_design")) NULL else "template"
}

#' Save phenotypic data to Breedbase
#'
#' Uploads trait values of plots that are already in Breedbase (for example
#' saved with [save_breedbase()]) through BrAPI. Give the plot name and the
#' trait name with each value; the plots and traits are looked up in
#' Breedbase first, and every problem is listed before anything is written.
#'
#' The data can be wide, one row per plot and one column per trait (as made
#' by [phenotype_template()]):
#'
#' | plot_name                 | fresh root yield | dry matter content percentage |
#' |---------------------------|------------------|-------------------------------|
#' | AYT_2026_Ibadan_G31_1_1   | 25.4             | 31.2                          |
#'
#' or long, one row per value, with `trait` and `value` columns
#' (`trait = "trait", value = "value"`).
#'
#' Trait names are the Breedbase names, with or without the ontology id:
#' `"fresh root yield"` or `"fresh root yield|CO_334:0000013"` (the column
#' names of a Breedbase phenotype download). Empty values (`NA` or `""`) are
#' skipped. Values are checked by Breedbase against the trait (e.g. numbers
#' within the allowed range).
#'
#' @param data A data frame with the plot names and the values.
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns the observations without checking them.
#' @param plot Name of the column with the plot names.
#' @param traits For wide data: the trait columns. By default every column
#'   except `plot`, `timestamp` and the design columns of
#'   [phenotype_template()] (`trial_name`, `accession`, `plot_number`,
#'   `block`, `rep`, `row`, `col`, `is_check`).
#' @param trait,value For long data: the names of the trait and value columns.
#' @param timestamp Optional name of a column with the date and time of each
#'   measurement (e.g. `"2026-08-15 10:30:00"`). Breedbase uses it to tell
#'   repeated measurements apart.
#' @param dry_run Logical. Check everything but write nothing.
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#' @param refresh Logical. After saving, refresh the phenotype cache of
#'   Breedbase so the values show on the trial pages (see
#'   [breedbase_refresh()]).
#'
#' @inheritSection breedbase_check Slow servers
#'
#' @return Invisibly, a tibble with one row per trait: `trait`,
#'   `observation_variable_db_id`, `values_sent` and `values_found` (read back
#'   from Breedbase after saving).
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' sheet <- phenotype_template(rc, "AYT_2026_Ibadan", traits = "fresh root yield")
#' sheet$`fresh root yield` <- round(rnorm(nrow(sheet), 25, 4), 1)
#'
#' # What would be sent (no connection needed):
#' obs <- save_breedbase_phenotypes(sheet, dry_run = TRUE)
#' head(obs)
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase_phenotypes(sheet, con, dry_run = TRUE)   # check first
#' save_breedbase_phenotypes(sheet, con)
#'
#' # Long data: one row per value
#' long <- data.frame(plot_name = sheet$plot_name, trait = "fresh root yield",
#'                    value = sheet$`fresh root yield`)
#' save_breedbase_phenotypes(long, con, trait = "trait", value = "value")
#' }
#'
#' @export
save_breedbase_phenotypes <- function(data, con = NULL, plot = "plot_name", traits = NULL,
                                      trait = NULL, value = NULL, timestamp = NULL,
                                      dry_run = FALSE, confirm = interactive(), refresh = TRUE) {
  bb_save_phenotypes(data, con, plot, traits, trait, value, timestamp, dry_run, confirm, refresh = refresh)
}

# `trait_ids`: trait ids already known, named by trait (e.g. composed traits,
# which the BrAPI variable search may not list).
bb_save_phenotypes <- function(data, con, plot = "plot_name", traits = NULL, trait = NULL, value = NULL,
                               timestamp = NULL, dry_run = FALSE, confirm = interactive(), trait_ids = NULL,
                               refresh = TRUE) {
  obs <- bb_phenotypes_long(data, plot, traits, trait, value, timestamp)

  if (is.null(con)) {
    if (!isTRUE(dry_run)) stop("`con` is needed to save; use breedbase_login() first, or `dry_run = TRUE`.", call. = FALSE)
    message("Dry run without a connection: returning the observations that would be sent.")
    return(obs)
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)

  cat("Checking phenotypes before saving to ", con$url, "\n", sep = "")
  cat("  [OK]   ", nrow(obs), " values of ", length(unique(obs$trait)), " trait",
      if (length(unique(obs$trait)) > 1) "s" else "", " on ", length(unique(obs$plot_name)), " plots\n", sep = "")

  # Plots and traits are both checked before stopping, so every problem is
  # listed at once.
  plot_names <- unique(obs$plot_name)
  plot_ids <- bb_search_names(con, "observationunits", "observationUnitNames", "observationUnitName", plot_names,
                              batch = getOption("breedtools.plot_batch", 200), label = "plots",
                              id_field = "observationUnitDbId")
  known_ids <- trait_ids
  trait_ids <- bb_find_variables(con, setdiff(unique(obs$trait), names(known_ids)))
  trait_ids <- c(trait_ids, known_ids)[unique(obs$trait)]
  problems <- character(0)
  missing_plots <- setdiff(plot_names, names(plot_ids))
  if (length(missing_plots)) {
    problems <- c(problems, paste0(
      length(missing_plots), " of ", length(plot_names), " plots are not in Breedbase: ",
      paste(utils::head(missing_plots, 10), collapse = ", "), if (length(missing_plots) > 10) ", ..." else "",
      ". Save the trial first (save_breedbase()) and use its plot names (phenotype_template())."))
  } else {
    cat("  [OK]   all ", length(plot_names), " plots are in Breedbase\n", sep = "")
  }
  missing_traits <- names(trait_ids)[is.na(trait_ids)]
  if (length(missing_traits)) {
    problems <- c(problems, paste0(
      "Trait", if (length(missing_traits) > 1) "s" else "", " not found in Breedbase: ",
      paste0("'", missing_traits, "'", collapse = ", "),
      ". Use the trait name as shown in Breedbase, e.g. \"fresh root yield\" or ",
      "\"fresh root yield|CO_334:0000013\"."))
  } else {
    for (t in names(trait_ids)) cat("  [OK]   trait '", t, "' (id ", trait_ids[[t]], ")\n", sep = "")
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)

  obs$observation_unit_db_id <- unname(plot_ids[obs$plot_name])
  obs$observation_variable_db_id <- unname(trait_ids[obs$trait])
  records <- lapply(seq_len(nrow(obs)), function(i) {
    r <- list(observationUnitDbId = obs$observation_unit_db_id[i],
              observationVariableDbId = obs$observation_variable_db_id[i],
              value = obs$value[i])
    if (!is.na(obs$timestamp[i])) r$observationTimeStamp <- obs$timestamp[i]
    r
  })

  plan <- sprintf("%d values of %d trait%s on %d plots at %s", nrow(obs), length(trait_ids),
                  if (length(trait_ids) > 1) "s" else "", length(plot_names), con$url)
  if (isTRUE(dry_run)) {
    cat("Dry run: everything is ready to save ", plan, ". Nothing was written.\n", sep = "")
    return(invisible(bb_phenotype_summary(obs, found = NULL)))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }

  batch <- getOption("breedtools.observation_batch", 500)
  groups <- split(seq_along(records), ceiling(seq_along(records) / batch))
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    bb_write_checked(
      what = sprintf("save values %d-%d of %d", min(idx), max(idx), length(records)),
      write = function() {
        bb_post(con, "observations", records[idx], what = "save the phenotypes")
        "done"
      },
      find = function() {
        n <- bb_count_observations(con, obs[idx, , drop = FALSE])
        if (n >= length(idx)) "done" else NA_character_
      }
    )
    if (length(groups) > 1) {
      cat(sprintf("  ...    batch %d of %d saved: %d/%d values\n", g, length(groups), max(idx), length(records)))
    }
  }

  found <- bb_observations_found(con, obs)
  out <- bb_phenotype_summary(obs, found)
  for (i in seq_len(nrow(out))) {
    status <- if (out$values_found[i] >= out$values_sent[i]) "[OK]  " else "[FAIL]"
    cat(sprintf("  %s trait '%s': %d of %d values saved\n", status, out$trait[i], out$values_found[i],
                out$values_sent[i]))
  }
  # Breedbase's BrAPI upload does not refresh the phenotype cache that the
  # trial pages and the search wizard read; refresh it so the values show.
  if (isTRUE(refresh)) breedbase_refresh(con, quiet = FALSE)
  invisible(out)
}

#' Refresh the Breedbase phenotype cache
#'
#' Values saved through BrAPI ([save_breedbase_phenotypes()],
#' [save_breedbase_first_stage()]) are in the database at once (the BrAPI
#' search finds them), but the trial pages and the search wizard of
#' Breedbase read a cache ("materialized views") that the BrAPI upload does
#' not refresh. This asks Breedbase to refresh it, like the "Update" button
#' of the search wizard. The save functions call it after saving. If it
#' cannot be reached, the values are still saved: they appear after the next
#' refresh of the site (usually overnight), or when you click "Update" on the
#' Search Wizard page.
#'
#' @param con A connection from [breedbase_login()].
#' @param views Which cache to refresh: `"phenotypes"` (default) or
#'   `"fullview"`.
#' @param quiet Logical. If `FALSE`, print what happened.
#'
#' @return Invisibly, `TRUE` if Breedbase accepted the refresh, else `FALSE`.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_refresh(con)
#' }
#'
#' @export
breedbase_refresh <- function(con, views = c("phenotypes", "fullview"), quiet = FALSE) {
  views <- match.arg(views)
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  res <- tryCatch(bb_ajax_get(con, "ajax/breeder/refresh", what = "refresh the phenotype cache",
                              query = list(matviews = views)),
                  breedbase_error = function(e) e)
  ok <- !inherits(res, "breedbase_error")
  if (!isTRUE(quiet)) {
    if (ok) {
      cat("  [OK]   Breedbase is refreshing its phenotype cache: the values appear on the trial pages",
          "and in the search wizard in a few minutes\n")
    } else {
      cat("  [NOTE] The values are saved, but the phenotype cache could not be refreshed (",
          conditionMessage(res), ").\n         They appear on the trial pages after the next refresh of the site, ",
          "or click \"Update\" on the Search Wizard page.\n", sep = "")
    }
  }
  invisible(ok)
}

# Data (wide or long) to one row per value: plot_name, trait, value, timestamp.
bb_phenotypes_long <- function(data, plot, traits, trait, value, timestamp) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c(plot, trait, value, timestamp)
  absent <- setdiff(need, names(data))
  if (length(absent)) {
    stop("Column", if (length(absent) > 1) "s" else "", " not found in `data`: ",
         paste0("'", absent, "'", collapse = ", "), ". Columns are: ", paste(names(data), collapse = ", "), ".",
         call. = FALSE)
  }
  ts <- if (is.null(timestamp)) rep(NA_character_, nrow(data)) else bb_format_time(data[[timestamp]])
  if (!is.null(trait) || !is.null(value)) {
    if (is.null(trait) || is.null(value)) stop("Long data needs both `trait` and `value` columns.", call. = FALSE)
    long <- data.frame(plot_name = as.character(data[[plot]]), trait = as.character(data[[trait]]),
                       value = bb_value_text(data[[value]]), timestamp = ts, stringsAsFactors = FALSE)
  } else {
    design_cols <- c("trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
    if (is.null(traits)) traits <- setdiff(names(data), c(plot, timestamp, design_cols))
    absent <- setdiff(traits, names(data))
    if (length(absent)) stop("Trait columns not found in `data`: ", paste(absent, collapse = ", "), call. = FALSE)
    if (!length(traits)) stop("No trait columns in `data`. Give them with `traits = c(...)`.", call. = FALSE)
    long <- do.call(rbind, lapply(traits, function(t) {
      data.frame(plot_name = as.character(data[[plot]]), trait = t, value = bb_value_text(data[[t]]),
                 timestamp = ts, stringsAsFactors = FALSE)
    }))
  }
  long$plot_name <- trimws(long$plot_name)
  long$trait <- trimws(long$trait)
  long <- long[!is.na(long$value) & nzchar(long$value), , drop = FALSE]
  if (!nrow(long)) stop("There are no values to save (all are empty).", call. = FALSE)
  if (anyNA(long$plot_name) || any(!nzchar(long$plot_name))) stop("Some rows have no plot name.", call. = FALSE)
  key <- paste(long$plot_name, long$trait, long$timestamp, sep = "\r")
  if (anyDuplicated(key)) {
    d <- long[duplicated(key), , drop = FALSE][1, ]
    stop("Plot '", d$plot_name, "' has more than one value of '", d$trait, "'",
         if (is.na(d$timestamp)) "; give a `timestamp` column to save repeated measurements" else " at the same time",
         ".", call. = FALSE)
  }
  rownames(long) <- NULL
  tibble::as_tibble(long)
}

bb_value_text <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  trimws(out)
}

bb_format_time <- function(x) {
  if (inherits(x, "POSIXt")) return(format(x, "%Y-%m-%d %H:%M:%S"))
  if (inherits(x, "Date")) return(format(x, "%Y-%m-%d 00:00:00"))
  out <- trimws(as.character(x))
  out[out == ""] <- NA_character_
  out
}

# Trait ids by name. Breedbase searches the plain trait name and shows it as
# "name|ONTOLOGY:id"; both forms are accepted.
bb_find_variables <- function(con, traits) {
  if (!length(traits)) {
    out <- stats::setNames(character(0), character(0))
    attr(out, "full_name") <- out
    return(out)
  }
  base <- trimws(sub("\\|.*$", "", traits))
  data <- bb_search(con, "variables", list(observationVariableNames = as.list(unique(base)),
                                           pageSize = 10 * length(unique(base))))
  full <- vapply(data, function(x) as.character(x$observationVariableName %||% NA_character_), character(1))
  ids <- vapply(data, function(x) as.character(x$observationVariableDbId %||% NA_character_), character(1))
  short <- trimws(sub("\\|.*$", "", full))
  hits <- lapply(seq_along(traits), function(i) {
    if (grepl("|", traits[i], fixed = TRUE)) which(full == traits[i]) else which(short == base[i])
  })
  pick <- function(v) vapply(hits, function(h) if (length(h) == 1) v[h] else NA_character_, character(1))
  out <- stats::setNames(pick(ids), traits)
  # The full name "trait|ONTOLOGY:id" (what Breedbase analyses need).
  attr(out, "full_name") <- stats::setNames(pick(full), traits)
  out
}

# Stored observations of these plots and traits, searched a few plots at a
# time so each answer stays small.
bb_search_observations <- function(con, obs) {
  units <- unique(obs$observation_unit_db_id)
  vars <- unique(obs$observation_variable_db_id)
  per_search <- max(1L, floor(getOption("breedtools.observation_batch", 500) / length(vars)))
  out <- list()
  for (chunk in split(units, ceiling(seq_along(units) / per_search))) {
    out <- c(out, bb_search(con, "observations", list(
      observationUnitDbIds = as.list(chunk),
      observationVariableDbIds = as.list(vars),
      pageSize = 2 * length(chunk) * length(vars)
    )))
  }
  out
}

# How many of these values are stored (same plot, trait and value).
bb_observations_found <- function(con, obs) {
  data <- bb_search_observations(con, obs)
  stored <- vapply(data, function(x) paste(x$observationUnitDbId, x$observationVariableDbId, x$value, sep = "\r"),
                   character(1))
  paste(obs$observation_unit_db_id, obs$observation_variable_db_id, obs$value, sep = "\r") %in% stored
}

bb_count_observations <- function(con, obs) sum(bb_observations_found(con, obs))

bb_phenotype_summary <- function(obs, found) {
  traits <- unique(obs$trait)
  tibble::tibble(
    trait = traits,
    observation_variable_db_id = vapply(traits, function(t) obs$observation_variable_db_id[obs$trait == t][1],
                                        character(1), USE.NAMES = FALSE),
    values_sent = vapply(traits, function(t) sum(obs$trait == t), integer(1), USE.NAMES = FALSE),
    values_found = if (is.null(found)) NA_integer_ else
      vapply(traits, function(t) sum(found[obs$trait == t]), integer(1), USE.NAMES = FALSE)
  )
}
