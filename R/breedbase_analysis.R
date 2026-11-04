# Analysis results to Breedbase.
#
# First stage (one value per field plot, e.g. spatially corrected plot
# values): saved as a trait on the plots of the field trial, through BrAPI
# observations (save_breedbase_phenotypes()).
#
# Second stage (one value per genotype, e.g. adjusted means across
# locations): saved as a Breedbase analysis (an "analysis" project), like the
# Breedbase mixed-models tool does, through POST /ajax/analysis/store/json
# (not BrAPI; it accepts the BrAPI login token as `sgn_session_id`). Each
# trait is saved composed with a statistic term, e.g. "fresh root
# yield|CO_334:0000013" + "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036".

#' First-stage plot values from an analysis
#'
#' Turns analyses of a trial into one value per field plot, ready for
#' [save_breedbase_first_stage()]: the genotype's adjusted mean plus the plot's
#' residual. This removes the design effects (blocks, rows, columns,
#' replicates) and the spatial trend, and keeps the plot-to-plot variation, so
#' the result has as many plots as the field trial. With SpATS
#' (`spatial = TRUE` or `stage1 = "spatial"`) these are spatially corrected
#' plot values.
#'
#' @param fits One analysis, or a named list of analyses, one per trait (the
#'   names are the Breedbase trait names). Analyses from
#'   [analyze_augmented_row_column_design()], [analyze_augmented_prep_design()],
#'   [analyze_augmented_prep_two_stage()] (stage-1 fits per location) or
#'   [analyze_resolvable_row_column()].
#' @param design The design that was saved to Breedbase.
#' @param trial_name The trial name used when saving the design.
#' @param trait Trait name when `fits` is a single analysis.
#' @param plot_name_pattern,plot_names The plot naming used in [save_breedbase()].
#'
#' @return A tibble with `plot_name`, `trial_name`, `accession`,
#'   `plot_number`, `block`, `rep`, `row`, `col`, `is_check` and one column
#'   per trait. Its attribute `"method"` is `"spatial"` when any analysis used
#'   SpATS, otherwise `"mixed"`.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield")
#' head(values)
#'
#' @export
first_stage_values <- function(fits, design, trial_name, trait = NULL,
                               plot_name_pattern = "{trial}_{accession}_{block}_{rep}", plot_names = NULL) {
  if (bb_is_analysis(fits)) {
    if (is.null(trait) || length(trait) != 1 || !nzchar(trait)) {
      stop("Give the Breedbase `trait` name of this analysis, e.g. trait = \"fresh root yield\".", call. = FALSE)
    }
    fits <- stats::setNames(list(fits), trait)
  }
  if (!is.list(fits) || is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("`fits` must be one analysis (with `trait`) or a list of analyses named by trait.", call. = FALSE)
  }
  for (nm in names(fits)) {
    if (!bb_is_analysis(fits[[nm]])) stop("`fits$", nm, "` is not an analysis from breedtools.", call. = FALSE)
  }

  table <- bb_design_plot_table(design, trial_name, plot_name_pattern, plot_names)
  is_prep <- inherits(design, "augmented_prep_design")
  key_table <- paste(if (is_prep) table$.design_location else "", table$.design_plots, sep = "\r")
  spatial <- FALSE
  for (nm in names(fits)) {
    plots <- bb_first_stage_plots(fits[[nm]])
    spatial <- spatial || isTRUE(attr(plots, "spatial"))
    if (!"plots" %in% names(plots)) stop("The analysis of '", nm, "' has no `plots` column to match.", call. = FALSE)
    key <- paste(if (is_prep) as.character(plots$location) else "", plots$plots, sep = "\r")
    unknown <- setdiff(key, key_table)
    if (length(unknown)) {
      stop("The analysis of '", nm, "' has ", length(unknown), " plots that are not in the design. ",
           "Was it run on this design?", call. = FALSE)
    }
    table[[nm]] <- plots$value[match(key_table, key)]
  }
  table$.design_location <- NULL
  table$.design_plots <- NULL
  attr(table, "method") <- if (spatial) "spatial" else "mixed"
  table
}

bb_is_analysis <- function(x) {
  inherits(x, c("augmented_row_column_analysis", "augmented_prep_analysis", "augmented_prep_two_stage",
                "resolvable_row_column_analysis")) ||
    (is.list(x) && !is.null(x$plot_results) && !is.null(x$adjusted_means))
}

# Plot values (adjusted mean of the genotype + plot residual) of one analysis.
bb_first_stage_plots <- function(fit) {
  if (inherits(fit, "augmented_prep_two_stage") || (!is.null(fit$stage1_fits) && is.null(fit$plot_results))) {
    parts <- lapply(fit$stage1_fits, bb_first_stage_plots)
    out <- dplyr::bind_rows(parts)
    attr(out, "spatial") <- any(vapply(parts, function(p) isTRUE(attr(p, "spatial")), logical(1))) ||
      ("stage1" %in% names(fit$diagnostics) && identical(as.character(fit$diagnostics$stage1[1]), "spatial"))
    return(out)
  }
  pr <- fit$plot_results
  if (is.null(pr) || !all(c("residual") %in% names(pr))) {
    stop("This analysis has no plot results (residuals).", call. = FALSE)
  }
  gcol <- fit$input$genotype_col %||% "all_entries"
  if (!gcol %in% names(pr)) stop("The plot results have no genotype column '", gcol, "'.", call. = FALSE)
  am <- fit$adjusted_means
  means <- stats::setNames(am$adjusted_mean, as.character(am$genotype))
  out <- tibble::tibble(
    plots = as.integer(pr$plots),
    location = if ("location" %in% names(pr)) as.character(pr$location) else NA_character_,
    value = unname(means[as.character(pr[[gcol]])]) + pr$residual
  )
  diag_value <- function(name) {
    dg <- fit$diagnostics
    if (is.data.frame(dg) && name %in% names(dg)) as.character(dg[[name]][1]) else ""
  }
  attr(out, "spatial") <- isTRUE(fit$input$spatial) ||
    grepl("spats|spatial", paste(diag_value("model_type"), diag_value("stage1_model")), ignore.case = TRUE)
  out
}

#' Find or create the Breedbase trait for adjusted values
#'
#' First-stage values need their own trait, so raw and adjusted values stay
#' apart. Breedbase makes such traits by composing the measured trait with a
#' statistic term, e.g. "fresh root yield" + "Adjusted Means from Spatial
#' Correction using SpATS R" gives the composed trait
#' `"fresh root yield|Adjusted Means from Spatial Correction using SpATS R"`
#' (ontology `COMP`). This function finds that trait, or creates it with
#' `create = TRUE` through the same calls as the "Compose a new trait" page
#' of Breedbase.
#'
#' @param con A connection from [breedbase_login()].
#' @param traits Measured trait names, as in Breedbase.
#' @param statistic The statistic term (from the SGNSTAT ontology of
#'   Breedbase), with or without its id. The default is the term of the
#'   Breedbase spatial correction tool; for non-spatial models
#'   `"Adjusted means from BLUEs using LMER R"` is another choice.
#' @param create Logical. Create the traits that do not exist yet.
#'
#' @return A tibble with `trait`, `adjusted_trait`, `adjusted_trait_db_id`
#'   and `status` (`"exists"`, `"created"` or `"missing"`).
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_adjusted_trait(con, "fresh root yield")                  # look only
#' breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)   # create if needed
#' }
#'
#' @export
breedbase_adjusted_trait <- function(con, traits,
                                     statistic = "Adjusted Means from Spatial Correction using SpATS R",
                                     create = FALSE) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  ids <- bb_find_variables(con, traits)
  full <- attr(ids, "full_name")
  if (anyNA(ids)) {
    stop("Trait(s) not found in Breedbase: ", paste0("'", traits[is.na(ids)], "'", collapse = ", "),
         ". Use the trait name as shown in Breedbase.", call. = FALSE)
  }
  stat <- bb_find_statistic(con, statistic)
  base <- sub("\\|[^|]*$", "", full)          # trait name without its ontology id
  out <- tibble::tibble(trait = traits, adjusted_trait = paste(base, stat$name, sep = "|"),
                        adjusted_trait_db_id = NA_character_, status = "missing")
  for (i in seq_along(traits)) {
    components <- c(ids[[i]], stat$id)
    id <- bb_composed_trait_id(con, components)
    if (is.na(id) && isTRUE(create)) {
      bb_ajax_post(con, "ajax/onto/store_composed_term",
                   list(new_trait_names = stats::setNames(list(paste(components, collapse = ",")),
                                                          out$adjusted_trait[i])),
                   what = paste0("create trait '", out$adjusted_trait[i], "'"), json = "new_trait_names")
      id <- bb_composed_trait_id(con, components)
      if (is.na(id)) stop("Breedbase did not create the trait '", out$adjusted_trait[i], "'.", call. = FALSE)
      out$status[i] <- "created"
    } else if (!is.na(id)) {
      out$status[i] <- "exists"
    }
    out$adjusted_trait_db_id[i] <- id
  }
  out
}

# The SGNSTAT term id of a statistic (name with or without "|SGNSTAT:id").
bb_find_statistic <- function(con, statistic) {
  name <- trimws(sub("\\|[^|]*$", "", statistic))
  db <- if (grepl("|", statistic, fixed = TRUE)) sub(":.*$", "", sub("^.*\\|", "", statistic)) else "SGNSTAT"
  hits <- bb_ajax_get(con, "ajax/onto/match", what = "find the statistic term",
                      query = list(db_name = db, term_name = name))
  rows <- if (is.null(names(hits))) hits else list()
  match <- Filter(function(h) identical(tolower(trimws(h$cvterm_name %||% "")), tolower(name)), rows)
  if (!length(match)) {
    stop("The statistic '", name, "' was not found in this Breedbase (ontology ", db, ").", call. = FALSE)
  }
  list(id = as.character(match[[1]]$cvterm_id), name = match[[1]]$cvterm_name)
}

# Id of the composed trait made of exactly these components, or NA.
bb_composed_trait_id <- function(con, components) {
  body <- bb_ajax_get(con, "ajax/onto/get_trait_from_exact_components", what = "look up the composed trait",
                      query = list(`ids[]` = components), allow_error = TRUE)
  if (is.null(body$trait_id)) NA_character_ else as.character(body$trait_id)
}

#' Save first-stage results to Breedbase as a trait of the field plots
#'
#' Saves plot-level first-stage results (from [first_stage_values()]) on the
#' plots of the field trial, as values of a trait made for adjusted values,
#' so raw and adjusted values stay apart. By default that trait is the
#' measured trait composed with a statistic, e.g. "fresh root yield|Adjusted
#' Means from Spatial Correction using SpATS R": it is found, or created
#' (after you confirm) when it does not exist yet; see
#' [breedbase_adjusted_trait()]. Or name your own traits with `as_trait`.
#' The values are saved through BrAPI, like [save_breedbase_phenotypes()].
#'
#' @param values A table from [first_stage_values()] (or any table with a
#'   `plot_name` column and one numeric column per trait).
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns what would be sent.
#' @param as_trait Optional named vector: for each trait column of `values`,
#'   an existing Breedbase trait to save it under, e.g.
#'   `c("fresh root yield" = "fresh root yield adjusted|COMP:0000123")`.
#' @param statistic The statistic composed with each trait when `as_trait`
#'   is not given (see [breedbase_adjusted_trait()]).
#' @param create_trait Logical. Create the composed traits that do not exist
#'   yet (you are asked first when `confirm = TRUE`).
#' @param dry_run,confirm As in [save_breedbase_phenotypes()].
#'
#' @inheritSection breedbase_check Slow servers
#'
#' @return As [save_breedbase_phenotypes()]: one row per trait with the
#'   values sent and found.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield")
#' obs <- save_breedbase_first_stage(values, dry_run = TRUE)
#' head(obs)
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase_first_stage(values, con, dry_run = TRUE)   # check; shows the trait to create
#' save_breedbase_first_stage(values, con)                   # creates the trait if needed, saves
#' }
#'
#' @export
save_breedbase_first_stage <- function(values, con = NULL, as_trait = NULL,
                                       statistic = "Adjusted Means from Spatial Correction using SpATS R",
                                       create_trait = TRUE, dry_run = FALSE, confirm = interactive()) {
  if (!is.data.frame(values) || !"plot_name" %in% names(values)) {
    stop("`values` must be a table with a `plot_name` column, as from first_stage_values().", call. = FALSE)
  }
  design_cols <- c("plot_name", "trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
  traits <- setdiff(names(values), design_cols)
  if (!length(traits)) stop("No trait columns in `values`.", call. = FALSE)
  values <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  stat_name <- trimws(sub("\\|[^|]*$", "", statistic))
  trait_ids <- NULL

  if (!is.null(as_trait)) {
    if (is.null(names(as_trait)) || any(!nzchar(names(as_trait)))) {
      stop("`as_trait` must be named by the trait columns of `values`.", call. = FALSE)
    }
    absent <- setdiff(names(as_trait), traits)
    if (length(absent)) {
      stop("Not trait columns of `values`: ", paste(absent, collapse = ", "), ". Trait columns are: ",
           paste(traits, collapse = ", "), ".", call. = FALSE)
    }
    same <- names(as_trait)[trimws(sub("\\|.*$", "", as_trait)) == trimws(sub("\\|.*$", "", names(as_trait)))]
    if (length(same)) {
      stop("First-stage values of '", same[1], "' would be saved as the measured trait itself, mixing raw and ",
           "adjusted values on the same plots. Choose a trait made for adjusted values.", call. = FALSE)
    }
    targets <- as_trait[traits[traits %in% names(as_trait)]]
  } else if (is.null(con)) {
    targets <- stats::setNames(paste(sub("\\|.*$", "", traits), stat_name, sep = "|"), traits)
  } else {
    if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
    cat("Finding the traits for the adjusted values in ", con$url, "\n", sep = "")
    adj <- breedbase_adjusted_trait(con, traits, statistic = statistic, create = FALSE)
    for (i in seq_len(nrow(adj))) {
      if (adj$status[i] == "exists") {
        cat("  [OK]   trait '", adj$adjusted_trait[i], "' (id ", adj$adjusted_trait_db_id[i], ")\n", sep = "")
      } else {
        cat("  [NEW]  trait '", adj$adjusted_trait[i], "' will be created\n", sep = "")
      }
    }
    missing <- adj$status == "missing"
    if (any(missing) && !isTRUE(dry_run)) {
      if (!isTRUE(create_trait)) {
        stop("The trait(s) for adjusted values do not exist: ",
             paste0("'", adj$adjusted_trait[missing], "'", collapse = ", "),
             ". Use create_trait = TRUE, or name existing traits with `as_trait`.", call. = FALSE)
      }
      if (isTRUE(confirm)) {
        answer <- readline(sprintf("Create %d new trait%s in Breedbase (%s)? Type 'yes' to continue: ",
                                   sum(missing), if (sum(missing) > 1) "s" else "",
                                   paste(adj$adjusted_trait[missing], collapse = ", ")))
        if (!identical(tolower(trimws(answer)), "yes")) {
          message("Nothing was created or saved.")
          return(invisible(NULL))
        }
      }
      adj <- breedbase_adjusted_trait(con, traits, statistic = statistic, create = TRUE)
      for (i in which(adj$status == "created")) {
        cat("  [OK]   trait '", adj$adjusted_trait[i], "' created (id ", adj$adjusted_trait_db_id[i], ")\n", sep = "")
      }
    }
    targets <- stats::setNames(adj$adjusted_trait, traits)
    trait_ids <- stats::setNames(ifelse(is.na(adj$adjusted_trait_db_id), "(new)", adj$adjusted_trait_db_id),
                                 adj$adjusted_trait)
  }

  out <- values[, c("plot_name", names(targets)), drop = FALSE]
  names(out) <- c("plot_name", unname(targets))
  bb_save_phenotypes(out, con, plot = "plot_name", traits = unname(targets), dry_run = dry_run,
                     confirm = confirm, trait_ids = trait_ids)
}

#' Second-stage genotype values from analyses
#'
#' Collects the adjusted means of analyses (one per trait) into one table
#' with one row per accession, ready for [save_breedbase_analysis()].
#'
#' @param fits One analysis, or a named list of analyses, one per trait (the
#'   names are the Breedbase trait names). Any breedtools analysis with
#'   `adjusted_means`, e.g. [analyze_augmented_prep_two_stage()].
#' @param trait Trait name when `fits` is a single analysis.
#'
#' @return A tibble with `accession` and one column per trait. Its attribute
#'   `"estimate"` is `"blup"` when the genotypes were random, else `"blue"`.
#'
#' @examples
#' prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
#' d <- prep$design
#' set.seed(1)
#' d$yield <- rnorm(nrow(d), 25, 3)
#' fit <- analyze_augmented_prep_two_stage(d, "yield")
#' head(second_stage_values(fit, trait = "fresh root yield"))
#'
#' @export
second_stage_values <- function(fits, trait = NULL) {
  if (bb_is_analysis(fits) || (is.list(fits) && !is.null(fits$adjusted_means))) {
    if (is.null(trait) || length(trait) != 1 || !nzchar(trait)) {
      stop("Give the Breedbase `trait` name of this analysis, e.g. trait = \"fresh root yield\".", call. = FALSE)
    }
    fits <- stats::setNames(list(fits), trait)
  }
  if (!is.list(fits) || is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("`fits` must be one analysis (with `trait`) or a list of analyses named by trait.", call. = FALSE)
  }
  out <- NULL
  random <- FALSE
  for (nm in names(fits)) {
    am <- fits[[nm]]$adjusted_means
    if (is.null(am) || !all(c("genotype", "adjusted_mean") %in% names(am))) {
      stop("`fits$", nm, "` has no adjusted means.", call. = FALSE)
    }
    random <- random || identical(fits[[nm]]$input$entries, "random")
    part <- tibble::tibble(accession = as.character(am$genotype), value = am$adjusted_mean)
    names(part)[2] <- nm
    out <- if (is.null(out)) part else dplyr::full_join(out, part, by = "accession")
  }
  attr(out, "estimate") <- if (random) "blup" else "blue"
  out
}

#' Save second-stage results to Breedbase as an analysis
#'
#' Saves genotype-level results (for example adjusted means across
#' locations from [second_stage_values()]) as a Breedbase analysis, the way
#' the Breedbase mixed-models tool does: a new project of type analysis with
#' one entry per accession, and each trait saved together with a statistic
#' term. The analysis then appears on the Analyses page of Breedbase.
#'
#' Every accession and trait is checked first, and all problems are listed
#' before anything is written.
#'
#' @param data A data frame with one row per accession and one numeric column
#'   per trait (as from [second_stage_values()]).
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns what would be sent.
#' @param analysis_name Name of the new analysis (must not exist yet).
#' @param program Breeding program of the analysis.
#' @param description Description of the analysis (required).
#' @param year Year of the analysis (default: this year).
#' @param traits Trait columns (default: all except `accession`). Their
#'   names are the Breedbase trait names.
#' @param accession Name of the accession column.
#' @param estimate `"blue"` (statistic "Adjusted means from BLUEs using LMER
#'   R") or `"blup"` (statistic "Adjusted means from BLUPs using LMER R"). By
#'   default taken from [second_stage_values()], else `"blue"`.
#' @param statistic,model_type To use another Breedbase statistic term
#'   (`"name|SGNSTAT:0000000"`) or model type (a `protocol_type` term; default
#'   `"mixed_model_lmer"`).
#' @param protocol Text describing the model, saved with the analysis
#'   (default: the description).
#' @param dry_run Logical. Check everything but write nothing.
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#'
#' @return Invisibly, a tibble with the analysis name and id, the numbers of
#'   accessions and values sent, the accessions read back, and the analysis
#'   page.
#'
#' @examples
#' prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
#' d <- prep$design
#' set.seed(1)
#' d$yield <- rnorm(nrow(d), 25, 3)
#' fit <- analyze_augmented_prep_two_stage(d, "yield")
#' means <- second_stage_values(fit, trait = "fresh root yield")
#' sent <- save_breedbase_analysis(means, analysis_name = "PREP_2026_stage2", program = "IITA",
#'                                 description = "Two-stage analysis, weighted", dry_run = TRUE)
#' str(sent$analysis_result_values[1])
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase_analysis(means, con, analysis_name = "PREP_2026_stage2", program = "IITA",
#'                         description = "Two-stage analysis, weighted")
#' }
#'
#' @export
save_breedbase_analysis <- function(data, con = NULL, analysis_name, program, description,
                                    year = format(Sys.Date(), "%Y"), traits = NULL, accession = "accession",
                                    estimate = NULL, statistic = NULL, model_type = "mixed_model_lmer",
                                    protocol = NULL, dry_run = FALSE, confirm = interactive()) {
  if (missing(analysis_name) || length(analysis_name) != 1 || is.na(analysis_name) || !nzchar(trimws(analysis_name))) {
    stop("Give the `analysis_name` of the new analysis, e.g. \"PREP_2026_stage2\".", call. = FALSE)
  }
  if (missing(program) || length(program) != 1 || !nzchar(program)) stop("Give the breeding `program`.", call. = FALSE)
  if (missing(description) || length(description) != 1 || is.na(description) || !nzchar(trimws(description))) {
    stop("Give a `description` of the analysis, e.g. the model.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  estimate <- match.arg(estimate %||% attr(data, "estimate") %||% "blue", c("blue", "blup"))
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (!accession %in% names(data)) {
    stop("Column '", accession, "' not found in `data`. Columns are: ", paste(names(data), collapse = ", "), ".",
         call. = FALSE)
  }
  if (is.null(traits)) traits <- setdiff(names(data), accession)
  if (!length(traits)) stop("No trait columns in `data`. Give them with `traits = c(...)`.", call. = FALSE)
  absent <- setdiff(traits, names(data))
  if (length(absent)) stop("Trait columns not found in `data`: ", paste(absent, collapse = ", "), call. = FALSE)
  for (t in traits) if (!is.numeric(data[[t]])) stop("Trait column '", t, "' must be numeric.", call. = FALSE)
  acc <- trimws(as.character(data[[accession]]))
  if (anyNA(acc) || any(!nzchar(acc))) stop("Some rows have no accession name.", call. = FALSE)
  if (anyDuplicated(acc)) {
    stop("Each accession must appear once; repeated: ", paste(utils::head(unique(acc[duplicated(acc)]), 5),
                                                               collapse = ", "), call. = FALSE)
  }
  statistic <- statistic %||% switch(estimate,
    blue = "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036",
    blup = "Adjusted means from BLUPs using LMER R|SGNSTAT:0000034")
  if (!grepl("|", statistic, fixed = TRUE)) {
    stop("`statistic` must be a Breedbase statistic term \"name|SGNSTAT:0000000\".", call. = FALSE)
  }
  protocol <- protocol %||% description

  build <- function(trait_full, user) {
    stamp <- format(Sys.time(), "%Y-%m-%d_%H:%M:%S")
    values <- list()
    for (i in seq_len(nrow(data))) {
      entry <- list()
      for (t in traits) {
        v <- data[[t]][i]
        if (!is.na(v)) entry[[trait_full[[t]]]] <- list(as.character(v), stamp, user, "", "")
      }
      if (length(entry)) values[[acc[i]]] <- entry
    }
    if (!length(values)) stop("There are no values to save (all are empty).", call. = FALSE)
    list(
      analysis_to_save_boolean = "yes",
      analysis_name = analysis_name,
      analysis_description = description,
      analysis_year = as.character(year),
      analysis_protocol = protocol,
      analysis_dataset_id = "",
      analysis_accession_names = acc,
      analysis_trait_names = unname(unlist(trait_full)),
      analysis_statistical_ontology_term = statistic,
      analysis_precomputed_design_optional = "",
      analysis_result_values = values,
      analysis_result_values_type = "analysis_result_values_match_accession_names",
      analysis_result_summary = list(method = paste("breedtools:", protocol)),
      analysis_model_name = paste(analysis_name, "model"),
      analysis_model_description = description,
      analysis_model_is_public = "yes",
      analysis_model_language = "R",
      analysis_model_type = model_type,
      analysis_model_properties = list(protocol = protocol, traits = unname(unlist(trait_full))),
      analysis_model_application_name = "breedtools",
      analysis_model_application_version = as.character(utils::packageVersion("breedtools"))
    )
  }

  if (is.null(con)) {
    if (!isTRUE(dry_run)) stop("`con` is needed to save; use breedbase_login() first, or `dry_run = TRUE`.", call. = FALSE)
    message("Dry run without a connection: returning the analysis that would be sent ",
            "(trait names are completed with their ontology ids when connected).")
    return(build(stats::setNames(as.list(traits), traits), "<your username>"))
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)

  cat("Checking the analysis before saving to ", con$url, "\n", sep = "")
  prog <- bb_find_program(con, program)
  cat("  [OK]   breeding program '", prog$program, "'\n", sep = "")
  problems <- character(0)
  if (!is.na(bb_find_analysis(con, analysis_name)) || !is.na(bb_find_study(con, analysis_name))) {
    problems <- c(problems, paste0("An analysis or trial named '", analysis_name,
                                   "' already exists. Choose another `analysis_name`."))
  } else {
    cat("  [OK]   analysis name '", analysis_name, "' is free\n", sep = "")
  }
  missing_acc <- setdiff(acc, bb_existing_germplasm(con, acc))
  if (length(missing_acc)) {
    problems <- c(problems, paste0(length(missing_acc), " of ", length(acc), " accessions are not in Breedbase: ",
                                   paste(utils::head(missing_acc, 10), collapse = ", ")))
  } else {
    cat("  [OK]   all ", length(acc), " accessions are in Breedbase\n", sep = "")
  }
  full <- attr(bb_find_variables(con, traits), "full_name")
  if (anyNA(full)) {
    problems <- c(problems, paste0("Trait(s) not found in Breedbase: ",
                                   paste0("'", traits[is.na(full)], "'", collapse = ", "),
                                   ". Use the trait name as shown in Breedbase."))
  } else {
    for (t in traits) cat("  [OK]   trait '", full[[t]], "'\n", sep = "")
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)
  cat("  [OK]   statistic '", statistic, "', model type '", model_type, "'\n", sep = "")

  payload <- build(as.list(full), con$username)
  payload$analysis_breeding_program_id <- prog$program_db_id
  n_values <- sum(vapply(payload$analysis_result_values, length, integer(1)))
  plan <- sprintf("analysis '%s' with %d accessions and %d values of %d trait%s in program '%s' at %s",
                  analysis_name, length(acc), n_values, length(traits), if (length(traits) > 1) "s" else "",
                  prog$program, con$url)
  if (isTRUE(dry_run)) {
    cat("Dry run: everything is ready to save ", plan, ". Nothing was written.\n", sep = "")
    return(invisible(payload))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }

  analysis_id <- bb_write_checked(
    what = paste0("save analysis '", analysis_name, "'"),
    write = function() {
      res <- bb_ajax_post(con, "ajax/analysis/store/json", payload, what = paste0("save analysis '", analysis_name, "'"),
                          json = c("analysis_accession_names", "analysis_trait_names", "analysis_result_values",
                                   "analysis_result_summary", "analysis_model_properties"))
      as.character(res$analysis_id %||% NA)
    },
    find = function() bb_find_analysis(con, analysis_name),
    attempts = 1
  )
  n_back <- bb_analysis_row_count(con, analysis_id)
  page <- paste0(con$url, "/analyses/", analysis_id)
  status <- if (!is.na(n_back) && n_back >= length(acc)) "[OK]  " else "[CHECK]"
  cat(sprintf("  %s analysis '%s' saved with %s of %d accessions\n", status, analysis_name,
              if (is.na(n_back)) "?" else as.character(n_back), length(acc)))
  cat("Open the analysis in Breedbase: ", page, "\n", sep = "")
  invisible(tibble::tibble(analysis_name = analysis_name, analysis_db_id = analysis_id,
                           accessions_sent = length(acc), values_sent = n_values,
                           accessions_found = n_back, analysis_page = page))
}

# POST to a Breedbase (non-BrAPI) AJAX endpoint as a form, like the web page
# does: lists and objects are sent as JSON text, and the login token as
# sgn_session_id. Errors come back as {"error": "..."} with HTTP 200.
bb_ajax_post <- function(con, path, params, what, json = character(0)) {
  bb_require_httr2()
  form <- stats::setNames(lapply(names(params), function(nm) {
    v <- params[[nm]]
    if (nm %in% json) {
      # Character vectors stay JSON arrays even with one element.
      if (is.list(v)) as.character(jsonlite::toJSON(v, auto_unbox = TRUE, null = "null", digits = NA))
      else as.character(jsonlite::toJSON(as.character(v)))
    } else {
      as.character(v)
    }
  }), names(params))
  form$sgn_session_id <- con$token
  req <- httr2::request(paste0(con$url, "/", path))
  req <- do.call(httr2::req_body_form, c(list(req), form))
  body <- bb_perform(req, what = what)
  if (!is.null(body$error)) {
    bb_stop(paste0("Could not ", what, ": ", paste(unlist(body$error), collapse = "; ")), code = 400L)
  }
  body
}

bb_ajax_get <- function(con, path, what, query = list(), allow_error = FALSE) {
  bb_require_httr2()
  req <- httr2::request(paste0(con$url, "/", path))
  req <- do.call(httr2::req_url_query, c(list(req), query, list(sgn_session_id = con$token, .multi = "explode")))
  body <- bb_perform(req, what = what)
  if (!isTRUE(allow_error) && is.list(body) && !is.null(names(body)) && !is.null(body$error)) {
    bb_stop(paste0("Could not ", what, ": ", paste(unlist(body$error), collapse = "; ")), code = 400L)
  }
  body
}

# Id of the user's analysis with this name (from the Analyses list), or NA.
bb_find_analysis <- function(con, name) {
  body <- tryCatch(bb_ajax_get(con, "ajax/analyses/by_user", what = "list your analyses"),
                   breedbase_error = function(e) NULL)
  if (is.null(body)) return(NA_character_)
  for (row in body$data) {
    link <- as.character(row[[1]])
    label <- sub("^.*>([^<]*)</a>.*$", "\\1", link)
    if (identical(label, name)) return(sub("^.*/analyses/([0-9]+).*$", "\\1", link))
  }
  NA_character_
}

# Entries of a saved analysis, read back from Breedbase.
bb_analysis_row_count <- function(con, analysis_id) {
  body <- tryCatch(bb_ajax_get(con, paste0("ajax/analysis/", analysis_id, "/retrieve"),
                               what = "read the saved analysis"),
                   breedbase_error = function(e) NULL)
  if (is.null(body) || is.null(body$data)) return(NA_integer_)
  max(length(body$data) - 1L, 0L)
}
