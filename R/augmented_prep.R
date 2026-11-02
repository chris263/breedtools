#' Generate an augmented p-rep design
#'
#' Creates an augmented partially replicated (p-rep) design inspired by the
#' Williams, Piepho, and Whitaker approach. Instead of replicated checks, a
#' subset of test entries is replicated at each location, while all entries
#' appear once at each location and the additional replicated entries are placed
#' into incomplete blocks. Candidate designs are scored to reduce pairwise
#' block concurrence, especially avoiding pairs that occur together more than
#' once.
#'
#' The incomplete block structure is independent from the field geometry:
#' `block_size` and `n_blocks_per_location` define block membership. By
#' default, `field_cols` lays the ordered plots into rows and columns. For a
#' rectangular block layout, use `rows_per_block` and `blocks_per_row` to place
#' multiple blocks side by side within each location.
#'
#' @param treatments Character vector of entry names.
#' @param locations Character vector of location names.
#' @param block_size Maximum final incomplete block size. Can be a single
#'   integer or a named vector indexed by location. Actual block sizes may
#'   differ by one plot when totals do not divide evenly.
#' @param prep_prop Optional proportion of entries to duplicate per location.
#'   If `NULL`, entries are split across locations so each entry is duplicated
#'   at approximately one location.
#' @param n_duplicate_per_location Optional number of entries duplicated at
#'   each location. Can be a single integer, a vector matching `locations`, or a
#'   named vector indexed by location.
#' @param replications_per_repeated_entry Total number of plots for entries
#'   selected for p-rep replication at a location. The default `2` gives one
#'   base plot plus one extra plot. Can be a single integer or a named vector
#'   indexed by location.
#' @param n_blocks_per_location Optional number of incomplete blocks at each
#'   location. If `NULL`, the number of blocks is the minimum needed so no block
#'   is larger than `block_size`. Can be a single integer or a named vector
#'   indexed by location.
#' @param n_candidates Maximum number of candidate designs to generate and
#'   score. The search stops early when a candidate reaches the best possible
#'   score (no entry repeated within a block and no pair of entries sharing more
#'   than one block); `summary$n_candidates_evaluated` reports how many were
#'   generated.
#' @param seed Optional random seed.
#' @param field_cols Optional number of field columns for row/column layout.
#'   Can be a single integer or a named vector indexed by location.
#' @param rows_per_block Optional number of physical field rows each incomplete
#'   block should occupy. Can be a single integer or a named vector for selected
#'   locations. Must be used with `blocks_per_row` for those locations.
#' @param blocks_per_row Optional number of incomplete blocks to place side by
#'   side in each field row band. Can be a single integer or a named vector for
#'   selected locations. Must be used with `rows_per_block` for those locations.
#' @param serpentine Logical. If `TRUE`, row/column coordinates use serpentine
#'   ordering within each location.
#' @param output_file Optional tab-delimited output file for the complete
#'   design.
#'
#' @return A list with design inputs, summary diagnostics, complete plot-level
#' design, Breedbase-compatible design columns, block summary, replication
#' summary, repeated-entry groups, and pairwise concurrence matrix. In the
#' complete design, `repeated_at` is the location where the entry receives its
#' extra p-rep plot(s), or `"not_repeated"` when the entry is not selected for
#' extra replication.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2", "Loc3"),
#'   block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 6),
#'   n_blocks_per_location = c(Loc1 = 4, Loc2 = 4, Loc3 = 4),
#'   n_duplicate_per_location = c(Loc1 = 4, Loc2 = 2, Loc3 = 4),
#'   replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3, Loc3 = 2),
#'   n_candidates = 5,
#'   seed = 123,
#'   field_cols = 8
#' )
#'
#' design$summary
#' head(design$design)
#'
#' @export
augmented_prep_design <- function(
    treatments,
    locations,
    block_size = 8,
    prep_prop = NULL,
    n_duplicate_per_location = NULL,
    replications_per_repeated_entry = 2,
    n_blocks_per_location = NULL,
    n_candidates = 1000,
    seed = NULL,
    field_cols = NULL,
    rows_per_block = NULL,
    blocks_per_row = NULL,
    serpentine = TRUE,
    output_file = NULL
) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  inputs <- prep_validate_inputs(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    replications_per_repeated_entry = replications_per_repeated_entry,
    n_blocks_per_location = n_blocks_per_location,
    n_candidates = n_candidates,
    field_cols = field_cols,
    rows_per_block = rows_per_block,
    blocks_per_row = blocks_per_row
  )

  setup <- prep_candidate_setup(inputs)

  # Candidates are generated on integer entry ids only; the plot-level design
  # is built once, for the best candidate. The score has a known lower bound
  # (no pair of entries shares a block more than once and no entry repeats in
  # a block), so the search stops as soon as a candidate reaches it.
  best_candidate <- NULL
  n_evaluated <- 0L
  for (i in seq_len(inputs$n_candidates)) {
    candidate <- prep_make_candidate(setup)
    n_evaluated <- i

    if (is.null(best_candidate) || candidate$score < best_candidate$score) {
      best_candidate <- candidate
    }

    if (best_candidate$score <= setup$score_lower_bound) {
      break
    }
  }

  best <- prep_build_candidate_design(
    candidate = best_candidate,
    inputs = inputs,
    setup = setup,
    serpentine = serpentine
  )
  best$summary$n_candidates_evaluated <- n_evaluated

  block_summary <- best$design |>
    dplyr::group_by(.data$location, .data$block, .data$block_uid) |>
    dplyr::summarise(n_plots = dplyr::n(), .groups = "drop")

  replication_summary <- best$design |>
    dplyr::group_by(.data$all_entries) |>
    dplyr::summarise(
      total_plots = dplyr::n(),
      locations = paste(unique(.data$location), collapse = ", "),
      repeated_at = paste(unique(.data$repeated_at[.data$repeated_at != "not_repeated"]), collapse = ", "),
      .groups = "drop"
    )
  replication_summary$repeated_at[replication_summary$repeated_at == ""] <- "not_repeated"

  breedbase_design <- best$design |>
    dplyr::select(dplyr::all_of(c("plots", "block", "all_entries", "rep", "is_control")))

  result <- list(
    input = inputs,
    summary = best$summary,
    design = best$design,
    breedbase_design = breedbase_design,
    block_summary = block_summary,
    replication_summary = replication_summary,
    repeated_groups = best$repeated_groups,
    concurrence_matrix = best$concurrence_matrix
  )
  class(result) <- c("augmented_prep_design", class(result))

  if (!is.null(output_file)) {
    utils::write.table(result$design, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  }

  result
}

#' Calculate pairwise concurrence for a p-rep design
#'
#' Counts how many times each pair of entries appears together in the same
#' incomplete block across all locations.
#'
#' @param design A p-rep design data frame from `augmented_prep_design()`.
#' @param treatments Optional treatment names defining matrix order.
#'
#' @return A symmetric integer matrix of pairwise block concurrences.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2"),
#'   block_size = 4,
#'   n_candidates = 5,
#'   seed = 1
#' )
#' conc <- augmented_prep_concurrence(design$design)
#' table(conc[upper.tri(conc)])
#'
#' @export
augmented_prep_concurrence <- function(design, treatments = NULL) {
  if (is.null(treatments)) {
    treatments <- sort(unique(as.character(design$all_entries)))
  }

  treatments <- as.character(treatments)
  entry_idx <- match(as.character(design$all_entries), treatments)
  if (anyNA(entry_idx)) {
    stop("`design` contains entries not listed in `treatments`.", call. = FALSE)
  }
  blocks <- interaction(design$location, design$block, drop = TRUE)

  # Block x entry presence matrix; crossprod counts shared blocks per pair.
  presence <- matrix(0, nrow = nlevels(blocks), ncol = length(treatments))
  presence[cbind(as.integer(blocks), entry_idx)] <- 1
  mat <- crossprod(presence)
  diag(mat) <- 0
  storage.mode(mat) <- "integer"
  dimnames(mat) <- list(treatments, treatments)
  mat
}


#' Analyze an augmented p-rep design
#'
#' Fits the fixed-effects incomplete-block model:
#' `y = mean + location + block(location) + genotype + error`.
#' Sums of squares are sequential, in the order above. Degrees of freedom in
#' the formatted ANOVA are those estimable from the observed (non-missing)
#' plots, so F-tests stay valid when plots are missing. For a complete,
#' connected design they equal the design formulas reported in `anova_df`,
#' which use total plots (`n`), locations (`l`), incomplete blocks (`b`), and
#' genotypes (`v`). Terms with a single level in the data (for example the
#' location term in a single-location trial) are left out of the model and
#' shown with 0 df.
#'
#' @param data A complete augmented p-rep design data frame, usually
#'   `design$design` from `augmented_prep_design()`, or an object returned by
#'   `augmented_prep_design()`, with a response column added.
#' @param response Name of the numeric response column.
#' @param genotype_col Name of the genotype/treatment column.
#' @param location_col Name of the location column.
#' @param block_col Name of the incomplete block column. Block numbers may
#'   repeat across locations; blocks are treated as nested within locations.
#' @param entry_type_col Optional entry type column. Used only for genotype
#'   summaries.
#' @param extra_repeat_col Optional extra-repeat indicator column. Used only
#'   for genotype summaries.
#' @param repeated_at_col Optional column identifying where each entry is
#'   repeated. Used only for genotype summaries.
#' @param control_col Optional control flag column. Used only for genotype
#'   summaries.
#'
#' @return An object of class `augmented_prep_analysis`: a list with
#' \describe{
#'   \item{input, formula, model}{The column mapping, model formula, and the
#'     `lm` fit.}
#'   \item{anova, anova_df, anova_formatted}{Sequential ANOVA (estimable df),
#'     design df formulas, and a formatted ANOVA table.}
#'   \item{plot_results}{The input data with fitted values and residuals.}
#'   \item{genotype_summary}{One row per genotype with raw, fitted, and
#'     adjusted means and the standard error of the adjusted mean.}
#'   \item{adjusted_means}{Genotypes ranked by adjusted mean (best first).}
#'   \item{diagnostics}{Error df, residual standard deviation, R-squared,
#'     `mean_sed` (mean standard error of a difference between two adjusted
#'     means), `lsd_5pct` (the matching least significant difference), and
#'     flags for rank deficiency or a saturated fit.}
#'   \item{message}{A note when F-tests or adjusted means are not estimable.}
#' }
#'
#' @section Adjusted means:
#' Genotype effects are fixed, so the adjusted means are best linear unbiased
#' estimates (BLUEs): the genotype effect plus the average location and block
#' effect over all observed plots, i.e. the expected value of the genotype in
#' an average plot of the trial. They remove the advantage or penalty of the
#' blocks and locations where each entry happened to be planted, so use
#' `adjusted_mean` (not `raw_mean`) to rank and select entries. Differences
#' between adjusted means are the same as those from `emmeans` on the fitted
#' model.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(8)),
#'   locations = c("Loc1", "Loc2"),
#'   block_size = 4,
#'   n_duplicate_per_location = c(Loc1 = 2, Loc2 = 2),
#'   n_candidates = 3,
#'   seed = 123
#' )
#'
#' design_data <- design$design
#' design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
#' analysis <- analyze_augmented_prep_design(design_data, response = "yield")
#' analysis
#' head(analysis$adjusted_means)
#'
#' @export
analyze_augmented_prep_design <- function(
    data,
    response,
    genotype_col = "all_entries",
    location_col = "location",
    block_col = "block",
    entry_type_col = "entry_type",
    extra_repeat_col = "is_extra_repeat",
    repeated_at_col = "repeated_at",
    control_col = "is_control"
) {
  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) {
    data <- data$design
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or an augmented p-rep design object.", call. = FALSE)
  }

  required_cols <- c(response, genotype_col, location_col, block_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  design_structure_data <- tibble::as_tibble(data) |>
    dplyr::transmute(
      genotype = as.factor(.data[[genotype_col]]),
      location = as.factor(.data[[location_col]]),
      block = as.factor(.data[[block_col]]),
      block_nested = interaction(.data[[location_col]], .data[[block_col]], drop = TRUE),
      entry_type = if (entry_type_col %in% names(data)) {
        as.character(.data[[entry_type_col]])
      } else {
        NA_character_
      },
      is_extra_repeat = if (extra_repeat_col %in% names(data)) {
        as.integer(.data[[extra_repeat_col]])
      } else {
        NA_integer_
      },
      repeated_at = if (repeated_at_col %in% names(data)) {
        as.character(.data[[repeated_at_col]])
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(data)) {
        as.integer(.data[[control_col]])
      } else {
        NA_integer_
      }
    ) |>
    dplyr::filter(
      !is.na(.data$genotype),
      !is.na(.data$location),
      !is.na(.data$block),
      !is.na(.data$block_nested)
    )

  analysis_data <- tibble::as_tibble(data) |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::transmute(
      .row_id = .data$.row_id,
      response = suppressWarnings(as.numeric(.data[[response]])),
      genotype = as.factor(.data[[genotype_col]]),
      location = as.factor(.data[[location_col]]),
      block = as.factor(.data[[block_col]]),
      block_nested = interaction(.data[[location_col]], .data[[block_col]], drop = TRUE),
      entry_type = if (entry_type_col %in% names(data)) {
        as.character(.data[[entry_type_col]])
      } else {
        NA_character_
      },
      is_extra_repeat = if (extra_repeat_col %in% names(data)) {
        as.integer(.data[[extra_repeat_col]])
      } else {
        NA_integer_
      },
      repeated_at = if (repeated_at_col %in% names(data)) {
        as.character(.data[[repeated_at_col]])
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(data)) {
        as.integer(.data[[control_col]])
      } else {
        NA_integer_
      }
    ) |>
    dplyr::filter(
      !is.na(.data$response),
      !is.na(.data$genotype),
      !is.na(.data$location),
      !is.na(.data$block),
      !is.na(.data$block_nested)
    )

  if (nrow(analysis_data) == 0) {
    stop("No complete observations are available for analysis.", call. = FALSE)
  }

  if (dplyr::n_distinct(analysis_data$genotype) < 2) {
    stop("At least two genotypes are required for analysis.", call. = FALSE)
  }

  analysis_data <- droplevels(analysis_data)

  # Only fit terms that vary in the observed data: a single-location trial has
  # no location term, and a location with one block has no block term.
  model_terms <- c(
    if (nlevels(analysis_data$location) > 1) "location",
    if (nlevels(analysis_data$block_nested) > nlevels(analysis_data$location)) "block_nested",
    "genotype"
  )
  model_formula <- stats::as.formula(paste("response ~", paste(model_terms, collapse = " + ")))
  design_df <- prep_augmented_prep_design_df(design_structure_data)
  fit <- arc_lm_drop_aliased_terms(model_formula, analysis_data)
  anova_result <- arc_lm_anova_table(fit)
  anova_table <- anova_result$table
  anova_formatted <- prep_format_augmented_prep_anova(
    anova_table = anova_table,
    response = analysis_data$response
  )

  plot_results <- tibble::as_tibble(data) |>
    dplyr::mutate(
      .analysis_response = suppressWarnings(as.numeric(.data[[response]]))
    )

  plot_results$fitted <- NA_real_
  plot_results$residual <- NA_real_
  plot_results$fitted[analysis_data$.row_id] <- stats::fitted(fit)
  plot_results$residual[analysis_data$.row_id] <- stats::residuals(fit)

  genotype_summary <- plot_results |>
    dplyr::group_by(.data[[genotype_col]]) |>
    dplyr::summarise(
      genotype = dplyr::first(as.character(.data[[genotype_col]])),
      entry_type = if (entry_type_col %in% names(plot_results)) {
        paste(sort(unique(as.character(.data[[entry_type_col]]))), collapse = ", ")
      } else {
        NA_character_
      },
      repeated_at = if (repeated_at_col %in% names(plot_results)) {
        paste(sort(unique(as.character(.data[[repeated_at_col]]))), collapse = ", ")
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(plot_results)) {
        dplyr::first(as.integer(.data[[control_col]]))
      } else {
        NA_integer_
      },
      n_observations = sum(!is.na(.data$.analysis_response)),
      n_locations = if (location_col %in% names(plot_results)) {
        dplyr::n_distinct(.data[[location_col]][!is.na(.data$.analysis_response)])
      } else {
        NA_integer_
      },
      n_extra_repeats = if (extra_repeat_col %in% names(plot_results)) {
        sum(as.integer(.data[[extra_repeat_col]]) == 1L & !is.na(.data$.analysis_response), na.rm = TRUE)
      } else {
        NA_integer_
      },
      raw_mean = mean(.data$.analysis_response, na.rm = TRUE),
      fitted_mean = mean(.data$fitted, na.rm = TRUE),
      residual_mean = mean(.data$residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::select(-dplyr::all_of(genotype_col))

  genotype_summary$raw_mean[is.nan(genotype_summary$raw_mean)] <- NA_real_
  genotype_summary$fitted_mean[is.nan(genotype_summary$fitted_mean)] <- NA_real_
  genotype_summary$residual_mean[is.nan(genotype_summary$residual_mean)] <- NA_real_

  adjusted_fit <- prep_adjusted_means(fit, analysis_data)
  adjusted_parts <- arc_attach_adjusted_means(genotype_summary, adjusted_fit$adjusted_means)
  genotype_summary <- adjusted_parts$genotype_summary
  precision <- adjusted_fit$precision

  anova_warning <- anova_result$warning
  if (is.null(anova_warning)) {
    anova_warning <- NA_character_
  }

  aliased_coefficients_dropped <- attr(fit, "aliased_coefficients_dropped")
  if (is.null(aliased_coefficients_dropped)) {
    aliased_coefficients_dropped <- 0L
  }

  diagnostics <- tibble::tibble(
    n_observations = stats::nobs(fit),
    n_locations = dplyr::n_distinct(analysis_data$location),
    n_blocks = dplyr::n_distinct(analysis_data$block_nested),
    n_genotypes = dplyr::n_distinct(analysis_data$genotype),
    model_rank = fit$rank,
    residual_df = stats::df.residual(fit),
    anova_error_df = stats::df.residual(fit),
    design_error_df = design_df$display_df[design_df$term == "Residuals"],
    aliased_coefficients_dropped = aliased_coefficients_dropped,
    sigma = summary(fit)$sigma,
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    rank_deficient = fit$rank < length(stats::coef(fit)),
    saturated_or_near_perfect_fit = isTRUE(stats::df.residual(fit) == 0) ||
      isTRUE(!is.na(summary(fit)$sigma) && summary(fit)$sigma < sqrt(.Machine$double.eps)),
    mean_sed = precision$mean_sed,
    lsd_5pct = precision$lsd_5pct,
    anova_warning = anova_warning
  )

  analysis_message <- arc_augmented_analysis_message(diagnostics)
  analysis_message <- arc_append_estimability_message(analysis_message, adjusted_fit$adjusted_means)

  result <- list(
    input = list(
      response = response,
      genotype_col = genotype_col,
      location_col = location_col,
      block_col = block_col,
      entry_type_col = entry_type_col,
      extra_repeat_col = extra_repeat_col,
      repeated_at_col = repeated_at_col,
      control_col = control_col
    ),
    formula = model_formula,
    model = fit,
    anova = anova_table,
    anova_df = design_df,
    anova_formatted = anova_formatted,
    plot_results = plot_results |> dplyr::select(-dplyr::all_of(".analysis_response")),
    genotype_summary = genotype_summary,
    adjusted_means = adjusted_parts$adjusted_means,
    diagnostics = diagnostics,
    message = analysis_message
  )

  class(result) <- c("augmented_prep_analysis", class(result))
  result
}

#' Print augmented p-rep analysis
#'
#' @param x An object returned by `analyze_augmented_prep_design()`.
#' @param ... Additional arguments passed to `print()`.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.augmented_prep_analysis <- function(x, ...) {
  cat("Augmented p-rep design analysis\n\n")
  print(x$anova_formatted, ...)
  arc_print_adjusted_means(x)
  if (!is.null(x$message) && !is.na(x$message) && nzchar(x$message)) {
    cat("\n", x$message, "\n", sep = "")
  }
  invisible(x)
}


#' Plot an augmented p-rep design grid
#'
#' Draws one field grid per location for an augmented p-rep design, with
#' plot cells, row and column labels, and thick outlines around incomplete
#' blocks.
#'
#' @param design An object returned by `augmented_prep_design()` or its
#'   `design` data frame.
#' @param locations Optional character vector of locations to plot. Defaults
#'   to all locations in design order.
#' @param label Cell label to show. Use `"all_entries"`, `"entry_type"`,
#'   `"repeated_at"`, `"block"`, or `"none"`.
#' @param fill Cell coloring variable. Use `"entry_type"`,
#'   `"is_extra_repeat"`, `"repeated_at"`, or `"block"`.
#' @param show_plot_labels Logical. If `TRUE`, labels are printed inside
#'   cells unless `label = "none"`.
#' @param show_block_labels Logical. If `TRUE`, block numbers are printed.
#'   Defaults to `FALSE` so block labels do not overlap accession labels.
#' @param block_label_position Position for block labels when
#'   `show_block_labels = TRUE`. Use `"background"` for a faint watermark
#'   behind accession labels or `"center"` for a regular center label.
#' @param label_cex Optional text size for plot-cell labels.
#' @param block_label_cex Optional text size for block labels.
#' @param block_label_col Color for block labels.
#' @param block_lwd Line width for block boundaries.
#' @param main Optional plot title. If multiple locations are plotted, a
#'   single value is combined with the location name; a vector matching the
#'   number of locations is used panel by panel.
#' @param legend Logical. If `TRUE`, draw a legend for cell colors.
#' @param ... Additional arguments passed to `graphics::plot.window()`.
#'
#' @return Invisibly returns the plotted design data frame with internal
#'   plotting coordinates.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2", "Loc3"),
#'   block_size = 5,
#'   n_candidates = 5,
#'   seed = 123,
#'   field_cols = 5
#' )
#'
#' plot_augmented_prep_design(design)
#'
#' @export
plot_augmented_prep_design <- function(
    design,
    locations = NULL,
    label = c("all_entries", "entry_type", "repeated_at", "block", "none"),
    fill = c("entry_type", "is_extra_repeat", "repeated_at", "block"),
    show_plot_labels = TRUE,
    show_block_labels = FALSE,
    block_label_position = c("background", "center"),
    label_cex = NULL,
    block_label_cex = NULL,
    block_label_col = NULL,
    block_lwd = 2,
    main = NULL,
    legend = TRUE,
    ...
) {
  label <- match.arg(label)
  fill <- match.arg(fill)
  block_label_position <- match.arg(block_label_position)

  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) {
    design <- design$design
  }

  if (!is.data.frame(design)) {
    stop("`design` must be a design data frame or an augmented p-rep design object.", call. = FALSE)
  }

  required_cols <- c("location", "row", "col", "block")
  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (label != "none" && !label %in% names(design)) {
    stop("Column `", label, "` was not found in `design`.", call. = FALSE)
  }
  if (!fill %in% names(design)) {
    stop("Column `", fill, "` was not found in `design`.", call. = FALSE)
  }

  design <- as.data.frame(design, stringsAsFactors = FALSE)
  design$location <- as.character(design$location)
  design$row <- suppressWarnings(as.numeric(design$row))
  design$col <- suppressWarnings(as.numeric(design$col))

  if (any(is.na(design$row)) || any(is.na(design$col))) {
    stop(
      "Plotting augmented p-rep designs requires row and column coordinates. ",
      "Call `augmented_prep_design(..., field_cols = ...)` before plotting.",
      call. = FALSE
    )
  }

  if (is.null(locations)) {
    locations <- unique(design$location)
  } else {
    locations <- as.character(locations)
  }

  missing_locations <- setdiff(locations, unique(design$location))
  if (length(missing_locations) > 0) {
    stop("Location(s) not found in `design`: ", paste(missing_locations, collapse = ", "), call. = FALSE)
  }

  design <- design[design$location %in% locations, , drop = FALSE]
  design$.panel_index <- match(design$location, locations)
  design <- design[order(design$.panel_index, design$row, design$col), , drop = FALSE]
  design$.fill_value <- prep_plot_fill_values(design, fill)

  fill_levels <- sort(unique(design$.fill_value))
  fill_colors <- grDevices::hcl.colors(length(fill_levels), palette = "Set 3")
  names(fill_colors) <- fill_levels
  plot_args <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_locations <- length(locations)
  panel_cols <- ceiling(sqrt(n_locations))
  panel_rows <- ceiling(n_locations / panel_cols)
  graphics::par(
    mfrow = c(panel_rows, panel_cols),
    mar = c(4.5, 4.5, 4, if (isTRUE(legend)) 7 else 2),
    xpd = NA
  )

  plotted <- vector("list", n_locations)
  for (i in seq_along(locations)) {
    location <- locations[i]
    location_design <- design[design$location == location, , drop = FALSE]
    panel_main <- if (is.null(main)) {
      location
    } else if (length(main) == n_locations) {
      main[i]
    } else {
      paste(main[1], location, sep = " - ")
    }

    plotted[[i]] <- prep_plot_prep_location(
      design = location_design,
      fill_colors = fill_colors,
      label = label,
      show_plot_labels = show_plot_labels,
      show_block_labels = show_block_labels,
      block_label_position = block_label_position,
      label_cex = label_cex,
      block_label_cex = block_label_cex,
      block_label_col = block_label_col,
      block_lwd = block_lwd,
      main = panel_main,
      legend = legend,
      fill_title = fill,
      plot_args = plot_args
    )
  }

  rownames(design) <- NULL
  invisible(tibble::as_tibble(do.call(rbind, plotted)))
}

prep_plot_fill_values <- function(design, fill) {
  value <- as.character(design[[fill]])
  if (fill == "entry_type") {
    value <- ifelse(
      value == "p_rep_entry_at_this_location",
      "p-rep entry here",
      "single entry here"
    )
  }
  if (fill == "is_extra_repeat") {
    value <- ifelse(as.integer(design[[fill]]) == 1L, "extra repeat", "base plot")
  }
  value
}

prep_plot_prep_location <- function(design, fill_colors, label, show_plot_labels,
                                    show_block_labels, block_label_position,
                                    label_cex, block_label_cex, block_label_col,
                                    block_lwd, main, legend, fill_title,
                                    plot_args) {
  rows <- sort(unique(design$row))
  cols <- sort(unique(design$col))
  design$.plot_row <- match(design$row, rows)
  design$.plot_col <- match(design$col, cols)

  if (is.null(label_cex)) {
    label_cex <- max(0.25, min(0.8, 7 / max(length(rows), length(cols))))
  }

  graphics::plot.new()
  do.call(
    graphics::plot.window,
    c(
      list(
        xlim = c(0.5, length(cols) + 0.5),
        ylim = c(length(rows) + 0.5, 0.5),
        xaxs = "i",
        yaxs = "i",
        asp = 1
      ),
      plot_args
    )
  )

  graphics::rect(
    xleft = design$.plot_col - 0.5,
    ybottom = design$.plot_row - 0.5,
    xright = design$.plot_col + 0.5,
    ytop = design$.plot_row + 0.5,
    col = fill_colors[design$.fill_value],
    border = "grey70",
    lwd = 0.7
  )

  graphics::abline(v = seq(0.5, length(cols) + 0.5, by = 1), col = "grey80", lwd = 0.5)
  graphics::abline(h = seq(0.5, length(rows) + 0.5, by = 1), col = "grey80", lwd = 0.5)
  prep_plot_prep_block_boundaries(design, block_lwd = block_lwd)

  if (isTRUE(show_block_labels)) {
    block_split <- split(design, design$block)
    for (block_design in block_split) {
      block_cex <- block_label_cex
      block_col <- block_label_col
      if (is.null(block_cex)) {
        block_cex <- if (block_label_position == "background") {
          max(1.1, min(3.5, 0.6 * sqrt(nrow(block_design))))
        } else {
          0.85
        }
      }
      if (is.null(block_col)) {
        block_col <- if (block_label_position == "background") {
          grDevices::adjustcolor("black", alpha.f = 0.18)
        } else {
          "black"
        }
      }
      graphics::text(
        x = mean(block_design$.plot_col),
        y = mean(block_design$.plot_row),
        labels = paste0("B", unique(block_design$block)[1]),
        cex = block_cex,
        font = 2,
        col = block_col
      )
    }
  }

  if (isTRUE(show_plot_labels) && label != "none") {
    graphics::text(
      x = design$.plot_col,
      y = design$.plot_row,
      labels = as.character(design[[label]]),
      cex = label_cex,
      col = "black"
    )
  }

  graphics::axis(3, at = seq_along(cols), labels = cols, las = 2, tick = FALSE)
  graphics::axis(2, at = seq_along(rows), labels = rows, las = 1, tick = FALSE)
  graphics::mtext("Column", side = 3, line = 2.7)
  graphics::mtext("Row", side = 2, line = 3)
  graphics::title(main = main)
  graphics::box()

  if (isTRUE(legend)) {
    graphics::legend(
      x = length(cols) + 1,
      y = 0.5,
      legend = names(fill_colors),
      fill = fill_colors,
      border = "grey60",
      title = fill_title,
      bty = "n",
      xjust = 0,
      yjust = 0,
      cex = 0.75
    )
  }

  design
}

prep_plot_prep_block_boundaries <- function(design, block_lwd) {
  block_lookup <- stats::setNames(
    as.character(design$block),
    paste(design$.plot_row, design$.plot_col, sep = ":")
  )
  same_block <- function(row, col, block) {
    key <- paste(row, col, sep = ":")
    !is.na(block_lookup[key]) && identical(unname(block_lookup[key]), as.character(block))
  }

  for (i in seq_len(nrow(design))) {
    row <- design$.plot_row[i]
    col <- design$.plot_col[i]
    block <- design$block[i]
    left <- col - 0.5
    right <- col + 0.5
    top <- row - 0.5
    bottom <- row + 0.5

    if (!same_block(row, col - 1L, block)) {
      graphics::segments(left, top, left, bottom, lwd = block_lwd, col = "black")
    }
    if (!same_block(row, col + 1L, block)) {
      graphics::segments(right, top, right, bottom, lwd = block_lwd, col = "black")
    }
    if (!same_block(row - 1L, col, block)) {
      graphics::segments(left, top, right, top, lwd = block_lwd, col = "black")
    }
    if (!same_block(row + 1L, col, block)) {
      graphics::segments(left, bottom, right, bottom, lwd = block_lwd, col = "black")
    }
  }
}


# Adjusted genotype means for the p-rep model: the genotype effect plus the
# average location/block effect over all observed plots (the prediction for
# an "average plot" of the trial). Estimability is checked against the full
# model matrix, so genotypes whose effect is aliased get NA.
prep_adjusted_means <- function(fit, analysis_data) {
  x_full <- stats::model.matrix(fit$terms, fit$model)
  y <- stats::model.response(fit$model)
  assign <- attr(x_full, "assign")
  genotype_term <- match("genotype", attr(fit$terms, "term.labels"))
  genotype_cols <- which(assign == genotype_term)
  genotype_levels <- levels(analysis_data$genotype)

  L <- matrix(colMeans(x_full), nrow = length(genotype_levels), ncol = ncol(x_full), byrow = TRUE)
  L[, genotype_cols] <- 0
  col_level <- match(colnames(x_full)[genotype_cols], paste0("genotype", genotype_levels))
  L[cbind(col_level, genotype_cols)] <- 1

  qr_full <- qr(x_full)
  est <- arc_linear_estimates(qr_full, qr.coef(qr_full, y), L)
  error_df <- stats::df.residual(fit)
  sigma2 <- if (error_df > 0) sum(stats::residuals(fit)^2) / error_df else NA_real_
  std_error <- sqrt(sigma2 * colSums(est$W^2))
  estimate <- est$estimate
  estimate[!est$estimable] <- NA_real_
  std_error[!est$estimable] <- NA_real_

  list(
    adjusted_means = tibble::tibble(
      genotype = genotype_levels,
      adjusted_mean = estimate,
      std_error = std_error,
      estimable = est$estimable
    ),
    precision = arc_mean_precision(
      est$W[, est$estimable, drop = FALSE],
      rep(0, sum(est$estimable)),
      sigma2,
      error_df
    )
  )
}

prep_format_augmented_prep_anova <- function(anova_table, response) {
  source_map <- data.frame(
    term = c(
      "location",
      "block_nested",
      "genotype",
      "Residuals"
    ),
    source = c(
      "Locations",
      "Blocks, nested within locations",
      "Genotypes",
      "Error"
    ),
    order = seq_len(4),
    stringsAsFactors = FALSE
  )

  # Degrees of freedom are those estimable from the observed plots. A term that
  # was not fitted (e.g. locations in a single-location trial) has 0 df.
  out <- source_map |>
    dplyr::left_join(anova_table, by = "term") |>
    dplyr::arrange(.data$order) |>
    dplyr::mutate(
      df = dplyr::coalesce(.data$df, 0),
      sum_sq = dplyr::coalesce(.data$sum_sq, 0),
      mean_sq = dplyr::if_else(.data$df > 0, .data$sum_sq / .data$df, NA_real_)
    )

  error_df <- out$df[out$term == "Residuals"]
  error_mean_sq <- out$mean_sq[out$term == "Residuals"]

  out <- out |>
    dplyr::mutate(
      f_value = dplyr::if_else(
        .data$term != "Residuals" &
          !is.na(error_mean_sq) &
          error_mean_sq > 0 &
          .data$df > 0,
        .data$mean_sq / error_mean_sq,
        NA_real_
      ),
      p_value = dplyr::if_else(
        !is.na(.data$f_value) &
          !is.na(error_df) &
          error_df > 0,
        stats::pf(.data$f_value, df1 = .data$df, df2 = error_df, lower.tail = FALSE),
        NA_real_
      )
    ) |>
    dplyr::transmute(
      Source = .data$source,
      Df = .data$df,
      `Sum Sq` = .data$sum_sq,
      `Mean Sq` = .data$mean_sq,
      `F value` = .data$f_value,
      `Pr(>F)` = .data$p_value,
      Significance = arc_significance_stars(.data$p_value)
    )

  out <- arc_clean_anova_numbers(out)

  corrected_total <- tibble::tibble(
    Source = "Corrected total",
    Df = length(response) - 1,
    `Sum Sq` = sum((response - mean(response))^2),
    `Mean Sq` = NA_real_,
    `F value` = NA_real_,
    `Pr(>F)` = NA_real_,
    Significance = ""
  )

  dplyr::bind_rows(out, corrected_total)
}

prep_augmented_prep_design_df <- function(data) {
  n <- nrow(data)
  l <- dplyr::n_distinct(data$location)
  b <- dplyr::n_distinct(data$block_nested)
  v <- dplyr::n_distinct(data$genotype)

  source_df <- c(
    l - 1,
    b - l,
    v - 1
  )
  corrected_total_df <- n - 1
  error_df <- corrected_total_df - sum(source_df)

  tibble::tibble(
    term = c(
      "location",
      "block_nested",
      "genotype",
      "Residuals"
    ),
    display_df = c(
      source_df,
      error_df
    ),
    n = n,
    l = l,
    b = b,
    v = v,
    corrected_total_df = corrected_total_df,
    error_df_formula = paste0(
      "(",
      n,
      " - 1) - [(",
      l,
      " - 1) + (",
      b,
      " - ",
      l,
      ") + (",
      v,
      " - 1)]"
    )
  )
}


prep_validate_inputs <- function(treatments, locations, block_size, prep_prop,
                                 n_duplicate_per_location,
                                 replications_per_repeated_entry,
                                 n_blocks_per_location, n_candidates,
                                 field_cols, rows_per_block,
                                 blocks_per_row) {
  treatments <- as.character(treatments)
  locations <- as.character(locations)
  treatments <- treatments[!is.na(treatments) & treatments != ""]
  locations <- locations[!is.na(locations) & locations != ""]
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

  if (length(treatments) < 2) {
    stop("At least two treatments are required.", call. = FALSE)
  }
  if (length(locations) < 1) {
    stop("At least one location is required.", call. = FALSE)
  }
  if (anyDuplicated(treatments)) {
    stop("`treatments` must be unique.", call. = FALSE)
  }
  if (anyDuplicated(locations)) {
    stop("`locations` must be unique.", call. = FALSE)
  }
  if (!is.null(prep_prop) && (length(prep_prop) != 1 || is.na(prep_prop) ||
    prep_prop <= 0 || prep_prop > 1)) {
    stop("`prep_prop` must be a single number greater than 0 and no larger than 1.", call. = FALSE)
  }

  n_duplicate_per_location <- prep_normalize_duplicate_counts(
    treatments, locations, prep_prop, n_duplicate_per_location
  )
  block_size <- prep_normalize_location_counts(block_size, locations, "block_size")
  replications_per_repeated_entry <- prep_normalize_location_counts(
    replications_per_repeated_entry,
    locations,
    "replications_per_repeated_entry"
  )
  if (any(replications_per_repeated_entry < 2L)) {
    stop("`replications_per_repeated_entry` must be at least 2.", call. = FALSE)
  }
  n_blocks_per_location <- prep_normalize_optional_location_counts(
    n_blocks_per_location,
    locations,
    "n_blocks_per_location"
  )
  field_cols <- prep_normalize_field_cols(field_cols, locations)
  rows_per_block <- prep_normalize_optional_partial_location_counts(
    rows_per_block,
    locations,
    "rows_per_block"
  )
  blocks_per_row <- prep_normalize_optional_partial_location_counts(
    blocks_per_row,
    locations,
    "blocks_per_row"
  )
  prep_validate_block_layout_counts(rows_per_block, blocks_per_row, locations)

  list(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    replications_per_repeated_entry = replications_per_repeated_entry,
    n_blocks_per_location = n_blocks_per_location,
    n_candidates = n_candidates,
    field_cols = field_cols,
    rows_per_block = rows_per_block,
    blocks_per_row = blocks_per_row
  )
}

prep_normalize_duplicate_counts <- function(treatments, locations, prep_prop, n_duplicate_per_location) {
  n_treatments <- length(treatments)
  n_locations <- length(locations)

  if (is.null(n_duplicate_per_location)) {
    if (is.null(prep_prop)) {
      out <- prep_balanced_counts(n_treatments, n_locations)
    } else {
      out <- rep(round(n_treatments * prep_prop), n_locations)
    }
  } else if (length(n_duplicate_per_location) == 1) {
    out <- rep(as.integer(n_duplicate_per_location), n_locations)
  } else if (!is.null(names(n_duplicate_per_location))) {
    out <- as.integer(n_duplicate_per_location[locations])
  } else {
    out <- as.integer(n_duplicate_per_location)
  }

  if (length(out) != n_locations || any(is.na(out)) || any(out < 0)) {
    stop("`n_duplicate_per_location` must be non-negative and match locations.", call. = FALSE)
  }
  if (sum(out) > n_treatments) {
    stop("Total duplicated entries cannot exceed number of treatments.", call. = FALSE)
  }

  out
}


prep_normalize_location_counts <- function(x, locations, label) {
  if (length(x) == 1) {
    out <- rep(arc_positive_integer(x, label), length(locations))
    names(out) <- locations
    return(out)
  }

  if (!is.null(names(x))) {
    out <- as.integer(x[locations])
  } else {
    out <- as.integer(x)
  }

  if (length(out) != length(locations) || any(is.na(out)) || any(out < 1)) {
    stop("`", label, "` must contain positive integers for every location.", call. = FALSE)
  }
  names(out) <- locations
  out
}

prep_normalize_optional_location_counts <- function(x, locations, label) {
  if (is.null(x)) {
    return(NULL)
  }
  prep_normalize_location_counts(x, locations, label)
}

prep_normalize_optional_partial_location_counts <- function(x, locations, label) {
  out <- stats::setNames(rep(NA_integer_, length(locations)), locations)
  if (is.null(x)) {
    return(out)
  }

  if (length(x) == 1 && is.null(names(x))) {
    out[] <- arc_positive_integer(x, label)
    return(out)
  }

  if (!is.null(names(x))) {
    unknown <- setdiff(names(x), locations)
    if (length(unknown) > 0) {
      stop("`", label, "` contains unknown location(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    values <- as.integer(x)
    if (any(is.na(values)) || any(values < 1)) {
      stop("`", label, "` must contain positive integers.", call. = FALSE)
    }
    out[names(x)] <- values
    return(out)
  }

  values <- as.integer(x)
  if (length(values) != length(locations) || any(is.na(values)) || any(values < 1)) {
    stop("`", label, "` must be a positive integer, a named vector by location, or match `locations`.", call. = FALSE)
  }
  out[] <- values
  out
}

prep_validate_block_layout_counts <- function(rows_per_block, blocks_per_row, locations) {
  rows_set <- !is.na(rows_per_block)
  blocks_set <- !is.na(blocks_per_row)
  incomplete <- locations[xor(rows_set, blocks_set)]
  if (length(incomplete) > 0) {
    stop(
      "`rows_per_block` and `blocks_per_row` must both be supplied for location(s): ",
      paste(incomplete, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


prep_normalize_field_cols <- function(field_cols, locations) {
  if (is.null(field_cols)) {
    return(NULL)
  }
  if (length(field_cols) == 1) {
    out <- rep(arc_positive_integer(field_cols, "field_cols"), length(locations))
    names(out) <- locations
    return(out)
  }
  if (is.null(names(field_cols))) {
    stop("If `field_cols` has length greater than 1, it must be named by location.", call. = FALSE)
  }
  out <- as.integer(field_cols[locations])
  if (any(is.na(out)) || any(out < 1)) {
    stop("`field_cols` must contain positive integers for every location.", call. = FALSE)
  }
  names(out) <- locations
  out
}

# Fixed, candidate-independent quantities for the candidate search.
prep_candidate_setup <- function(inputs) {
  n_treatments <- length(inputs$treatments)
  locations <- inputs$locations
  reps <- as.integer(inputs$replications_per_repeated_entry[locations])
  n_duplicate <- as.integer(inputs$n_duplicate_per_location)

  block_sizes <- lapply(seq_along(locations), function(l) {
    loc <- locations[l]
    n_plots <- n_treatments + n_duplicate[l] * (reps[l] - 1L)
    n_blocks_here <- if (is.null(inputs$n_blocks_per_location)) NULL else inputs$n_blocks_per_location[[loc]]
    sizes <- prep_balanced_block_sizes(n_plots, inputs$block_size[[loc]], n_blocks_here)
    if (n_duplicate[l] > 0 && reps[l] > length(sizes)) {
      stop(
        "Location `", loc, "` has ", reps[l],
        " total reps for selected entries but only ", length(sizes),
        " blocks. Increase `n_blocks_per_location` or reduce `replications_per_repeated_entry`.",
        call. = FALSE
      )
    }
    sizes
  })

  # Each entry has one plot per location plus (reps - 1) extra plots at the
  # location where it is repeated, so the multiset of replication counts (and
  # its variance term in the score) is the same for every candidate.
  replication_counts <- c(
    rep(length(locations), n_treatments - sum(n_duplicate)),
    rep(length(locations) + reps - 1L, n_duplicate)
  )
  replication_variance <- sum((replication_counts - mean(replication_counts))^2)
  any_pairs <- any(unlist(block_sizes) > 1L)

  list(
    n_treatments = n_treatments,
    n_locations = length(locations),
    reps = reps,
    n_duplicate = n_duplicate,
    block_sizes = block_sizes,
    replication_variance = replication_variance,
    score_lower_bound = replication_variance + if (any_pairs) 1000 else 0
  )
}

# One candidate: random repeated-entry groups, random plot order per location,
# and a greedy block allocation that avoids repeating an entry in a block and
# avoids pairs of entries meeting in more than one block (across locations).
# Works on integer entry ids and keeps the score statistics incrementally.
prep_make_candidate <- function(setup) {
  n <- setup$n_treatments
  n_dup <- setup$n_duplicate
  duplicated_ids <- sample.int(n, sum(n_dup))
  groups <- split(duplicated_ids, factor(rep.int(seq_along(n_dup), n_dup), levels = seq_along(n_dup)))

  concurrence <- matrix(0L, n, n)
  n_pairs_gt1 <- 0
  total_excess <- 0
  max_conc <- 0
  same_entry_same_block <- 0
  allocations <- vector("list", setup$n_locations)

  for (l in seq_len(setup$n_locations)) {
    g <- groups[[l]]
    reps <- setup$reps[l]
    n_extra <- reps - 1L
    occ_entry <- c(seq_len(n), rep.int(g, n_extra))
    occ_number <- c(rep.int(1L, n), rep(seq_len(n_extra) + 1L, each = length(g)))
    ord <- sample.int(length(occ_entry))
    occ_entry <- occ_entry[ord]
    occ_number <- occ_number[ord]

    block_sizes <- setup$block_sizes[[l]]
    n_blocks <- length(block_sizes)
    remaining <- block_sizes
    presence <- matrix(0, nrow = n_blocks, ncol = n)
    members <- vector("list", n_blocks)
    occ_block <- integer(length(occ_entry))

    for (i in seq_along(occ_entry)) {
      e <- occ_entry[i]
      open <- which(remaining > 0L)
      free <- presence[open, e] == 0
      if (any(free)) {
        open <- open[free]
      }

      penalty <- (block_sizes[open] - remaining[open]) / block_sizes[open] + stats::runif(length(open))
      ce <- concurrence[, e]
      met <- which(ce > 0L)
      if (length(met) > 0) {
        # sum over block members m of: 1e4 * [conc >= 1] + 100 * conc^2
        penalty <- penalty + as.vector(presence[open, met, drop = FALSE] %*% (1e4 + 100 * ce[met]^2))
      }
      b <- open[which.min(penalty)]

      if (presence[b, e] == 0) {
        mem <- members[[b]]
        if (length(mem) > 0) {
          new_conc <- concurrence[mem, e] + 1L
          concurrence[mem, e] <- new_conc
          concurrence[e, mem] <- new_conc
          n_pairs_gt1 <- n_pairs_gt1 + sum(new_conc == 2L)
          total_excess <- total_excess + sum(new_conc >= 2L)
          max_conc <- max(max_conc, new_conc)
        }
        presence[b, e] <- 1
        members[[b]] <- c(mem, e)
      } else {
        same_entry_same_block <- same_entry_same_block + 1
      }

      remaining[b] <- remaining[b] - 1L
      occ_block[i] <- b
    }

    allocations[[l]] <- list(entry = occ_entry, occurrence = occ_number, block = occ_block)
  }

  list(
    groups = groups,
    allocations = allocations,
    score = same_entry_same_block * 1e12 +
      n_pairs_gt1 * 1e7 +
      total_excess * 1e5 +
      max_conc * 1000 +
      setup$replication_variance
  )
}

# Plot-level design, concurrence matrix, and summary for the chosen candidate.
prep_build_candidate_design <- function(candidate, inputs, setup, serpentine) {
  treatments <- inputs$treatments
  locations <- inputs$locations

  groups <- stats::setNames(
    lapply(candidate$groups, function(ids) treatments[ids]),
    locations
  )
  repeat_map <- rep("not_repeated", length(treatments))
  for (l in seq_along(locations)) {
    repeat_map[candidate$groups[[l]]] <- locations[l]
  }

  rows <- lapply(seq_along(locations), function(l) {
    loc <- locations[l]
    alloc <- candidate$allocations[[l]]
    block_sizes <- setup$block_sizes[[l]]
    # Plots within a block keep the order in which they were allocated.
    o <- order(alloc$block, seq_along(alloc$block))
    entry <- alloc$entry[o]
    block <- alloc$block[o]
    repeated_at <- repeat_map[entry]
    data.frame(
      all_entries = treatments[entry],
      occurrence = alloc$occurrence[o],
      repeated_at = repeated_at,
      is_extra_repeat = as.integer(alloc$occurrence[o] > 1L),
      entry_type = ifelse(
        repeated_at == loc,
        "p_rep_entry_at_this_location",
        "single_entry_at_this_location"
      ),
      location = loc,
      block = block,
      block_uid = paste0(loc, "_B", sprintf("%03d", block)),
      block_size = block_sizes[block],
      plot_in_block = sequence(tabulate(block, nbins = length(block_sizes))),
      stringsAsFactors = FALSE
    )
  })

  design <- do.call(rbind, rows)
  design <- prep_add_field_coordinates(
    design,
    locations,
    inputs$field_cols,
    inputs$rows_per_block,
    inputs$blocks_per_row,
    serpentine
  )
  design <- prep_finalize_design(design)
  concurrence <- augmented_prep_concurrence(design, treatments)

  list(
    design = design,
    summary = prep_score_design(design, concurrence),
    concurrence_matrix = concurrence,
    repeated_groups = groups
  )
}

prep_add_field_coordinates <- function(design, locations, field_cols,
                                       rows_per_block, blocks_per_row,
                                       serpentine) {
  design <- design[order(match(design$location, locations), design$block, design$plot_in_block), ]
  design$row <- NA_integer_
  design$col <- NA_integer_

  for (loc in locations) {
    idx <- which(design$location == loc)
    if (!is.na(rows_per_block[[loc]]) && !is.na(blocks_per_row[[loc]])) {
      design <- prep_add_block_layout_coordinates(
        design = design,
        idx = idx,
        rows_per_block = rows_per_block[[loc]],
        blocks_per_row = blocks_per_row[[loc]],
        serpentine = serpentine
      )
    } else if (!is.null(field_cols)) {
      plot_number <- seq_along(idx)
      n_cols <- field_cols[[loc]]
      row <- ceiling(plot_number / n_cols)
      col <- ((plot_number - 1L) %% n_cols) + 1L
      if (isTRUE(serpentine)) {
        col <- ifelse(row %% 2 == 0, n_cols - col + 1L, col)
      }
      design$row[idx] <- row
      design$col[idx] <- col
    }
  }

  rownames(design) <- NULL
  design
}

prep_add_block_layout_coordinates <- function(design, idx, rows_per_block,
                                              blocks_per_row, serpentine) {
  location_design <- design[idx, , drop = FALSE]
  block_ids <- unique(location_design$block)
  block_sizes <- vapply(
    block_ids,
    function(block) sum(location_design$block == block),
    integer(1)
  )
  block_width <- ceiling(max(block_sizes) / rows_per_block)

  for (block_index in seq_along(block_ids)) {
    block <- block_ids[block_index]
    block_idx <- idx[design$block[idx] == block]
    plot_number <- seq_along(block_idx)
    block_band <- (block_index - 1L) %/% blocks_per_row
    block_column <- (block_index - 1L) %% blocks_per_row
    local_row <- ceiling(plot_number / block_width)
    local_col <- ((plot_number - 1L) %% block_width) + 1L

    if (isTRUE(serpentine)) {
      local_col <- ifelse(local_row %% 2 == 0, block_width - local_col + 1L, local_col)
    }

    design$row[block_idx] <- (block_band * rows_per_block) + local_row
    design$col[block_idx] <- (block_column * block_width) + local_col
  }

  design
}

prep_finalize_design <- function(design) {
  design |>
    dplyr::arrange(.data$location, .data$block, .data$plot_in_block) |>
    dplyr::mutate(
      plots = dplyr::row_number(),
      rep = .data$occurrence,
      is_control = 0L
    ) |>
    dplyr::select(dplyr::all_of(c(
      "plots",
      "location",
      "row",
      "col",
      "block",
      "block_uid",
      "block_size",
      "plot_in_block",
      "all_entries",
      "occurrence",
      "entry_type",
      "is_extra_repeat",
      "repeated_at",
      "rep",
      "is_control"
    ))) |>
    tibble::as_tibble()
}

prep_score_design <- function(design, concurrence) {
  upper <- concurrence[upper.tri(concurrence)]
  block_entries <- split(design$all_entries, interaction(design$location, design$block, drop = TRUE))
  same_entry_same_block <- sum(vapply(block_entries, function(x) sum(duplicated(x)), numeric(1)))
  replication_counts <- table(design$all_entries)
  max_conc <- if (length(upper) == 0) 0 else max(upper)
  n_pairs_gt1 <- sum(upper > 1)
  total_excess <- sum(pmax(upper - 1, 0))
  score <- same_entry_same_block * 1e12 +
    n_pairs_gt1 * 1e7 +
    total_excess * 1e5 +
    max_conc * 1000 +
    sum((as.numeric(replication_counts) - mean(as.numeric(replication_counts)))^2)

  tibble::tibble(
    score = score,
    max_pairwise_concurrence = max_conc,
    n_pairs_with_concurrence_gt_1 = n_pairs_gt1,
    total_excess_concurrence = total_excess,
    same_entry_same_block = same_entry_same_block,
    min_replication = min(replication_counts),
    max_replication = max(replication_counts),
    mean_replication = mean(replication_counts),
    n_locations = dplyr::n_distinct(design$location),
    n_entries = dplyr::n_distinct(design$all_entries),
    n_plots = nrow(design),
    n_blocks = dplyr::n_distinct(design$block_uid)
  )
}

prep_balanced_counts <- function(total, n_groups) {
  base <- total %/% n_groups
  rem <- total %% n_groups
  out <- rep(base, n_groups)
  if (rem > 0) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  out
}

prep_balanced_block_sizes <- function(n_plots, max_block_size, n_blocks = NULL) {
  if (n_plots <= 0) stop("`n_plots` must be positive.", call. = FALSE)
  if (max_block_size <= 1) stop("`block_size` must be greater than 1.", call. = FALSE)

  if (is.null(n_blocks)) {
    n_blocks <- ceiling(n_plots / max_block_size)
  } else {
    n_blocks <- arc_positive_integer(n_blocks, "n_blocks_per_location")
  }

  if (n_blocks > n_plots) {
    stop("`n_blocks_per_location` cannot exceed the number of plots at a location.", call. = FALSE)
  }

  out <- prep_balanced_counts(n_plots, n_blocks)
  if (max(out) > max_block_size) {
    stop(
      "`n_blocks_per_location` is too small for `block_size`; ",
      "at least one block would contain ", max(out), " plots.",
      call. = FALSE
    )
  }

  out
}

