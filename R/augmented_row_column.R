#' Calculate augmented row-column design capacity
#'
#' @param rows_in_field Number of rows in the whole field.
#' @param cols_in_field Number of columns in the whole field.
#' @param rows_per_block Number of rows in each row-column block.
#' @param cols_per_block Number of columns in each row-column block.
#' @param n_controls Number of replicated checks/controls. Each block receives
#'   one copy of each control.
#'
#' @return A one-row tibble with number of blocks, total plots, check plots, and
#' entry plots available for unreplicated treatments.
#'
#' @examples
#' augmented_row_column_capacity(
#'   rows_in_field = 6,
#'   cols_in_field = 6,
#'   rows_per_block = 3,
#'   cols_per_block = 3,
#'   n_controls = 2
#' )
#'
#' @export
augmented_row_column_capacity <- function(
    rows_in_field,
    cols_in_field,
    rows_per_block,
    cols_per_block,
    n_controls
) {
  rows_in_field <- arc_positive_integer(rows_in_field, "rows_in_field")
  cols_in_field <- arc_positive_integer(cols_in_field, "cols_in_field")
  rows_per_block <- arc_positive_integer(rows_per_block, "rows_per_block")
  cols_per_block <- arc_positive_integer(cols_per_block, "cols_per_block")
  n_controls <- arc_positive_integer(n_controls, "n_controls")

  if (rows_in_field %% rows_per_block != 0) {
    stop("`rows_in_field` must be divisible by `rows_per_block`.", call. = FALSE)
  }

  if (cols_in_field %% cols_per_block != 0) {
    stop("`cols_in_field` must be divisible by `cols_per_block`.", call. = FALSE)
  }

  if (n_controls > rows_per_block || n_controls > cols_per_block) {
    stop(
      "`n_controls` must fit in each block without duplicated check rows or columns.",
      call. = FALSE
    )
  }

  n_block_rows <- rows_in_field / rows_per_block
  n_block_cols <- cols_in_field / cols_per_block
  n_blocks <- n_block_rows * n_block_cols
  total_plots <- rows_in_field * cols_in_field
  check_plots <- n_blocks * n_controls
  entry_plots <- total_plots - check_plots

  tibble::tibble(
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    n_block_rows = n_block_rows,
    n_block_cols = n_block_cols,
    n_blocks = n_blocks,
    n_controls = n_controls,
    total_plots = total_plots,
    check_plots = check_plots,
    entry_plots = entry_plots
  )
}

#' Generate an augmented row-column field design
#'
#' Randomizes unreplicated treatments and replicated controls in an augmented
#' row-column design. Each row-column block receives one copy of each control,
#' and controls are constrained so they do not share a row or column within a
#' block. Candidate layouts are scored and the best candidate is returned.
#'
#' @param treatments Character vector of unreplicated treatment/entry names.
#' @param controls Character vector of replicated check/control names.
#' @param rows_in_field Number of rows in the whole field.
#' @param cols_in_field Number of columns in the whole field.
#' @param rows_per_block Number of rows in each row-column block.
#' @param cols_per_block Number of columns in each row-column block.
#' @param plot_type Plot numbering order. Use `"serpentine"` or `"cartesian"`.
#' @param n_candidates Number of candidate layouts to generate and score.
#' @param seed Optional random seed.
#' @param output_file Optional path for writing the complete design as a
#'   tab-delimited file.
#'
#' @return A list with:
#' \describe{
#'   \item{input}{Validated design inputs.}
#'   \item{summary}{Design size and optimization summary.}
#'   \item{design}{Complete plot-level design with row, column, block, entry,
#'   control flag, and plot order.}
#'   \item{breedbase_design}{Breedbase-compatible columns:
#'   `plots`, `block`, `all_entries`, `rep`, and `is_control`.}
#' }
#'
#' @examples
#' treatments <- paste0("G", seq_len(28))
#' controls <- c("Check1", "Check2")
#'
#' design <- augmented_row_column_design(
#'   treatments = treatments,
#'   controls = controls,
#'   rows_in_field = 6,
#'   cols_in_field = 6,
#'   rows_per_block = 3,
#'   cols_per_block = 3,
#'   n_candidates = 5,
#'   seed = 123
#' )
#'
#' head(design$design)
#' design$summary
#'
#' @export
augmented_row_column_design <- function(
    treatments,
    controls,
    rows_in_field,
    cols_in_field,
    rows_per_block,
    cols_per_block,
    plot_type = c("serpentine", "cartesian"),
    n_candidates = 1000,
    seed = NULL,
    output_file = NULL
) {
  plot_type <- match.arg(plot_type)
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

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

  inputs <- arc_validate_inputs(
    treatments = treatments,
    controls = controls,
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    n_candidates = n_candidates,
    plot_type = plot_type
  )

  field_template <- arc_make_field_template(
    rows_in_field = inputs$rows_in_field,
    cols_in_field = inputs$cols_in_field,
    rows_per_block = inputs$rows_per_block,
    cols_per_block = inputs$cols_per_block
  )

  best_design <- NULL
  best_score <- Inf
  n_valid_candidates <- 0L

  for (i in seq_len(inputs$n_candidates)) {
    candidate <- try(
      arc_allocate_augmented_row_column(
        field_template = field_template,
        treatments = inputs$treatments,
        controls = inputs$controls,
        plot_type = inputs$plot_type
      ),
      silent = TRUE
    )

    if (inherits(candidate, "try-error")) {
      next
    }

    n_valid_candidates <- n_valid_candidates + 1L
    candidate_score <- arc_score_augmented_design(candidate)

    if (candidate_score < best_score) {
      best_score <- candidate_score
      best_design <- candidate
    }
  }

  if (is.null(best_design)) {
    stop(
      "Unable to generate a valid augmented row-column design. ",
      "Try increasing `n_candidates` or changing the field/block dimensions.",
      call. = FALSE
    )
  }

  capacity <- augmented_row_column_capacity(
    rows_in_field = inputs$rows_in_field,
    cols_in_field = inputs$cols_in_field,
    rows_per_block = inputs$rows_per_block,
    cols_per_block = inputs$cols_per_block,
    n_controls = length(inputs$controls)
  )

  design <- best_design |>
    dplyr::transmute(
      plots = .data$plots,
      row = .data$row,
      col = .data$col,
      block = .data$block,
      rowgroup = .data$rowgroup,
      colgroup = .data$colgroup,
      all_entries = .data$trt,
      type = .data$type,
      rep = .data$rep,
      is_control = .data$is_control
    )

  breedbase_design <- design |>
    dplyr::select(dplyr::all_of(c("plots", "block", "all_entries", "rep", "is_control")))

  summary <- capacity |>
    dplyr::mutate(
      plot_type = inputs$plot_type,
      n_candidates = inputs$n_candidates,
      n_valid_candidates = n_valid_candidates,
      design_score = best_score
    )

  result <- list(
    input = list(
      treatments = inputs$treatments,
      controls = inputs$controls,
      rows_in_field = inputs$rows_in_field,
      cols_in_field = inputs$cols_in_field,
      rows_per_block = inputs$rows_per_block,
      cols_per_block = inputs$cols_per_block,
      plot_type = inputs$plot_type,
      n_candidates = inputs$n_candidates,
      seed = seed
    ),
    summary = summary,
    design = design,
    breedbase_design = breedbase_design
  )

  class(result) <- c("augmented_row_column_design", class(result))

  if (!is.null(output_file)) {
    utils::write.table(
      design,
      file = output_file,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
  }

  result
}

#' Plot an augmented row-column design grid
#'
#' Draws a field grid for an augmented row-column design, with plot cells,
#' row and column labels, and thick outlines around row-column blocks.
#'
#' @param design An object returned by `augmented_row_column_design()` or a
#'   design data frame with row, column, and block columns.
#' @param label Cell label to show. Use `"all_entries"`, `"type"`, `"block"`,
#'   or `"none"`.
#' @param fill Cell coloring variable. Use `"type"`, `"block"`, `"rowgroup"`,
#'   `"colgroup"`, or `"is_control"`.
#' @param show_plot_labels Logical. If `TRUE`, labels are printed inside cells
#'   unless `label = "none"`.
#' @param show_block_labels Logical. If `TRUE`, block numbers are printed at
#'   the center of each block.
#' @param label_cex Optional text size for plot-cell labels.
#' @param block_label_cex Text size for block labels.
#' @param block_lwd Line width for block boundaries.
#' @param main Optional plot title.
#' @param legend Logical. If `TRUE`, draw a legend for cell colors.
#' @param ... Additional arguments passed to `graphics::plot.window()`.
#'
#' @return Invisibly returns the plotted design data frame with internal
#'   plotting coordinates.
#'
#' @examples
#' design <- augmented_row_column_design(
#'   treatments = paste0("G", seq_len(28)),
#'   controls = c("Check1", "Check2"),
#'   rows_in_field = 6,
#'   cols_in_field = 6,
#'   rows_per_block = 3,
#'   cols_per_block = 3,
#'   n_candidates = 3,
#'   seed = 123
#' )
#'
#' plot_augmented_row_column_design(design, label = "all_entries")
#'
#' @export
plot_augmented_row_column_design <- function(
    design,
    label = c("all_entries", "type", "block", "none"),
    fill = c("type", "block", "rowgroup", "colgroup", "is_control"),
    show_plot_labels = TRUE,
    show_block_labels = TRUE,
    label_cex = NULL,
    block_label_cex = 0.9,
    block_lwd = 2,
    main = NULL,
    legend = TRUE,
    ...
) {
  label <- match.arg(label)
  fill <- match.arg(fill)

  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) {
    design <- design$design
  }

  if (!is.data.frame(design)) {
    stop("`design` must be a design data frame or an augmented row-column design object.", call. = FALSE)
  }

  required_cols <- c("row", "col", "block")
  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (label != "none" && !label %in% names(design)) {
    stop("Column `", label, "` was not found in `design`.", call. = FALSE)
  }

  if (!fill %in% names(design)) {
    fill <- "block"
  }

  design <- as.data.frame(design, stringsAsFactors = FALSE)
  design$row <- suppressWarnings(as.numeric(design$row))
  design$col <- suppressWarnings(as.numeric(design$col))

  if (any(is.na(design$row)) || any(is.na(design$col))) {
    stop("`row` and `col` columns must be numeric or coercible to numeric.", call. = FALSE)
  }

  rows <- sort(unique(design$row))
  cols <- sort(unique(design$col))
  design$.plot_row <- match(design$row, rows)
  design$.plot_col <- match(design$col, cols)
  design$.fill_value <- as.character(design[[fill]])
  fill_levels <- sort(unique(design$.fill_value))
  fill_colors <- grDevices::hcl.colors(length(fill_levels), palette = "Set 3")
  names(fill_colors) <- fill_levels

  if (is.null(label_cex)) {
    label_cex <- max(0.25, min(0.8, 7 / max(length(rows), length(cols))))
  }

  if (is.null(main)) {
    main <- "Augmented row-column design"
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 4, if (isTRUE(legend)) 7 else 2), xpd = NA)

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0.5, length(cols) + 0.5),
    ylim = c(length(rows) + 0.5, 0.5),
    xaxs = "i",
    yaxs = "i",
    asp = 1,
    ...
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

  block_split <- split(design, design$block)
  for (block_design in block_split) {
    min_col <- min(block_design$.plot_col)
    max_col <- max(block_design$.plot_col)
    min_row <- min(block_design$.plot_row)
    max_row <- max(block_design$.plot_row)

    graphics::rect(
      xleft = min_col - 0.5,
      ybottom = min_row - 0.5,
      xright = max_col + 0.5,
      ytop = max_row + 0.5,
      border = "black",
      lwd = block_lwd
    )

    if (isTRUE(show_block_labels)) {
      graphics::text(
        x = mean(c(min_col, max_col)),
        y = mean(c(min_row, max_row)),
        labels = paste0("B", unique(block_design$block)[1]),
        cex = block_label_cex,
        font = 2,
        col = "black"
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
      legend = fill_levels,
      fill = fill_colors,
      border = "grey60",
      title = fill,
      bty = "n",
      xjust = 0,
      yjust = 0,
      cex = 0.8
    )
  }

  invisible(tibble::as_tibble(design))
}

#' Analyze an augmented row-column design
#'
#' Fits the fixed-effects model:
#' `y = mean + rowgroup + colgroup + rowgroup:colgroup + rowgroup:row +
#' colgroup:col + genotype + error`.
#' The formatted ANOVA table reports degrees of freedom from augmented
#' row-column design formulas using field rows (`k`), field columns (`s`), row
#' groups (`gk`), column groups (`gs`), check cultivars (`vc`), and
#' unreplicated entries (`ve`).
#'
#' @param data A complete augmented row-column design data frame, usually
#'   `design$design` from `augmented_row_column_design()`, with a response
#'   column added.
#' @param response Name of the numeric response column.
#' @param genotype_col Name of the genotype/treatment column.
#' @param rowgroup_col Name of the row-group column.
#' @param colgroup_col Name of the column-group column.
#' @param row_col Name of the field row column.
#' @param col_col Name of the field column column.
#' @param type_col Optional entry type column. Used only for genotype summaries.
#' @param control_col Optional control flag column. Used only for genotype
#'   summaries.
#'
#' @return A list with the analysis input, fitted model, model formula, ANOVA
#' table, formatted ANOVA table, plot-level fitted values and residuals,
#' genotype summaries, and model diagnostics.
#'
#' @examples
#' design <- augmented_row_column_design(
#'   treatments = paste0("G", seq_len(108)),
#'   controls = paste0("Check", seq_len(4)),
#'   rows_in_field = 12,
#'   cols_in_field = 12,
#'   rows_per_block = 4,
#'   cols_per_block = 4,
#'   n_candidates = 3,
#'   seed = 123
#' )
#'
#' design_data <- design$design
#' design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
#' analysis <- analyze_augmented_row_column_design(design_data, response = "yield")
#' analysis$anova
#'
#' @export
analyze_augmented_row_column_design <- function(
    data,
    response,
    genotype_col = "all_entries",
    rowgroup_col = "rowgroup",
    colgroup_col = "colgroup",
    row_col = "row",
    col_col = "col",
    type_col = "type",
    control_col = "is_control"
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  required_cols <- c(response, genotype_col, rowgroup_col, colgroup_col, row_col, col_col)
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
      rowgroup = as.factor(.data[[rowgroup_col]]),
      colgroup = as.factor(.data[[colgroup_col]]),
      row = as.factor(.data[[row_col]]),
      col = as.factor(.data[[col_col]]),
      row_nested = interaction(.data[[rowgroup_col]], .data[[row_col]], drop = TRUE),
      col_nested = interaction(.data[[colgroup_col]], .data[[col_col]], drop = TRUE),
      type = if (type_col %in% names(data)) {
        as.character(.data[[type_col]])
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
      !is.na(.data$rowgroup),
      !is.na(.data$colgroup),
      !is.na(.data$row),
      !is.na(.data$col)
    )

  analysis_data <- tibble::as_tibble(data) |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::transmute(
      .row_id = .data$.row_id,
      response = suppressWarnings(as.numeric(.data[[response]])),
      genotype = as.factor(.data[[genotype_col]]),
      rowgroup = as.factor(.data[[rowgroup_col]]),
      colgroup = as.factor(.data[[colgroup_col]]),
      row = as.factor(.data[[row_col]]),
      col = as.factor(.data[[col_col]]),
      row_nested = interaction(.data[[rowgroup_col]], .data[[row_col]], drop = TRUE),
      col_nested = interaction(.data[[colgroup_col]], .data[[col_col]], drop = TRUE),
      type = if (type_col %in% names(data)) {
        as.character(.data[[type_col]])
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
      !is.na(.data$rowgroup),
      !is.na(.data$colgroup),
      !is.na(.data$row),
      !is.na(.data$col)
    )

  if (nrow(analysis_data) == 0) {
    stop("No complete observations are available for analysis.", call. = FALSE)
  }

  if (dplyr::n_distinct(analysis_data$genotype) < 2) {
    stop("At least two genotypes are required for analysis.", call. = FALSE)
  }

  model_formula <- stats::as.formula(
    "response ~ rowgroup + row_nested + colgroup + col_nested + rowgroup:colgroup + genotype"
  )

  design_df <- arc_augmented_design_df(design_structure_data)
  fit <- arc_lm_drop_aliased_terms(model_formula, analysis_data, design_df = design_df)
  anova_result <- arc_lm_anova_table(fit)
  anova_table <- anova_result$table
  anova_formatted <- arc_format_augmented_anova(
    anova_table = anova_table,
    response = analysis_data$response,
    design_df = design_df
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
      type = if (type_col %in% names(plot_results)) {
        dplyr::first(as.character(.data[[type_col]]))
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(plot_results)) {
        dplyr::first(as.integer(.data[[control_col]]))
      } else {
        NA_integer_
      },
      n_observations = sum(!is.na(.data$.analysis_response)),
      raw_mean = mean(.data$.analysis_response, na.rm = TRUE),
      fitted_mean = mean(.data$fitted, na.rm = TRUE),
      residual_mean = mean(.data$residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::select(-dplyr::all_of(genotype_col))

  genotype_summary$raw_mean[is.nan(genotype_summary$raw_mean)] <- NA_real_
  genotype_summary$fitted_mean[is.nan(genotype_summary$fitted_mean)] <- NA_real_
  genotype_summary$residual_mean[is.nan(genotype_summary$residual_mean)] <- NA_real_

  diagnostics <- tibble::tibble(
    n_observations = stats::nobs(fit),
    n_genotypes = dplyr::n_distinct(analysis_data$genotype),
    model_rank = fit$rank,
    residual_df = stats::df.residual(fit),
    anova_error_df = design_df$display_df[design_df$term == "Residuals"],
    aliased_coefficients_dropped = attr(fit, "aliased_coefficients_dropped") %||% 0L,
    sigma = summary(fit)$sigma,
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    rank_deficient = fit$rank < length(stats::coef(fit)),
    saturated_or_near_perfect_fit = isTRUE(stats::df.residual(fit) == 0) ||
      isTRUE(!is.na(summary(fit)$sigma) && summary(fit)$sigma < sqrt(.Machine$double.eps)),
    anova_warning = anova_result$warning %||% NA_character_
  )

  analysis_message <- arc_augmented_analysis_message(diagnostics)

  result <- list(
    input = list(
      response = response,
      genotype_col = genotype_col,
      rowgroup_col = rowgroup_col,
      colgroup_col = colgroup_col,
      row_col = row_col,
      col_col = col_col,
      type_col = type_col,
      control_col = control_col
    ),
    formula = model_formula,
    model = fit,
    anova = anova_table,
    anova_df = design_df,
    anova_formatted = anova_formatted,
    plot_results = plot_results |> dplyr::select(-dplyr::all_of(".analysis_response")),
    genotype_summary = genotype_summary,
    diagnostics = diagnostics,
    message = analysis_message
  )

  class(result) <- c("augmented_row_column_analysis", class(result))
  result
}

#' Print augmented row-column analysis
#'
#' @param x An object returned by `analyze_augmented_row_column_design()`.
#' @param ... Additional arguments passed to `print()`.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.augmented_row_column_analysis <- function(x, ...) {
  cat("Augmented row-column design analysis\n\n")
  print(x$anova_formatted, ...)
  if (!is.null(x$message) && !is.na(x$message) && nzchar(x$message)) {
    cat("\n", x$message, "\n", sep = "")
  }
  invisible(x)
}

arc_positive_integer <- function(x, label) {
  if (length(x) != 1 || is.na(x)) {
    stop("`", label, "` must be a single positive integer.", call. = FALSE)
  }

  x <- as.integer(x)

  if (is.na(x) || x < 1) {
    stop("`", label, "` must be a single positive integer.", call. = FALSE)
  }

  x
}

arc_validate_inputs <- function(
    treatments,
    controls,
    rows_in_field,
    cols_in_field,
    rows_per_block,
    cols_per_block,
    n_candidates,
    plot_type
) {
  treatments <- as.character(treatments)
  controls <- as.character(controls)

  treatments <- treatments[!is.na(treatments) & treatments != ""]
  controls <- controls[!is.na(controls) & controls != ""]

  rows_in_field <- arc_positive_integer(rows_in_field, "rows_in_field")
  cols_in_field <- arc_positive_integer(cols_in_field, "cols_in_field")
  rows_per_block <- arc_positive_integer(rows_per_block, "rows_per_block")
  cols_per_block <- arc_positive_integer(cols_per_block, "cols_per_block")
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

  if (length(controls) == 0) {
    stop("At least one control/check is required.", call. = FALSE)
  }

  if (length(treatments) == 0) {
    stop("At least one treatment/entry is required.", call. = FALSE)
  }

  if (anyDuplicated(controls)) {
    stop("`controls` contains duplicated names.", call. = FALSE)
  }

  if (anyDuplicated(treatments)) {
    stop("`treatments` contains duplicated names.", call. = FALSE)
  }

  overlap <- intersect(treatments, controls)
  if (length(overlap) > 0) {
    stop(
      "`treatments` and `controls` must not contain the same name: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }

  capacity <- augmented_row_column_capacity(
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    n_controls = length(controls)
  )

  expected_entries <- capacity$entry_plots

  if (length(treatments) != expected_entries) {
    stop(
      "The number of treatments does not match the available entry plots. ",
      "Expected ",
      expected_entries,
      " treatments, but received ",
      length(treatments),
      ".",
      call. = FALSE
    )
  }

  n_block_rows <- rows_in_field / rows_per_block
  n_block_cols <- cols_in_field / cols_per_block
  if (n_block_cols > rows_per_block) {
    stop(
      "`rows_per_block` must be at least the number of column groups ",
      "to avoid repeating the same control in a field row.",
      call. = FALSE
    )
  }
  if (n_block_rows > cols_per_block) {
    stop(
      "`cols_per_block` must be at least the number of row groups ",
      "to avoid repeating the same control in a field column.",
      call. = FALSE
    )
  }

  list(
    treatments = treatments,
    controls = controls,
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    n_candidates = n_candidates,
    plot_type = plot_type
  )
}

arc_make_field_template <- function(rows_in_field, cols_in_field, rows_per_block, cols_per_block) {
  block_cols <- cols_in_field / cols_per_block

  field <- expand.grid(
    row = seq_len(rows_in_field),
    col = seq_len(cols_in_field),
    KEEP.OUT.ATTRS = FALSE
  )

  field$rowgroup <- ((field$row - 1L) %/% rows_per_block) + 1L
  field$colgroup <- ((field$col - 1L) %/% cols_per_block) + 1L
  field$block <- ((field$rowgroup - 1L) * block_cols) + field$colgroup
  field$trt <- NA_character_
  field$type <- "entry"
  field
}

arc_lm_anova_table <- function(fit) {
  if (!is.null(fit$x) && !is.null(fit$y)) {
    out <- arc_sequential_anova_from_fit(fit)
    warning_text <- NULL
  } else {
    warning_text <- NULL
    tab <- withCallingHandlers(
      stats::anova(fit),
      warning = function(w) {
        warning_text <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )

    out <- tibble::tibble(
      term = rownames(tab),
      df = as.numeric(tab[["Df"]]),
      sum_sq = as.numeric(tab[["Sum Sq"]]),
      mean_sq = as.numeric(tab[["Mean Sq"]]),
      f_value = as.numeric(tab[["F value"]]),
      p_value = as.numeric(tab[["Pr(>F)"]])
    )
    rownames(out) <- NULL
  }

  out$term[out$term == "row_nested"] <- "rowgroup:row"
  out$term[out$term == "col_nested"] <- "colgroup:col"
  out <- arc_clean_anova_numbers(out)

  list(
    table = out,
    warning = warning_text
  )
}

arc_sequential_anova_from_fit <- function(fit) {
  x <- fit$x
  y <- fit$y
  assign <- fit$assign
  term_labels <- attr(fit$terms, "term.labels")
  current_cols <- which(assign == 0)

  current_x <- x[, current_cols, drop = FALSE]
  current_rank <- qr(current_x)$rank
  current_rss <- sum(stats::lm.fit(current_x, y)$residuals^2)
  rows <- vector("list", length(term_labels))

  for (i in seq_along(term_labels)) {
    term_cols <- which(assign == i)

    if (length(term_cols) == 0) {
      rows[[i]] <- tibble::tibble(
        term = term_labels[i],
        df = 0,
        sum_sq = 0,
        mean_sq = NA_real_,
        f_value = NA_real_,
        p_value = NA_real_
      )
      next
    }

    candidate_x <- cbind(current_x, x[, term_cols, drop = FALSE])
    candidate_fit <- stats::lm.fit(candidate_x, y)
    candidate_rank <- qr(candidate_x)$rank
    candidate_rss <- sum(candidate_fit$residuals^2)
    term_df <- candidate_rank - current_rank
    term_ss <- current_rss - candidate_rss

    rows[[i]] <- tibble::tibble(
      term = term_labels[i],
      df = term_df,
      sum_sq = term_ss,
      mean_sq = if (term_df > 0) term_ss / term_df else NA_real_,
      f_value = NA_real_,
      p_value = NA_real_
    )

    current_x <- candidate_x
    current_rank <- candidate_rank
    current_rss <- candidate_rss
  }

  residual_df <- length(y) - current_rank
  residual_mean_sq <- if (residual_df > 0) current_rss / residual_df else NA_real_
  out <- dplyr::bind_rows(rows)

  if (!is.na(residual_mean_sq) && residual_mean_sq > 0) {
    out$f_value <- out$mean_sq / residual_mean_sq
    out$p_value <- stats::pf(out$f_value, df1 = out$df, df2 = residual_df, lower.tail = FALSE)
    out$f_value[out$df <= 0] <- NA_real_
    out$p_value[out$df <= 0] <- NA_real_
  }

  dplyr::bind_rows(
    out,
    tibble::tibble(
      term = "Residuals",
      df = residual_df,
      sum_sq = current_rss,
      mean_sq = residual_mean_sq,
      f_value = NA_real_,
      p_value = NA_real_
    )
  )
}

arc_lm_drop_aliased_terms <- function(formula, data, design_df = NULL) {
  model_frame <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  response <- stats::model.response(model_frame)
  terms_obj <- stats::terms(formula, data = data, keep.order = TRUE)
  model_matrix <- stats::model.matrix(terms_obj, model_frame)
  full_model_matrix <- model_matrix
  full_assign <- attr(full_model_matrix, "assign")
  term_labels <- attr(terms_obj, "term.labels")
  full_term_names <- ifelse(full_assign == 0, "(Intercept)", term_labels[full_assign])

  if (!is.null(design_df)) {
    selected <- arc_select_augmented_model_columns(
      model_matrix = full_model_matrix,
      assign = full_assign,
      term_labels = term_labels,
      design_df = design_df
    )
  } else {
    qr_obj <- qr(full_model_matrix)
    selected <- sort(qr_obj$pivot[seq_len(qr_obj$rank)])
  }

  model_matrix <- full_model_matrix[, selected, drop = FALSE]
  dropped <- setdiff(seq_len(ncol(full_model_matrix)), selected)

  fit <- stats::lm.fit(x = model_matrix, y = response)
  fit$terms <- terms_obj
  fit$model <- model_frame
  fit$x <- model_matrix
  fit$y <- response
  fit$call <- match.call()
  fit$assign <- full_assign[selected]
  fit$contrasts <- attr(full_model_matrix, "contrasts")
  fit$xlevels <- stats::.getXlevels(terms_obj, model_frame)
  fit$formula <- formula
  fit$na.action <- attr(model_frame, "na.action")
  class(fit) <- "lm"
  attr(fit, "aliased_coefficients_dropped") <- length(dropped)
  attr(fit, "dropped_coefficients") <- colnames(full_model_matrix)[dropped]

  fit
}

arc_select_augmented_model_columns <- function(model_matrix, assign, term_labels, design_df) {
  term_for_column <- ifelse(assign == 0, "(Intercept)", term_labels[assign])
  df_lookup <- stats::setNames(design_df$display_df, design_df$term)
  term_limits <- c(
    rowgroup = df_lookup[["rowgroup"]],
    row_nested = df_lookup[["rowgroup:row"]],
    colgroup = df_lookup[["colgroup"]],
    col_nested = df_lookup[["colgroup:col"]],
    `rowgroup:colgroup` = df_lookup[["rowgroup:colgroup"]],
    genotype = df_lookup[["genotype"]]
  )

  selected <- which(term_for_column == "(Intercept)")
  selected <- selected[seq_len(min(length(selected), 1L))]
  current <- model_matrix[, selected, drop = FALSE]
  current_rank <- qr(current)$rank

  for (term in names(term_limits)) {
    limit <- term_limits[[term]]
    if (is.na(limit) || limit <= 0) {
      next
    }

    candidates <- which(term_for_column == term)
    kept_for_term <- 0L

    for (candidate in candidates) {
      if (kept_for_term >= limit) {
        break
      }

      candidate_matrix <- cbind(current, model_matrix[, candidate, drop = FALSE])
      candidate_rank <- qr(candidate_matrix)$rank

      if (candidate_rank > current_rank) {
        selected <- c(selected, candidate)
        current <- candidate_matrix
        current_rank <- candidate_rank
        kept_for_term <- kept_for_term + 1L
      }
    }
  }

  selected
}

arc_format_augmented_anova <- function(anova_table, response, design_df) {
  source_map <- data.frame(
    term = c(
      "rowgroup",
      "rowgroup:row",
      "colgroup",
      "colgroup:col",
      "rowgroup:colgroup",
      "genotype",
      "Residuals"
    ),
    source = c(
      "Row groups",
      "Rows, nested within row groups",
      "Column groups",
      "Columns, nested within column groups",
      "Row groups x column groups (blocks)",
      "Genotypes",
      "Error"
    ),
    order = seq_len(7),
    stringsAsFactors = FALSE
  )

  out <- source_map |>
    dplyr::left_join(anova_table, by = "term") |>
    dplyr::left_join(design_df, by = "term") |>
    dplyr::arrange(.data$order) |>
    dplyr::mutate(
      df = .data$display_df,
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
    Df = design_df$corrected_total_df[1],
    `Sum Sq` = sum((response - mean(response, na.rm = TRUE))^2, na.rm = TRUE),
    `Mean Sq` = NA_real_,
    `F value` = NA_real_,
    `Pr(>F)` = NA_real_,
    Significance = ""
  )

  dplyr::bind_rows(out, corrected_total)
}

arc_augmented_design_df <- function(data) {
  k <- dplyr::n_distinct(data$row)
  s <- dplyr::n_distinct(data$col)
  gk <- dplyr::n_distinct(data$rowgroup)
  gs <- dplyr::n_distinct(data$colgroup)

  if ("is_control" %in% names(data) && any(!is.na(data$is_control))) {
    vc <- dplyr::n_distinct(data$genotype[data$is_control == 1])
    ve <- dplyr::n_distinct(data$genotype[data$is_control != 1 | is.na(data$is_control)])
  } else if ("type" %in% names(data) && any(!is.na(data$type))) {
    vc <- dplyr::n_distinct(data$genotype[data$type == "check"])
    ve <- dplyr::n_distinct(data$genotype[data$type != "check" | is.na(data$type)])
  } else {
    vc <- 0
    ve <- dplyr::n_distinct(data$genotype)
  }

  genotype_df <- vc + ve - 1
  block_df <- (gk - 1) * (gs - 1)
  corrected_total_df <- k * s - 1
  source_df <- c(
    gk - 1,
    k - gk - 1,
    gs - 1,
    s - gs - 1,
    block_df,
    genotype_df
  )
  error_df <- corrected_total_df - sum(source_df)

  tibble::tibble(
    term = c(
      "rowgroup",
      "rowgroup:row",
      "colgroup",
      "colgroup:col",
      "rowgroup:colgroup",
      "genotype",
      "Residuals"
    ),
    display_df = c(
      source_df,
      error_df
    ),
    k = k,
    s = s,
    gk = gk,
    gs = gs,
    vc = vc,
    ve = ve,
    corrected_total_df = corrected_total_df,
    error_df_formula = paste0(
      "(",
      k,
      " * ",
      s,
      " - 1) - [(",
      gk,
      " - 1) + (",
      k,
      " - ",
      gk,
      " - 1) + (",
      gs,
      " - 1) + (",
      s,
      " - ",
      gs,
      " - 1) + (",
      gk,
      " - 1)(",
      gs,
      " - 1) + (",
      vc,
      " + ",
      ve,
      " - 1)]"
    ),
    previous_error_df_formula = paste0(
      "(",
      k,
      " - 1)(",
      s,
      " - 1) - (",
      gk,
      " - 1)(",
      gs,
      " - 1) - (",
      vc,
      " + ",
      ve,
      " - 1)"
    )
  )
}

arc_significance_stars <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

arc_clean_anova_numbers <- function(x) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], function(col) {
    col[is.nan(col)] <- NA_real_
    col
  })
  x
}

arc_augmented_analysis_message <- function(diagnostics) {
  if (isTRUE(diagnostics$anova_error_df <= 0)) {
    return(
      "F-tests are not estimable because the formatted ANOVA has no positive error degrees of freedom. The design/model is saturated for these data."
    )
  }

  if (isTRUE(diagnostics$saturated_or_near_perfect_fit)) {
    return(
      "F-tests should be interpreted cautiously because the model has an essentially perfect fit."
    )
  }

  NA_character_
}

arc_make_plot_order <- function(df, plot_type = "serpentine") {
  if (plot_type == "serpentine") {
    df$plot_order_col <- ifelse(df$row %% 2 == 1, df$col, -df$col)
    df <- df[order(df$row, df$plot_order_col), ]
    df$plot_order_col <- NULL
  } else {
    df <- df[order(df$row, df$col), ]
  }

  rownames(df) <- NULL
  df$plots <- seq_len(nrow(df))
  df
}

arc_validate_check_layout <- function(df) {
  check_df <- df[df$type == "check", c("block", "row", "col", "trt")]

  for (b in sort(unique(check_df$block))) {
    x <- check_df[check_df$block == b, ]

    if (anyDuplicated(x$row)) {
      stop("Invalid layout: block ", b, " has more than one check in the same row.", call. = FALSE)
    }

    if (anyDuplicated(x$col)) {
      stop("Invalid layout: block ", b, " has more than one check in the same column.", call. = FALSE)
    }

    if (anyDuplicated(x$trt)) {
      stop("Invalid layout: block ", b, " has duplicated check labels.", call. = FALSE)
    }
  }

  if (any(table(check_df$trt, check_df$row) > 1)) {
    stop("Invalid layout: the same check appears more than once in the same field row.", call. = FALSE)
  }

  if (any(table(check_df$trt, check_df$col) > 1)) {
    stop("Invalid layout: the same check appears more than once in the same field column.", call. = FALSE)
  }

  invisible(TRUE)
}

arc_control_permutations <- function(controls) {
  if (length(controls) == 1) {
    return(list(controls))
  }

  perms <- list()
  for (i in seq_along(controls)) {
    rest_perms <- arc_control_permutations(controls[-i])
    for (p in rest_perms) {
      perms[[length(perms) + 1L]] <- c(controls[i], p)
    }
  }

  perms
}

arc_assign_check_names <- function(check_slots, controls) {
  check_slots <- check_slots[order(check_slots$block, check_slots$row, check_slots$col), ]
  check_slots$trt <- NA_character_

  blocks <- split(seq_len(nrow(check_slots)), check_slots$block)
  control_perms <- arc_control_permutations(controls)
  used_rows <- stats::setNames(vector("list", length(controls)), controls)
  used_cols <- stats::setNames(vector("list", length(controls)), controls)

  assign_block <- function(block_number) {
    if (block_number > length(blocks)) {
      return(TRUE)
    }

    idx <- blocks[[block_number]]
    perms <- sample(control_perms)

    for (perm in perms) {
      ok <- TRUE

      for (i in seq_along(idx)) {
        control <- perm[i]

        if (check_slots$row[idx[i]] %in% used_rows[[control]] ||
          check_slots$col[idx[i]] %in% used_cols[[control]]) {
          ok <- FALSE
          break
        }
      }

      if (!ok) {
        next
      }

      old_rows <- used_rows
      old_cols <- used_cols

      for (i in seq_along(idx)) {
        control <- perm[i]
        check_slots$trt[idx[i]] <<- control
        used_rows[[control]] <<- c(used_rows[[control]], check_slots$row[idx[i]])
        used_cols[[control]] <<- c(used_cols[[control]], check_slots$col[idx[i]])
      }

      if (assign_block(block_number + 1L)) {
        return(TRUE)
      }

      check_slots$trt[idx] <<- NA_character_
      used_rows <<- old_rows
      used_cols <<- old_cols
    }

    FALSE
  }

  if (!assign_block(1L)) {
    stop("Could not assign check names without row or column repeats.", call. = FALSE)
  }

  check_slots
}

arc_make_structured_check_slots <- function(field, controls) {
  n_controls <- length(controls)
  rows_per_block <- max(vapply(split(field$row, field$block), function(x) length(unique(x)), integer(1)))
  cols_per_block <- max(vapply(split(field$col, field$block), function(x) length(unique(x)), integer(1)))
  n_block_rows <- dplyr::n_distinct(field$rowgroup)
  n_block_cols <- dplyr::n_distinct(field$colgroup)
  row_offsets <- sample(seq_len(rows_per_block), n_block_rows, replace = TRUE) - 1L
  col_offsets <- sample(seq_len(cols_per_block), n_block_cols, replace = TRUE) - 1L
  control_order <- sample(controls)
  field_key <- paste(field$row, field$col)
  rows <- vector("list", n_block_rows * n_block_cols * n_controls)
  counter <- 1L

  for (rowgroup in seq_len(n_block_rows)) {
    for (colgroup in seq_len(n_block_cols)) {
      block <- unique(field$block[field$rowgroup == rowgroup & field$colgroup == colgroup])

      for (control_idx in seq_len(n_controls)) {
        local_row <- ((control_idx + colgroup + row_offsets[rowgroup] - 2L) %% rows_per_block) + 1L
        local_col <- ((control_idx + rowgroup + col_offsets[colgroup] - 2L) %% cols_per_block) + 1L
        row <- ((rowgroup - 1L) * rows_per_block) + local_row
        col <- ((colgroup - 1L) * cols_per_block) + local_col

        rows[[counter]] <- data.frame(
          check_idx = match(paste(row, col), field_key),
          block = block,
          row = row,
          col = col,
          trt = control_order[control_idx],
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }

  dplyr::bind_rows(rows)
}

arc_score_count_balance <- function(x, weight) {
  if (length(x) == 0) {
    return(0)
  }

  weight * sum((as.numeric(x) - mean(as.numeric(x)))^2)
}

arc_score_augmented_design <- function(df) {
  check_df <- df[df$type == "check", ]

  if (nrow(check_df) == 0) {
    return(Inf)
  }

  score <- 0
  coords <- as.matrix(check_df[, c("row", "col")])
  distances <- as.matrix(stats::dist(coords, method = "manhattan"))
  distances <- distances[upper.tri(distances)]

  if (length(distances) > 0) {
    score <- score + sum(distances == 0) * 100000
    score <- score + sum(distances == 1) * 500
    score <- score + sum(distances == 2) * 100
    score <- score + sum(1 / distances[distances > 0])
  }

  score <- score + arc_score_count_balance(table(check_df$row), 50)
  score <- score + arc_score_count_balance(table(check_df$col), 50)
  score <- score + arc_score_count_balance(table(check_df$trt, check_df$rowgroup), 20)
  score <- score + arc_score_count_balance(table(check_df$trt, check_df$colgroup), 20)

  same_check_row <- table(check_df$trt, check_df$row)
  same_check_col <- table(check_df$trt, check_df$col)
  score <- score + sum(pmax(same_check_row - 1, 0)^2) * 100
  score <- score + sum(pmax(same_check_col - 1, 0)^2) * 100

  score
}

arc_allocate_augmented_row_column <- function(field_template, treatments, controls, plot_type) {
  field <- field_template
  check_slots <- arc_make_structured_check_slots(field, controls)
  field$trt[check_slots$check_idx] <- check_slots$trt
  field$type[check_slots$check_idx] <- "check"

  entry_idx <- which(field$type == "entry")
  field$trt[entry_idx] <- sample(treatments, length(entry_idx), replace = FALSE)

  arc_validate_check_layout(field)

  field <- arc_make_plot_order(field, plot_type = plot_type)
  field$rep <- ifelse(field$type == "check", field$block, 1L)
  field$is_control <- ifelse(field$type == "check", 1L, 0L)

  field
}
