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

#' Analyze an augmented row-column design
#'
#' Fits the fixed-effects model:
#' `y = mean + rowgroup + colgroup + rowgroup:colgroup + rowgroup:row +
#' colgroup:col + genotype + error`.
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
#' table, plot-level fitted values and residuals, genotype summaries, and model
#' diagnostics.
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
    "response ~ rowgroup + colgroup + rowgroup:colgroup + rowgroup:row + colgroup:col + genotype"
  )

  fit <- stats::lm(model_formula, data = analysis_data, na.action = stats::na.omit)
  anova_table <- arc_lm_anova_table(fit)

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
    sigma = summary(fit)$sigma,
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    rank_deficient = fit$rank < length(stats::coef(fit))
  )

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
    plot_results = plot_results |> dplyr::select(-dplyr::all_of(".analysis_response")),
    genotype_summary = genotype_summary,
    diagnostics = diagnostics
  )

  class(result) <- c("augmented_row_column_analysis", class(result))
  result
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
  tab <- stats::anova(fit)
  out <- tibble::tibble(
    term = rownames(tab),
    df = as.numeric(tab[["Df"]]),
    sum_sq = as.numeric(tab[["Sum Sq"]]),
    mean_sq = as.numeric(tab[["Mean Sq"]]),
    f_value = as.numeric(tab[["F value"]]),
    p_value = as.numeric(tab[["Pr(>F)"]])
  )
  rownames(out) <- NULL
  out
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
  n_checks <- length(controls)
  check_slots <- data.frame(
    check_idx = integer(0),
    block = integer(0),
    row = integer(0),
    col = integer(0)
  )

  for (b in sort(unique(field$block))) {
    block_idx <- which(field$block == b)
    block_rows <- sort(unique(field$row[block_idx]))
    block_cols <- sort(unique(field$col[block_idx]))

    check_rows <- sample(block_rows, n_checks, replace = FALSE)
    check_cols <- sample(block_cols, n_checks, replace = FALSE)

    local_idx <- match(
      paste(check_rows, check_cols),
      paste(field$row[block_idx], field$col[block_idx])
    )
    check_idx <- block_idx[local_idx]

    field$type[check_idx] <- "check"
    check_slots <- rbind(
      check_slots,
      data.frame(
        check_idx = check_idx,
        block = b,
        row = check_rows,
        col = check_cols
      )
    )
  }

  check_slots <- arc_assign_check_names(check_slots, controls)
  field$trt[check_slots$check_idx] <- check_slots$trt

  entry_idx <- which(field$type == "entry")
  field$trt[entry_idx] <- sample(treatments, length(entry_idx), replace = FALSE)

  arc_validate_check_layout(field)

  field <- arc_make_plot_order(field, plot_type = plot_type)
  field$rep <- ifelse(field$type == "check", field$block, 1L)
  field$is_control <- ifelse(field$type == "check", 1L, 0L)

  field
}
