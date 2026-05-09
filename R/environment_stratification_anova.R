has_factor_levels <- function(data, column) {
  column %in% names(data) &&
    dplyr::n_distinct(data[[column]][!is.na(data[[column]])]) > 1
}

has_numeric_levels <- function(data, column) {
  column %in% names(data) &&
    dplyr::n_distinct(data[[column]][!is.na(data[[column]])]) > 1
}

#' Detect field trial design for environment stratification
#'
#' Detects whether the data look like CRD, RCBD, row-column, alpha-lattice,
#' or augmented/block-only design.
#'
#' @param data Standardized phenotype data.
#'
#' @return A list with `design`, `formula`, and `message`.
#'
#' @export
detect_design <- function(data) {
  design_text <- ""

  if ("study_design" %in% names(data)) {
    design_text <- paste(
      unique(tolower(as.character(data$study_design))),
      collapse = " "
    )
  }

  has_environment <- has_factor_levels(data, "environment")
  has_accession <- has_factor_levels(data, "accession_name")
  has_rep <- has_factor_levels(data, "rep_number")
  has_block <- has_factor_levels(data, "block_number")
  has_row <- has_numeric_levels(data, "row_number")
  has_col <- has_numeric_levels(data, "col_number")

  block_differs_from_rep <- has_rep &&
    has_block &&
    any(
      as.character(data$block_number) != as.character(data$rep_number),
      na.rm = TRUE
    )

  is_row_column_design <- grepl(
    "row[- ]?column|row.*column|column.*row|spatial",
    design_text
  )

  if (is_row_column_design) {
    design_label <- "Row-column"
    design_message <- paste(
      "Detected row and column layout;",
      "row and column are fitted within environment when those terms have at least two levels."
    )

    design_terms <- c(
      if (has_rep) "environment:rep_number",
      if (has_row) "environment:row_number",
      if (has_col) "environment:col_number"
    )
  } else if (grepl("rcbd|randomized complete block|randomised complete block", design_text)) {
    design_label <- "RCBD"
    design_message <- paste(
      "Detected randomized complete block layout;",
      "blocks are fitted within environment when block has at least two levels."
    )

    design_terms <- c(
      if (has_block) {
        "environment:block_number"
      } else if (has_rep) {
        "environment:rep_number"
      }
    )
  } else if (grepl("alpha|lattice|incomplete", design_text) || block_differs_from_rep) {
    design_label <- "Incomplete block / alpha-lattice"
    design_message <- paste(
      "Detected replicate and block layout;",
      "blocks are fitted within replicate and environment when those terms have at least two levels."
    )

    design_terms <- c(
      if (has_rep) "environment:rep_number",
      if (has_rep && has_block) "environment:rep_number:block_number"
    )
  } else if (grepl("augmented", design_text) || (has_block && !has_rep)) {
    design_label <- "Augmented / block-only"
    design_message <- paste(
      "Detected block-only layout;",
      "blocks are fitted within environment when block has at least two levels."
    )

    design_terms <- c(
      if (has_block) "environment:block_number"
    )
  } else if (has_rep) {
    design_label <- "RCBD"
    design_message <- paste(
      "Detected randomized complete block layout;",
      "blocks are fitted within environment when available."
    )

    design_terms <- c(
      if (has_block) {
        "environment:block_number"
      } else {
        "environment:rep_number"
      }
    )
  } else {
    design_label <- "CRD"
    design_message <- paste(
      "No usable blocking, replicate, row, or column layout detected;",
      "using CRD model."
    )

    design_terms <- character()
  }

  model_terms <- c(
    if (has_environment) "environment",
    design_terms,
    if (has_accession) "accession_name",
    if (has_environment && has_accession) "environment:accession_name"
  )

  model_terms <- unique(model_terms)

  formula_text <- if (length(model_terms) > 0) {
    paste("phenotype ~", paste(model_terms, collapse = " + "))
  } else {
    "phenotype ~ 1"
  }

  list(
    design = design_label,
    formula = stats::as.formula(formula_text),
    message = design_message
  )
}

#' Prepare design factors
#'
#' @param data Standardized phenotype data.
#'
#' @return Data frame with design columns converted to factors.
#'
#' @export
prepare_design_factors <- function(data) {
  data |>
    dplyr::mutate(
      environment = factor(.data$environment),
      accession_name = factor(.data$accession_name),
      rep_number = factor(.data$rep_number),
      block_number = factor(.data$block_number),
      row_number = factor(.data$row_number),
      col_number = factor(.data$col_number)
    )
}

#' Calculate design-aware ANOVA
#'
#' @param data Standardized phenotype data.
#'
#' @return ANOVA table as a data frame.
#'
#' @export
calculate_environment_stratification_anova <- function(data) {
  d <- prepare_design_factors(data)
  design_info <- detect_design(d)

  fit <- tryCatch(
    stats::lm(design_info$formula, data = d, na.action = stats::na.omit),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      design = design_info$design,
      term = "ERROR",
      df = NA_real_,
      sum_sq = NA_real_,
      mean_sq = NA_real_,
      f_value = NA_real_,
      p_value = NA_real_,
      message = fit$message
    ))
  }

  tab <- as.data.frame(stats::anova(fit))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  tab |>
    dplyr::transmute(
      design = design_info$design,
      term = .data$term,
      df = .data$Df,
      sum_sq = .data$`Sum Sq`,
      mean_sq = .data$`Mean Sq`,
      f_value = .data$`F value`,
      p_value = .data$`Pr(>F)`,
      message = design_info$message
    )
}
