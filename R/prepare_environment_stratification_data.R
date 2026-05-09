#' Prepare phenotype data for Lin environment stratification
#'
#' Reads a Breedbase-style phenotype data frame and standardizes it for
#' environment stratification.
#'
#' @param pheno A phenotype data frame.
#' @param trait Trait column name.
#' @param normalize_names Logical. If `TRUE`, normalizes column names.
#'
#' @return A list with standardized phenotype data, environment metadata, and
#' selected trait column.
#'
#' @export
prepare_environment_stratification_data <- function(
    pheno,
    trait,
    normalize_names = TRUE
) {
  if (isTRUE(normalize_names)) {
    colnames(pheno) <- normalize_column_name(colnames(pheno))
    trait <- normalize_column_name(gsub("\\.", " ", trait))
  }

  trait_col <- find_column(pheno, c(trait), "selected trait")

  accession_col <- find_column(
    pheno,
    c("germplasmName", "accession_name", "accessionName", "stockName", "all_entries"),
    "accession"
  )

  location_col <- find_column(
    pheno,
    c("locationName", "location", "studyLocation"),
    "location"
  )

  trial_col <- c("studyName", "trialName", "trial_name", "projectName")
  year_col <- c("year", "Year", "season")
  design_col <- c("studyDesign", "study_design", "trialDesign", "design")
  rep_col <- c("replicate", "rep_number", "repNumber", "rep")
  block_col <- c("blockNumber", "block_number", "block")
  row_col <- c("rowNumber", "row_number", "row", "Y", "y")
  col_col <- c("colNumber", "col_number", "col", "X", "x")

  location <- as.character(pheno[[location_col]])

  trial <- rep("", nrow(pheno))
  year <- rep("", nrow(pheno))

  environment_parts <- list(location)

  if (any(trial_col %in% colnames(pheno))) {
    trial <- as.character(pheno[[trial_col[trial_col %in% colnames(pheno)][1]]])
    environment_parts <- c(environment_parts, list(trial))
  }

  if (any(year_col %in% colnames(pheno))) {
    year <- as.character(pheno[[year_col[year_col %in% colnames(pheno)][1]]])
    environment_parts <- c(environment_parts, list(year))
  }

  environment <- do.call(paste, c(environment_parts, sep = "_"))

  study_design <- if (any(design_col %in% colnames(pheno))) {
    pheno[[design_col[design_col %in% colnames(pheno)][1]]]
  } else {
    ""
  }

  rep_number <- if (any(rep_col %in% colnames(pheno))) {
    pheno[[rep_col[rep_col %in% colnames(pheno)][1]]]
  } else {
    "1"
  }

  block_number <- if (any(block_col %in% colnames(pheno))) {
    pheno[[block_col[block_col %in% colnames(pheno)][1]]]
  } else {
    "1"
  }

  row_number <- if (any(row_col %in% colnames(pheno))) {
    pheno[[row_col[row_col %in% colnames(pheno)][1]]]
  } else {
    "1"
  }

  col_number <- if (any(col_col %in% colnames(pheno))) {
    pheno[[col_col[col_col %in% colnames(pheno)][1]]]
  } else {
    "1"
  }

  df <- data.frame(
    environment = as.character(environment),
    location = location,
    trial = trial,
    year = year,
    study_design = as.character(study_design),
    accession_name = as.character(pheno[[accession_col]]),
    rep_number = as.character(rep_number),
    block_number = as.character(block_number),
    row_number = as.character(row_number),
    col_number = as.character(col_number),
    phenotype = as.numeric(gsub(",", ".", as.character(pheno[[trait_col]]))),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(.data$environment),
      !is.na(.data$accession_name),
      !is.na(.data$phenotype),
      .data$accession_name != ""
    )

  env_info <- df |>
    dplyr::distinct(.data$environment, .data$location, .data$trial, .data$year) |>
    dplyr::arrange(.data$location, .data$trial, .data$year) |>
    complete_environment_info()

  list(
    data = df,
    env_info = env_info,
    trait_col = trait_col
  )
}
