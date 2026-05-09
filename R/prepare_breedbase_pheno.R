#' Prepare Breedbase phenotype data for multi-environment analysis
#'
#' This function standardizes Breedbase phenotype data into five columns:
#' `environment`, `genotype`, `replicate`, `block`, and `value`.
#'
#' @param pheno_raw A raw Breedbase phenotype data frame, usually read from
#' a Breedbase Excel file.
#' @param trait_col Name of the trait column to analyze.
#' @param location_col Optional location/environment column name. If `NULL`,
#'   the function tries to guess it.
#' @param genotype_col Optional genotype/accession column name. If `NULL`,
#'   the function tries to guess it.
#' @param rep_col Optional replicate column name.
#' @param block_col Optional block column name.
#'
#' @return A tibble with standardized phenotype data.
#'
#' @examples
#' df <- data.frame(
#'   location = c("Loc1", "Loc1", "Loc2", "Loc2"),
#'   germplasmName = c("G1", "G2", "G1", "G2"),
#'   rep = c(1, 1, 1, 1),
#'   phenotype = c(10, 12, 15, 14)
#' )
#'
#' prepare_breedbase_pheno(df, trait_col = "phenotype")
#'
#' @export
prepare_breedbase_pheno <- function(
    pheno_raw,
    trait_col,
    location_col = NULL,
    genotype_col = NULL,
    rep_col = NULL,
    block_col = NULL
) {
  if (!trait_col %in% names(pheno_raw)) {
    stop("Trait column not found: ", trait_col, call. = FALSE)
  }

  if (is.null(location_col)) {
    location_col <- guess_column(
      pheno_raw,
      c(
        "location",
        "location_name",
        "locationName",
        "locationDbId",
        "location_db_id",
        "environment",
        "env",
        "trial_location",
        "studyLocation"
      )
    )
  }

  if (is.null(genotype_col)) {
    genotype_col <- guess_column(
      pheno_raw,
      c(
        "germplasmName",
        "germplasm_name",
        "germplasm",
        "accession",
        "accession_name",
        "stock_name",
        "genotype",
        "all_entries",
        "entry",
        "line",
        "variety"
      )
    )
  }

  if (is.null(rep_col)) {
    rep_col <- guess_column(
      pheno_raw,
      c(
        "rep",
        "replicate",
        "replication",
        "blockNumber",
        "block_number",
        "block",
        "plot_rep"
      )
    )
  }

  if (is.null(block_col)) {
    block_col <- guess_column(
      pheno_raw,
      c(
        "block",
        "blockNumber",
        "block_number",
        "rep",
        "replicate",
        "replication"
      )
    )
  }

  if (is.null(location_col)) {
    stop("Could not detect location/environment column.", call. = FALSE)
  }

  if (is.null(genotype_col)) {
    stop("Could not detect genotype/accession column.", call. = FALSE)
  }

  pheno <- pheno_raw |>
    dplyr::transmute(
      environment = as.factor(.data[[location_col]]),
      genotype = as.factor(.data[[genotype_col]]),
      replicate = if (!is.null(rep_col)) {
        as.factor(.data[[rep_col]])
      } else {
        factor(1)
      },
      block = if (!is.null(block_col)) {
        as.factor(.data[[block_col]])
      } else {
        factor(1)
      },
      value = suppressWarnings(as.numeric(.data[[trait_col]]))
    ) |>
    dplyr::filter(
      !is.na(.data$environment),
      !is.na(.data$genotype),
      !is.na(.data$value)
    )

  attr(pheno, "detected_columns") <- list(
    location_col = location_col,
    genotype_col = genotype_col,
    rep_col = rep_col,
    block_col = block_col,
    trait_col = trait_col
  )

  tibble::as_tibble(pheno)
}
