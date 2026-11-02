#' Run Breedbase multi-location phenotype variance check
#'
#' Reads a Breedbase phenotype file, runs individual ANOVA by environment,
#' performs Hartley Fmax tests, and identifies compatible environment groups for
#' joint analysis.
#'
#' @param file Path to a `.xlsx`, `.xls`, or `.csv` phenotype file, or an
#'   already loaded data frame.
#' @param trait_col Trait column to analyze.
#' @param alpha Significance level for Hartley Fmax test.
#' @param location_col Optional location/environment column.
#' @param genotype_col Optional genotype/accession column.
#' @param rep_col Optional replicate column.
#' @param block_col Optional block column.
#' @param max_group_size Maximum environment group size to test.
#' @param n_sim Number of simulations for Hartley critical value.
#' @param seed Random seed for the Hartley critical value simulation. The
#'   user's global random number generator state is not modified.
#' @param export Logical. If `TRUE`, writes CSV files to `output_dir`.
#' @param output_dir Directory where CSV files should be written. Defaults to
#'   `NULL`, in which case nothing is written; it must be supplied when
#'   `export = TRUE`.
#'
#' @return A list containing cleaned data, ANOVA results, Hartley test,
#'   all group tests, compatible groups, maximal compatible groups, and detected
#'   columns.
#'
#' @examples
#' pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
#'
#' result <- run_breedbase_variance_check(
#'   file = pheno_file,
#'   trait_col = "fresh root yield|CO_334:0000013",
#'   n_sim = 2000
#' )
#'
#' result$env_anova
#' result$maximal_compatible_groups
#'
#' # Export CSV files to a temporary directory
#' out_dir <- file.path(tempdir(), "variance_check")
#' result <- run_breedbase_variance_check(
#'   file = pheno_file,
#'   trait_col = "fresh root yield|CO_334:0000013",
#'   n_sim = 2000,
#'   export = TRUE,
#'   output_dir = out_dir
#' )
#' list.files(out_dir)
#' unlink(out_dir, recursive = TRUE)
#'
#' @export
run_breedbase_variance_check <- function(
    file,
    trait_col,
    alpha = 0.05,
    location_col = NULL,
    genotype_col = NULL,
    rep_col = NULL,
    block_col = NULL,
    max_group_size = Inf,
    n_sim = 100000,
    seed = 123,
    export = FALSE,
    output_dir = NULL
) {
  if (isTRUE(export) && is.null(output_dir)) {
    stop("`output_dir` must be supplied when `export = TRUE`.", call. = FALSE)
  }

  pheno_raw <- load_or_use_data_frame(file, "file")

  pheno <- prepare_breedbase_pheno(
    pheno_raw = pheno_raw,
    trait_col = trait_col,
    location_col = location_col,
    genotype_col = genotype_col,
    rep_col = rep_col,
    block_col = block_col
  )

  env_anova <- analyze_environment_anova(pheno)

  overall_hartley <- hartley_fmax_test(
    env_anova,
    alpha = alpha,
    n_sim = n_sim,
    seed = seed
  )

  group_tests <- make_environment_groups(
    env_anova = env_anova,
    alpha = alpha,
    max_group_size = max_group_size,
    n_sim = n_sim,
    seed = seed
  )

  compatible_groups <- get_compatible_groups(group_tests)

  maximal_compatible_groups <- get_maximal_groups(compatible_groups)

  detected_columns <- attr(pheno, "detected_columns")

  result <- list(
    pheno = pheno,
    env_anova = env_anova,
    overall_hartley = overall_hartley,
    group_tests = group_tests,
    compatible_groups = compatible_groups,
    maximal_compatible_groups = maximal_compatible_groups,
    detected_columns = detected_columns
  )

  if (isTRUE(export)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    utils::write.csv(
      env_anova,
      file.path(output_dir, paste0("individual_anova_", trait_col, ".csv")),
      row.names = FALSE
    )

    utils::write.csv(
      group_tests,
      file.path(output_dir, paste0("hartley_all_group_tests_", trait_col, ".csv")),
      row.names = FALSE
    )

    utils::write.csv(
      maximal_compatible_groups,
      file.path(output_dir, paste0("hartley_maximal_compatible_groups_", trait_col, ".csv")),
      row.names = FALSE
    )
  }

  result
}
