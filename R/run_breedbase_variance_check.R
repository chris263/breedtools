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
#' @param seed Random seed.
#' @param export Logical. If `TRUE`, writes CSV files.
#' @param output_dir Directory where CSV files should be written.
#'
#' @return A list containing cleaned data, ANOVA results, Hartley test,
#'   all group tests, compatible groups, maximal compatible groups, and detected
#'   columns.
#'
#' @examples
#' \dontrun{
#' result <- run_breedbase_variance_check(
#'   file = "breedbase_trials_pheno.xlsx",
#'   trait_col = "phenotype"
#' )
#'
#' result$env_anova
#' result$maximal_compatible_groups
#' }
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
    output_dir = "."
) {
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
