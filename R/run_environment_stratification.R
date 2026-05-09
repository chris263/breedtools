#' Run Lin-style environment stratification
#'
#' This function reads phenotype data, detects design structure, calculates
#' ANOVA, tests genotype-by-environment interaction for environment groups,
#' and returns compatible environment groups.
#'
#' @param phenotype_file Path to a `.xlsx`, `.xls`, or `.csv` phenotype file,
#'   or an already loaded data frame.
#' @param trait Trait column name.
#' @param alpha Significance level. Default is `0.05`.
#' @param sep Deprecated. File delimiters are detected from the extension by
#'   `read_breedbase_file()`.
#' @param normalize_names Logical. If `TRUE`, normalize Breedbase column names.
#' @param export_json Logical. If `TRUE`, writes JSON result files.
#' @param output_dir Output directory for JSON files.
#' @param prefix Prefix for JSON output files.
#'
#' @return A list with summary, ANOVA results, pairwise tests, group summary,
#' group membership, ungrouped environments, and message.
#'
#' @examples
#' \dontrun{
#' result <- run_environment_stratification(
#'   phenotype_file = "phenotype.tsv",
#'   trait = "plant height",
#'   alpha = 0.05
#' )
#'
#' result$group_summary
#' result$group_membership
#' }
#'
#' @export
run_environment_stratification <- function(
    phenotype_file,
    trait,
    alpha = 0.05,
    sep = NULL,
    normalize_names = TRUE,
    export_json = FALSE,
    output_dir = ".",
    prefix = "environment_stratification"
) {
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be a number between 0 and 1.", call. = FALSE)
  }

  pheno <- load_or_use_data_frame(phenotype_file, "phenotype_file")

  prepared <- prepare_environment_stratification_data(
    pheno = pheno,
    trait = trait,
    normalize_names = normalize_names
  )

  df <- prepared$data
  env_info <- prepared$env_info
  trait_col <- prepared$trait_col

  anova_results <- calculate_environment_stratification_anova(df)

  summary <- data.frame(
    trait = trait_col,
    alpha = alpha,
    n_environments = dplyr::n_distinct(df$environment),
    n_genotypes = dplyr::n_distinct(df$accession_name),
    n_observations = nrow(df),
    stringsAsFactors = FALSE
  )

  if (summary$n_environments < 2) {
    result <- list(
      pairwise = empty_pairwise(),
      group_summary = empty_group_summary(),
      group_membership = empty_group_membership(),
      ungrouped = env_info,
      summary = summary,
      anova = anova_results,
      message = "The selected trait must be measured in at least two environments."
    )

    if (isTRUE(export_json)) {
      write_environment_stratification_json(result, output_dir, prefix)
    }

    return(result)
  }

  lin_results <- lin_group_environments(
    data = df,
    alpha = alpha
  )

  lin_results$pairwise <- add_pairwise_environment_metadata(
    lin_results$pairwise,
    env_info
  )

  lin_results$group_membership <- add_environment_metadata(
    lin_results$group_membership,
    env_info
  )

  lin_results$ungrouped <- add_environment_metadata(
    lin_results$ungrouped,
    env_info
  )

  lin_results$group_summary <- environment_summary_by_group(
    lin_results$group_summary,
    lin_results$group_membership
  )

  group_count <- nrow(lin_results$group_summary)
  ungrouped_count <- nrow(lin_results$ungrouped)

  message_text <- paste0(
    "Environment stratification finished. ",
    group_count,
    " compatible group(s) found; ",
    ungrouped_count,
    " environment(s) left ungrouped."
  )

  result <- list(
    pairwise = lin_results$pairwise,
    group_summary = lin_results$group_summary,
    group_membership = lin_results$group_membership,
    ungrouped = lin_results$ungrouped,
    summary = summary,
    anova = anova_results,
    message = message_text
  )

  if (isTRUE(export_json)) {
    write_environment_stratification_json(result, output_dir, prefix)
  }

  result
}
