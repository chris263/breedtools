#' Run Lin-style environment stratification from a Breedbase phenotype file
#'
#' Reads a Breedbase phenotype file, prepares the phenotype data,
#' detects the experimental design, calculates ANOVA, and groups environments
#' using a Lin-style genotype-by-environment interaction test.
#'
#' @param file Path to a `.xlsx`, `.xls`, or `.csv` phenotype file, or an
#'   already loaded data frame.
#' @param trait_col Trait column name.
#' @param alpha Significance level. Default is `0.05`.
#' @param normalize_names Logical. If `TRUE`, normalize Breedbase-style column names.
#'
#' @return A list containing standardized phenotype data, environment metadata,
#' ANOVA results, pairwise environment tests, group summary, group membership,
#' ungrouped environments, summary, and message.
#'
#' @examples
#' \dontrun{
#' result <- run_breedbase_environment_stratification(
#'   file = "breedbase_trials_pheno.xlsx",
#'   trait_col = "phenotype",
#'   alpha = 0.05
#' )
#'
#' result$group_summary
#' result$group_membership
#' }
#'
#' @export
run_breedbase_environment_stratification <- function(
    file,
    trait_col,
    alpha = 0.05,
    normalize_names = TRUE
) {
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be a number between 0 and 1.", call. = FALSE)
  }

  pheno_raw <- load_or_use_data_frame(file, "file")

  prepared <- prepare_environment_stratification_data(
    pheno = pheno_raw,
    trait = trait_col,
    normalize_names = normalize_names
  )

  df <- prepared$data
  env_info <- prepared$env_info
  trait_col_detected <- prepared$trait_col

  anova_results <- calculate_environment_stratification_anova(df)

  summary <- data.frame(
    trait = trait_col_detected,
    alpha = alpha,
    n_environments = dplyr::n_distinct(df$environment),
    n_genotypes = dplyr::n_distinct(df$accession_name),
    n_observations = nrow(df),
    stringsAsFactors = FALSE
  )

  if (summary$n_environments < 2) {
    return(list(
      pheno = df,
      env_info = env_info,
      pairwise = empty_pairwise(),
      group_summary = empty_group_summary(),
      group_membership = empty_group_membership(),
      ungrouped = env_info,
      summary = summary,
      anova = anova_results,
      message = "The selected trait must be measured in at least two environments."
    ))
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

  list(
    pheno = df,
    env_info = env_info,
    pairwise = lin_results$pairwise,
    group_summary = lin_results$group_summary,
    group_membership = lin_results$group_membership,
    ungrouped = lin_results$ungrouped,
    summary = summary,
    anova = anova_results,
    message = message_text
  )
}
