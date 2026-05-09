#' Analyze one environment using ANOVA
#'
#' Runs an individual ANOVA for one environment. If more than one replicate is
#' present, the model is `value ~ replicate + genotype`; otherwise, the model is
#' `value ~ genotype`.
#'
#' @param df A data frame containing one environment with columns
#'   `environment`, `genotype`, `replicate`, and `value`.
#'
#' @return A tibble with ANOVA summary statistics.
#'
#' @export
analyze_one_environment <- function(df) {
  env_name <- unique(df$environment)

  if (length(env_name) != 1) {
    stop("Input data must contain exactly one environment.", call. = FALSE)
  }

  n_genotypes <- dplyr::n_distinct(df$genotype)

  if (n_genotypes < 2) {
    stop("Environment ", env_name, " has fewer than 2 genotypes.", call. = FALSE)
  }

  if (dplyr::n_distinct(df$replicate) > 1) {
    fit <- stats::aov(value ~ replicate + genotype, data = df)
  } else {
    fit <- stats::aov(value ~ genotype, data = df)
  }

  anova_tab <- stats::anova(fit)

  if (!"genotype" %in% rownames(anova_tab)) {
    stop("No genotype term found for environment: ", env_name, call. = FALSE)
  }

  if (!"Residuals" %in% rownames(anova_tab)) {
    stop("No residual term found for environment: ", env_name, call. = FALSE)
  }

  ms_genotype <- anova_tab["genotype", "Mean Sq"]
  df_genotype <- anova_tab["genotype", "Df"]

  ms_error <- anova_tab["Residuals", "Mean Sq"]
  df_error <- anova_tab["Residuals", "Df"]

  f_value <- anova_tab["genotype", "F value"]
  p_value <- anova_tab["genotype", "Pr(>F)"]

  trait_mean <- mean(df$value, na.rm = TRUE)

  cv <- if (is.na(trait_mean) || trait_mean == 0) {
    NA_real_
  } else {
    sqrt(ms_error) / trait_mean * 100
  }

  rep_table <- table(df$genotype)
  r_eff <- mean(rep_table)

  sigma2_g <- (ms_genotype - ms_error) / r_eff

  if (is.na(sigma2_g) || sigma2_g < 0) {
    sigma2_g <- 0
  }

  tibble::tibble(
    environment = as.character(env_name),
    n_genotypes = n_genotypes,
    r_eff = as.numeric(r_eff),
    df_genotype = df_genotype,
    ms_genotype = ms_genotype,
    df_error = df_error,
    ms_error = ms_error,
    mean = trait_mean,
    cv_percent = cv,
    sigma2_g = sigma2_g,
    f_value = f_value,
    p_value = p_value
  )
}

#' Run individual ANOVA by environment
#'
#' @param pheno Standardized phenotype data from [prepare_breedbase_pheno()].
#'
#' @return A tibble with one ANOVA summary row per environment.
#'
#' @examples
#' df <- data.frame(
#'   environment = factor(rep(c("Loc1", "Loc2"), each = 8)),
#'   genotype = factor(rep(rep(c("G1", "G2"), each = 4), times = 2)),
#'   replicate = factor(rep(rep(c(1, 2), each = 2), times = 4)),
#'   block = factor(1),
#'   value = c(10, 11, 12, 13, 14, 15, 16, 17,
#'             15, 16, 17, 18, 19, 20, 21, 22)
#' )
#'
#' analyze_environment_anova(df)
#'
#' @export
analyze_environment_anova <- function(pheno) {
  required_cols <- c("environment", "genotype", "replicate", "value")

  missing_cols <- setdiff(required_cols, names(pheno))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  pheno |>
    dplyr::group_by(.data$environment) |>
    dplyr::group_split() |>
    purrr::map_dfr(analyze_one_environment)
}
