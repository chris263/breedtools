#' Simulate Hartley Fmax critical value
#'
#' Estimates the Hartley Fmax critical value using simulation. This version
#' can handle different residual degrees of freedom across environments.
#'
#' @param df_error_vec Numeric vector of residual degrees of freedom.
#' @param alpha Significance level. Default is `0.05`.
#' @param n_sim Number of simulations. Default is `100000`.
#' @param seed Random seed.
#'
#' @return A numeric critical value.
#'
#' @export
hartley_critical_sim <- function(
    df_error_vec,
    alpha = 0.05,
    n_sim = 100000,
    seed = 123
) {
  if (length(df_error_vec) < 2) {
    stop("At least two residual degrees of freedom values are required.", call. = FALSE)
  }

  if (any(df_error_vec <= 0, na.rm = TRUE)) {
    stop("All residual degrees of freedom must be greater than zero.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  k <- length(df_error_vec)

  sim_fmax <- replicate(
    n_sim,
    {
      sim_vars <- stats::rchisq(k, df = df_error_vec) / df_error_vec
      max(sim_vars) / min(sim_vars)
    }
  )

  as.numeric(stats::quantile(sim_fmax, probs = 1 - alpha, na.rm = TRUE))
}

#' Hartley Fmax test for residual variance homogeneity
#'
#' @param env_anova A data frame from [analyze_environment_anova()].
#' @param alpha Significance level. Default is `0.05`.
#' @param n_sim Number of simulations for the critical value.
#' @param seed Random seed.
#'
#' @return A tibble with observed Fmax, critical Fmax, and homogeneity decision.
#'
#' @export
hartley_fmax_test <- function(
    env_anova,
    alpha = 0.05,
    n_sim = 100000,
    seed = 123
) {
  required_cols <- c("environment", "ms_error", "df_error")

  missing_cols <- setdiff(required_cols, names(env_anova))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(env_anova) < 2) {
    stop("At least two environments are required.", call. = FALSE)
  }

  if (any(env_anova$ms_error <= 0, na.rm = TRUE)) {
    stop("All residual mean squares must be greater than zero.", call. = FALSE)
  }

  fmax_obs <- max(env_anova$ms_error, na.rm = TRUE) /
    min(env_anova$ms_error, na.rm = TRUE)

  fmax_crit <- hartley_critical_sim(
    df_error_vec = env_anova$df_error,
    alpha = alpha,
    n_sim = n_sim,
    seed = seed
  )

  tibble::tibble(
    environments = paste(env_anova$environment, collapse = ", "),
    k = nrow(env_anova),
    min_ms_error = min(env_anova$ms_error, na.rm = TRUE),
    max_ms_error = max(env_anova$ms_error, na.rm = TRUE),
    fmax_observed = fmax_obs,
    fmax_critical = fmax_crit,
    alpha = alpha,
    homogeneous = fmax_obs <= fmax_crit
  )
}
