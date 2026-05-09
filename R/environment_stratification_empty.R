#' Empty pairwise Lin result table
#'
#' @return Empty data frame.
#'
#' @export
empty_pairwise <- function() {
  data.frame(
    env1 = character(),
    env2 = character(),
    environments = character(),
    n_env = integer(),
    n_genotypes = integer(),
    r_eff = numeric(),
    ss_ge = numeric(),
    df_ge = numeric(),
    ms_ge = numeric(),
    message = character(),
    mse_error = numeric(),
    df_error = numeric(),
    f_value = numeric(),
    p_value = numeric(),
    compatible = logical()
  )
}

#' Empty group summary table
#'
#' @return Empty data frame.
#'
#' @export
empty_group_summary <- function() {
  data.frame(
    group_id = character(),
    environments = character(),
    n_env = integer(),
    n_genotypes = integer(),
    r_eff = numeric(),
    ss_ge = numeric(),
    df_ge = numeric(),
    ms_ge = numeric(),
    message = character(),
    mse_error = numeric(),
    df_error = numeric(),
    f_value = numeric(),
    p_value = numeric(),
    compatible = logical()
  )
}

#' Empty group membership table
#'
#' @return Empty data frame.
#'
#' @export
empty_group_membership <- function() {
  data.frame(
    group_id = character(),
    environment = character()
  )
}

#' Empty ungrouped environment table
#'
#' @return Empty data frame.
#'
#' @export
empty_ungrouped <- function() {
  data.frame(
    environment = character(),
    location = character(),
    trial = character(),
    year = character()
  )
}

#' Empty ANOVA result table
#'
#' @return Empty data frame.
#'
#' @export
empty_anova <- function() {
  data.frame(
    design = character(),
    term = character(),
    df = numeric(),
    sum_sq = numeric(),
    mean_sq = numeric(),
    f_value = numeric(),
    p_value = numeric(),
    message = character()
  )
}
