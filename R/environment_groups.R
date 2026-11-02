#' Test all compatible environment groups using Hartley Fmax
#'
#' Tests all combinations of environments and identifies groups with homogeneous
#' residual variances.
#'
#' @param env_anova A data frame from [analyze_environment_anova()].
#' @param alpha Significance level.
#' @param max_group_size Maximum environment group size to test.
#' @param n_sim Number of simulations for Hartley critical value.
#' @param seed Random seed.
#'
#' @return A tibble with all tested groups.
#'
#' @examples
#' env_anova <- data.frame(
#'   environment = c("Loc1", "Loc2", "Loc3"),
#'   ms_error = c(10, 12, 30),
#'   df_error = c(20, 20, 20)
#' )
#'
#' make_environment_groups(env_anova, n_sim = 2000)
#'
#' @export
make_environment_groups <- function(
    env_anova,
    alpha = 0.05,
    max_group_size = Inf,
    n_sim = 100000,
    seed = 123
) {
  envs <- env_anova$environment
  n_env <- length(envs)

  if (n_env < 2) {
    stop("At least two environments are required.", call. = FALSE)
  }

  all_results <- list()
  counter <- 1

  for (group_size in 2:n_env) {
    if (group_size > max_group_size) {
      next
    }

    cmb <- utils::combn(envs, group_size, simplify = FALSE)

    for (x in cmb) {
      sub_df <- env_anova |>
        dplyr::filter(.data$environment %in% x)

      test <- hartley_fmax_test(
        sub_df,
        alpha = alpha,
        n_sim = n_sim,
        seed = seed
      )

      all_results[[counter]] <- test |>
        dplyr::mutate(
          group_id = paste0("G", counter),
          environment_list = paste(x, collapse = " | ")
        )

      counter <- counter + 1
    }
  }

  dplyr::bind_rows(all_results) |>
    dplyr::select(
      "group_id",
      "k",
      "environment_list",
      "min_ms_error",
      "max_ms_error",
      "fmax_observed",
      "fmax_critical",
      "alpha",
      "homogeneous"
    ) |>
    dplyr::arrange(dplyr::desc(.data$k), .data$fmax_observed)
}

#' Get compatible environment groups
#'
#' @param group_tests Output from [make_environment_groups()].
#'
#' @return A tibble containing only homogeneous groups.
#'
#' @examples
#' env_anova <- data.frame(
#'   environment = c("Loc1", "Loc2", "Loc3"),
#'   ms_error = c(10, 12, 30),
#'   df_error = c(20, 20, 20)
#' )
#'
#' group_tests <- make_environment_groups(env_anova, n_sim = 2000)
#' get_compatible_groups(group_tests)
#'
#' @export
get_compatible_groups <- function(group_tests) {
  group_tests |>
    dplyr::filter(.data$homogeneous) |>
    dplyr::arrange(dplyr::desc(.data$k), .data$fmax_observed)
}

#' Select maximal compatible environment groups
#'
#' Removes smaller compatible groups that are fully contained inside larger
#' compatible groups.
#'
#' @param compatible_groups Output from [get_compatible_groups()].
#'
#' @return A tibble with maximal compatible groups only.
#'
#' @examples
#' env_anova <- data.frame(
#'   environment = c("Loc1", "Loc2", "Loc3"),
#'   ms_error = c(10, 12, 30),
#'   df_error = c(20, 20, 20)
#' )
#'
#' group_tests <- make_environment_groups(env_anova, n_sim = 2000)
#' get_maximal_groups(get_compatible_groups(group_tests))
#'
#' @export
get_maximal_groups <- function(compatible_groups) {
  if (nrow(compatible_groups) == 0) {
    return(compatible_groups)
  }

  group_sets <- strsplit(compatible_groups$environment_list, " \\| ")

  keep <- rep(TRUE, length(group_sets))

  for (i in seq_along(group_sets)) {
    for (j in seq_along(group_sets)) {
      if (i == j) {
        next
      }

      set_i <- group_sets[[i]]
      set_j <- group_sets[[j]]

      if (length(set_i) < length(set_j) && all(set_i %in% set_j)) {
        keep[i] <- FALSE
      }
    }
  }

  compatible_groups[keep, ]
}
