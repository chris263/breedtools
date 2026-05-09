#' Calculate Lin genotype-by-environment sum of squares
#'
#' @param data Standardized phenotype data.
#' @param envs Character vector of environments to test.
#'
#' @return A list with a summary table and common genotype names.
#'
#' @export
lin_ss_ge <- function(data, envs) {
  d <- data |>
    dplyr::filter(.data$environment %in% envs)

  common_genotypes <- d |>
    dplyr::distinct(.data$accession_name, .data$environment) |>
    dplyr::count(.data$accession_name, name = "n_env") |>
    dplyr::filter(.data$n_env == length(envs)) |>
    dplyr::pull(.data$accession_name)

  d <- d |>
    dplyr::filter(.data$accession_name %in% common_genotypes)

  n_gen <- length(common_genotypes)
  n_env <- length(envs)

  if (n_gen < 2 || n_env < 2) {
    return(list(
      summary = data.frame(
        environments = paste(envs, collapse = ", "),
        n_env = n_env,
        n_genotypes = n_gen,
        r_eff = NA_real_,
        ss_ge = NA_real_,
        df_ge = NA_real_,
        ms_ge = NA_real_,
        message = "Not enough genotypes or environments"
      ),
      common_genotypes = common_genotypes
    ))
  }

  cell_means <- d |>
    dplyr::group_by(.data$accession_name, .data$environment) |>
    dplyr::summarise(
      mean_y = mean(.data$phenotype, na.rm = TRUE),
      n_rep = dplyr::n(),
      .groups = "drop"
    )

  wide <- cell_means |>
    dplyr::select(.data$accession_name, .data$environment, .data$mean_y) |>
    tidyr::pivot_wider(
      names_from = .data$environment,
      values_from = .data$mean_y
    ) |>
    dplyr::arrange(.data$accession_name)

  Y <- wide |>
    dplyr::select(dplyr::all_of(envs)) |>
    as.matrix()

  storage.mode(Y) <- "numeric"

  r_eff <- harmonic_mean(cell_means$n_rep)

  row_mean <- rowMeans(Y, na.rm = TRUE)
  col_mean <- colMeans(Y, na.rm = TRUE)
  grand_mean <- mean(Y, na.rm = TRUE)

  interaction_matrix <- sweep(Y, 1, row_mean, "-")
  interaction_matrix <- sweep(interaction_matrix, 2, col_mean, "-")
  interaction_matrix <- interaction_matrix + grand_mean

  ss_ge <- r_eff * sum(interaction_matrix^2, na.rm = TRUE)
  df_ge <- (n_gen - 1) * (n_env - 1)
  ms_ge <- ss_ge / df_ge

  list(
    summary = data.frame(
      environments = paste(envs, collapse = ", "),
      n_env = n_env,
      n_genotypes = n_gen,
      r_eff = r_eff,
      ss_ge = ss_ge,
      df_ge = df_ge,
      ms_ge = ms_ge,
      message = "OK"
    ),
    common_genotypes = common_genotypes
  )
}

#' Calculate residual error MSE for a Lin environment group
#'
#' @param data Standardized phenotype data.
#' @param envs Character vector of environments.
#' @param common_genotypes Character vector of common genotypes.
#'
#' @return A list with MSE, residual df, and message.
#'
#' @export
lin_error_mse <- function(data, envs, common_genotypes) {
  d <- data |>
    dplyr::filter(
      .data$environment %in% envs,
      .data$accession_name %in% common_genotypes
    ) |>
    prepare_design_factors()

  design <- detect_design(d)

  fit <- tryCatch(
    stats::lm(design$formula, data = d, na.action = stats::na.omit),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      mse = NA_real_,
      df_error = NA_real_,
      message = fit$message
    ))
  }

  df_error <- stats::df.residual(fit)

  if (df_error <= 0) {
    return(list(
      mse = NA_real_,
      df_error = df_error,
      message = "No residual degrees of freedom"
    ))
  }

  list(
    mse = sum(stats::residuals(fit)^2, na.rm = TRUE) / df_error,
    df_error = df_error,
    message = design$message
  )
}

#' Test one Lin environment group
#'
#' @param data Standardized phenotype data.
#' @param envs Character vector of environments.
#' @param alpha Significance level.
#'
#' @return Data frame with Lin GE test result.
#'
#' @export
lin_test_group <- function(data, envs, alpha = 0.05) {
  ss_obj <- lin_ss_ge(data, envs)
  ss_tab <- ss_obj$summary
  common_genotypes <- ss_obj$common_genotypes

  if (is.na(ss_tab$ss_ge)) {
    ss_tab$f_value <- NA_real_
    ss_tab$p_value <- NA_real_
    ss_tab$compatible <- NA
    return(ss_tab)
  }

  mse_obj <- lin_error_mse(
    data = data,
    envs = envs,
    common_genotypes = common_genotypes
  )

  if (is.na(mse_obj$mse)) {
    ss_tab$mse_error <- NA_real_
    ss_tab$df_error <- mse_obj$df_error
    ss_tab$f_value <- NA_real_
    ss_tab$p_value <- NA_real_
    ss_tab$compatible <- NA
    ss_tab$message <- mse_obj$message
    return(ss_tab)
  }

  f_value <- ss_tab$ms_ge / mse_obj$mse

  p_value <- stats::pf(
    q = f_value,
    df1 = ss_tab$df_ge,
    df2 = mse_obj$df_error,
    lower.tail = FALSE
  )

  ss_tab$mse_error <- mse_obj$mse
  ss_tab$df_error <- mse_obj$df_error
  ss_tab$f_value <- f_value
  ss_tab$p_value <- p_value
  ss_tab$compatible <- p_value >= alpha

  ss_tab
}

#' Group environments using Lin-style GE compatibility testing
#'
#' @param data Standardized phenotype data.
#' @param alpha Significance level.
#'
#' @return A list with pairwise tests, group summary, group membership, and
#' ungrouped environments.
#'
#' @export
lin_group_environments <- function(data, alpha = 0.05) {
  all_envs <- sort(unique(as.character(data$environment)))

  if (length(all_envs) < 2) {
    return(list(
      pairwise = empty_pairwise(),
      group_summary = empty_group_summary(),
      group_membership = empty_group_membership(),
      ungrouped = data.frame(environment = all_envs)
    ))
  }

  env_pairs <- utils::combn(all_envs, 2, simplify = FALSE)

  pairwise <- dplyr::bind_rows(lapply(env_pairs, function(x) {
    x <- as.character(x)

    res <- lin_test_group(
      data = data,
      envs = x,
      alpha = alpha
    )

    res$env1 <- x[1]
    res$env2 <- x[2]
    res
  })) |>
    dplyr::select(.data$env1, .data$env2, dplyr::everything()) |>
    dplyr::arrange(.data$ss_ge)

  remaining_envs <- all_envs
  groups <- list()
  group_tests <- list()
  group_id <- 1

  while (length(remaining_envs) >= 2) {
    candidate_pairs <- pairwise |>
      dplyr::filter(
        .data$env1 %in% remaining_envs,
        .data$env2 %in% remaining_envs,
        .data$compatible == TRUE
      ) |>
      dplyr::arrange(.data$ss_ge)

    if (nrow(candidate_pairs) == 0) {
      break
    }

    current_group <- c(candidate_pairs$env1[1], candidate_pairs$env2[1])

    repeat {
      candidates_to_add <- setdiff(remaining_envs, current_group)

      if (length(candidates_to_add) == 0) {
        break
      }

      add_tests <- dplyr::bind_rows(lapply(candidates_to_add, function(candidate_env) {
        res <- lin_test_group(
          data = data,
          envs = c(current_group, candidate_env),
          alpha = alpha
        )

        res$candidate_env <- candidate_env
        res
      })) |>
        dplyr::filter(.data$compatible == TRUE) |>
        dplyr::arrange(.data$ss_ge)

      if (nrow(add_tests) == 0) {
        break
      }

      current_group <- c(current_group, add_tests$candidate_env[1])
    }

    final_test <- lin_test_group(
      data = data,
      envs = current_group,
      alpha = alpha
    )

    final_test$group_id <- paste0("Group_", group_id)

    groups[[group_id]] <- current_group
    group_tests[[group_id]] <- final_test

    remaining_envs <- setdiff(remaining_envs, current_group)
    group_id <- group_id + 1
  }

  if (length(group_tests) > 0) {
    group_summary <- dplyr::bind_rows(group_tests) |>
      dplyr::select(.data$group_id, dplyr::everything())

    group_membership <- dplyr::bind_rows(lapply(seq_along(groups), function(i) {
      data.frame(
        group_id = paste0("Group_", i),
        environment = groups[[i]]
      )
    }))
  } else {
    group_summary <- empty_group_summary()
    group_membership <- empty_group_membership()
  }

  list(
    pairwise = pairwise,
    group_summary = group_summary,
    group_membership = group_membership,
    ungrouped = data.frame(environment = remaining_envs)
  )
}
