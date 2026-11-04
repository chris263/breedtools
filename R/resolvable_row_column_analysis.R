#' Analyze a resolvable row-column trial with RC, LV + RC or AR + RC models
#'
#' Fits the variance models of Williams, John and Whitaker (2006) by residual
#' maximum likelihood (REML), with replicates and treatments fixed:
#' \deqn{E(y) = R\pi + X\tau,\quad var(y) = \sigma^2 V.}
#' \describe{
#'   \item{`model = "rc"`}{The standard resolvable row-column mixed model:
#'     \eqn{V = I + \eta_R P_R + \eta_C P_C}, i.e. random rows and columns
#'     within replicates.}
#'   \item{`model = "lv"`}{The two-dimensional linear variance model
#'     (LV + RC, Eq. 2):
#'     \eqn{V = I + \eta_R P_R + \eta_C P_C - (\phi_R F_R + \phi_C F_C +
#'     \phi_{RC} F_{RC})}, with LV matrices for plots in the same row, in the
#'     same column, and in neither. `lv_terms` chooses which \eqn{\phi} terms
#'     are fitted and `equal_phi = TRUE` uses one \eqn{\phi} for rows and
#'     columns (square plots).}
#'   \item{`model = "ar"`}{Rows and columns plus a separable first-order
#'     autoregressive process (AR + RC, Eq. 4):
#'     \eqn{V = I + \eta_R P_R + \eta_C P_C + \psi D(\rho_C) \otimes
#'     D(\rho_R)} within replicates.}
#' }
#' In every model the identity term is the independent error ("nugget").
#' With `compare = TRUE` the nested sequence RC, LV with one \eqn{\phi},
#' LV with \eqn{\phi_R, \phi_C}, LV with \eqn{\phi_{RC}} as well, and AR + RC
#' are all fitted and compared by REML deviance (likelihood-ratio tests for
#' the nested LV terms, and AIC), which is the sequential fitting process the
#' paper recommends for deciding whether spatial terms are needed.
#'
#' For latinized designs the long columns (or long rows) can be fitted as
#' fixed effects with `long_col`.
#'
#' @param data A design data frame (for example `design$design` from
#'   [resolvable_row_column_design()]) or a design object, with a response.
#' @param response Name of the numeric response column.
#' @param model Variance model reported in `adjusted_means`: `"rc"`, `"lv"`,
#'   or `"ar"`.
#' @param lv_terms For `model = "lv"`, the spatial terms to fit: any of
#'   `"phi_r"`, `"phi_c"`, `"phi_rc"`.
#' @param equal_phi Logical. For `model = "lv"`, use a single \eqn{\phi} for
#'   rows and columns.
#' @param long_col Optional name of a column with long-column (or long-row)
#'   labels to fit as fixed effects, as in latinized designs.
#' @param compare Logical. If `TRUE`, also fit and compare the sequence of
#'   models.
#' @param genotype_col,replicate_col,row_col,col_col Column names.
#'
#' @return An object of class `resolvable_row_column_analysis`, a list with
#' \describe{
#'   \item{plot_results}{One row per plot: the data, the fitted value (fixed
#'     effects plus the predicted row, column and spatial effects) and the
#'     residual \eqn{V^{-1}(y - X\hat\beta)}.}
#'   \item{adjusted_means}{Treatment BLUEs (average over replicates) with
#'     standard errors, ranked best first.}
#'   \item{variance_parameters}{The REML estimates: \eqn{\sigma^2}, the
#'     \eqn{\eta}, \eqn{\phi}, \eqn{\psi} and \eqn{\rho} parameters, and the
#'     implied row and column variance components \eqn{\sigma^2\eta_R/s}
#'     and \eqn{\sigma^2\eta_C/k}.}
#'   \item{model_comparison}{With `compare = TRUE`, the deviance table.}
#'   \item{anova}{Intra-block ANOVA from the fixed model
#'     `replicate + rows(replicate) + columns(replicate) + treatment`.}
#'   \item{diagnostics}{REML deviance, number of variance parameters, mean
#'     SED and approximate 5% LSD.}
#' }
#'
#' @references Williams, E. R., John, J. A. and Whitaker, D. (2006).
#'   Construction of resolvable spatial row-column designs. Biometrics 62,
#'   103-108. \doi{10.1111/j.1541-0420.2005.00393.x}
#'
#' @examples
#' design <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5,
#'                                        n_starts = 1, n_iter = 300, seed = 1)
#' trial <- design$design
#' set.seed(2)
#' trial$y <- 10 + 0.3 * trial$row + stats::rnorm(nrow(trial))
#' fit <- analyze_resolvable_row_column(trial, "y", model = "lv")
#' fit$model_comparison
#' head(fit$adjusted_means)
#'
#' @export
analyze_resolvable_row_column <- function(
    data,
    response,
    model = c("rc", "lv", "ar"),
    lv_terms = c("phi_r", "phi_c"),
    equal_phi = FALSE,
    long_col = NULL,
    compare = TRUE,
    genotype_col = "all_entries",
    replicate_col = "replicate",
    row_col = "row_in_rep",
    col_col = "col_in_rep"
) {
  model <- match.arg(model)
  lv_terms <- match.arg(lv_terms, c("phi_r", "phi_c", "phi_rc"), several.ok = TRUE)
  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) data <- data$design
  d <- as.data.frame(data, stringsAsFactors = FALSE)
  required <- c(response, genotype_col, replicate_col, row_col, col_col, long_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d$.y <- suppressWarnings(as.numeric(d[[response]]))
  keep <- !is.na(d$.y) & !is.na(d[[genotype_col]])
  if ("is_filler" %in% names(d)) keep <- keep & d$is_filler != 1
  d <- d[keep, , drop = FALSE]
  if (nrow(d) < 3) stop("Not enough observations.", call. = FALSE)

  d$.rep <- factor(d[[replicate_col]])
  d$.row <- as.integer(d[[row_col]])
  d$.col <- as.integer(d[[col_col]])
  d$.trt <- factor(d[[genotype_col]])
  k <- max(d$.row)
  s <- max(d$.col)

  fixed_terms <- c(if (nlevels(d$.rep) > 1) ".rep", if (!is.null(long_col)) ".long", ".trt")
  if (!is.null(long_col)) d$.long <- factor(d[[long_col]])
  fit0 <- arc_lm_drop_aliased_terms(
    stats::as.formula(paste(".y ~", paste(fixed_terms, collapse = " + "))),
    d
  )
  x <- fit0$x
  x_full <- stats::model.matrix(fit0$terms, fit0$model)
  kept <- match(colnames(x), colnames(x_full))

  mats <- rrc_variance_matrices(d$.rep, d$.row, d$.col, k, s)
  engine <- rrc_reml_engine(d$.y, x, mats, d$.rep, d$.row, d$.col)

  specs <- list(
    rc = list(label = "RC", terms = c("eta_r", "eta_c"), ar = FALSE),
    lv1 = list(label = "RC + LV (single phi)", terms = c("eta_r", "eta_c", "phi"), ar = FALSE),
    lv2 = list(label = "RC + LV (phi_r, phi_c)", terms = c("eta_r", "eta_c", "phi_r", "phi_c"), ar = FALSE),
    lv3 = list(label = "RC + LV (phi_r, phi_c, phi_rc)", terms = c("eta_r", "eta_c", "phi_r", "phi_c", "phi_rc"), ar = FALSE),
    ar = list(label = "RC + AR1 x AR1", terms = c("eta_r", "eta_c"), ar = TRUE)
  )
  chosen <- switch(model,
    rc = specs$rc,
    ar = specs$ar,
    lv = list(
      label = "RC + LV",
      terms = c("eta_r", "eta_c", if (equal_phi) {
        c("phi", if ("phi_rc" %in% lv_terms) "phi_rc")
      } else {
        lv_terms
      }),
      ar = FALSE
    )
  )
  fit <- engine$fit(chosen)

  comparison <- NULL
  if (isTRUE(compare)) {
    fits <- lapply(specs, engine$fit)
    dev <- vapply(fits, `[[`, numeric(1), "deviance")
    npar <- vapply(fits, `[[`, integer(1), "n_par")
    comparison <- tibble::tibble(
      model = unname(vapply(specs, `[[`, character(1), "label")),
      n_variance_parameters = unname(npar),
      deviance = unname(dev),
      AIC = unname(dev + 2 * npar),
      change_vs_rc = unname(dev[["rc"]] - dev)
    )
    nested <- c(lv1 = "rc", lv2 = "rc", lv3 = "lv2")
    comparison$compared_with <- c(NA, "RC", "RC", "RC + LV (phi_r, phi_c)", NA)
    comparison$lrt_p_value <- NA_real_
    for (nm in names(nested)) {
      i <- match(nm, names(specs))
      j <- match(nested[[nm]], names(specs))
      df <- npar[i] - npar[j]
      stat <- max(dev[j] - dev[i], 0)
      comparison$lrt_p_value[i] <- stats::pchisq(stat, df = df, lower.tail = FALSE)
    }
  }

  # Treatment BLUEs: intercept + average replicate (and long-column) effect +
  # treatment effect.
  assign <- attr(x_full, "assign")
  trt_term <- match(".trt", attr(fit0$terms, "term.labels"))
  trt_cols <- which(assign == trt_term)
  levels_t <- levels(d$.trt)
  L <- matrix(colMeans(x_full), length(levels_t), ncol(x_full), byrow = TRUE)
  L[, trt_cols] <- 0
  hit <- match(colnames(x_full)[trt_cols], paste0(".trt", levels_t))
  L[cbind(hit, trt_cols)] <- 1
  qr_full <- qr(x_full)
  estimable <- arc_linear_estimates(qr_full, qr.coef(qr_full, fit0$y), L)$estimable
  Lk <- L[, kept, drop = FALSE]
  estimate <- drop(Lk %*% fit$beta)
  cov_means <- Lk %*% fit$cov_beta %*% t(Lk)
  std_error <- sqrt(pmax(diag(cov_means), 0))
  estimate[!estimable] <- NA_real_
  std_error[!estimable] <- NA_real_
  ok <- which(estimable)
  mean_sed <- prep_mean_sed(cov_means[ok, ok, drop = FALSE])

  adjusted <- tibble::tibble(
    genotype = levels_t,
    n_observations = as.integer(table(d$.trt)[levels_t]),
    adjusted_mean = estimate,
    std_error = std_error
  )
  adjusted <- adjusted[order(-adjusted$adjusted_mean, na.last = TRUE), , drop = FALSE]
  adjusted$rank <- ifelse(is.na(adjusted$adjusted_mean), NA_integer_,
                          rank(-adjusted$adjusted_mean, ties.method = "min", na.last = "keep"))
  adjusted <- adjusted[, c("rank", "genotype", "n_observations", "adjusted_mean", "std_error")]

  # Intra-block ANOVA (rows and columns fixed within replicates).
  d$.rowrep <- interaction(d$.rep, d$.row, drop = TRUE)
  d$.colrep <- interaction(d$.rep, d$.col, drop = TRUE)
  anova_fit <- arc_lm_drop_aliased_terms(
    stats::as.formula(paste(".y ~", paste(c(if (nlevels(d$.rep) > 1) ".rep", ".rowrep", ".colrep", ".trt"), collapse = " + "))),
    d
  )
  anova_table <- arc_lm_anova_table(anova_fit)$table
  labels <- c(.rep = "Replicates", .rowrep = "Rows within replicates",
              .colrep = "Columns within replicates", .trt = "Treatments", Residuals = "Error")
  anova_table$term <- ifelse(anova_table$term %in% names(labels), labels[anova_table$term], anova_table$term)
  anova_df_error <- stats::df.residual(anova_fit)

  par <- fit$parameters
  sigma2 <- fit$sigma2
  variance_parameters <- tibble::tibble(
    parameter = c("sigma2", names(par), "row_variance", "column_variance"),
    estimate = c(sigma2, unname(par),
                 if ("eta_r" %in% names(par)) sigma2 * par[["eta_r"]] / s else NA_real_,
                 if ("eta_c" %in% names(par)) sigma2 * par[["eta_c"]] / k else NA_real_)
  )

  # One row per plot: the observed value, the fitted value (fixed effects
  # plus predicted row, column and spatial effects) and the residual.
  plot_results <- tibble::as_tibble(d[, setdiff(names(d), c(".y", ".rep", ".row", ".col", ".trt", ".long",
                                                        ".rowrep", ".colrep")), drop = FALSE])
  plot_results$fitted <- d$.y - fit$residual
  plot_results$residual <- fit$residual

  result <- list(
    input = list(response = response, model = model, lv_terms = lv_terms, equal_phi = equal_phi,
                 long_col = long_col, genotype_col = genotype_col, replicate_col = replicate_col,
                 row_col = row_col, col_col = col_col),
    model_label = chosen$label,
    plot_results = plot_results,
    adjusted_means = adjusted,
    adjusted_means_vcov = {
      m <- cov_means[ok, ok, drop = FALSE]
      dimnames(m) <- list(levels_t[ok], levels_t[ok])
      m
    },
    variance_parameters = variance_parameters,
    model_comparison = comparison,
    anova = anova_table,
    diagnostics = tibble::tibble(
      model = chosen$label,
      n_observations = nrow(d),
      n_treatments = nlevels(d$.trt),
      n_variance_parameters = fit$n_par,
      reml_deviance = fit$deviance,
      converged = fit$converged,
      intrablock_error_df = anova_df_error,
      mean_sed = mean_sed,
      lsd_5pct = if (!is.na(mean_sed) && anova_df_error > 0) stats::qt(0.975, anova_df_error) * mean_sed else NA_real_
    )
  )
  class(result) <- c("resolvable_row_column_analysis", class(result))
  result
}

#' Print a resolvable row-column analysis
#'
#' @param x An object returned by [analyze_resolvable_row_column()].
#' @param n Number of top-ranked treatments to show.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.resolvable_row_column_analysis <- function(x, n = 10, ...) {
  cat("Resolvable row-column analysis:", x$model_label, "\n\n")
  cat("Variance parameters (REML):\n")
  print(as.data.frame(x$variance_parameters), row.names = FALSE, digits = 4)
  if (!is.null(x$model_comparison)) {
    cat("\nModel comparison (REML deviance):\n")
    print(as.data.frame(x$model_comparison), row.names = FALSE, digits = 4)
  }
  cat(sprintf("\nTop adjusted means (%d treatments):\n", nrow(x$adjusted_means)))
  print(as.data.frame(utils::head(x$adjusted_means, n)), row.names = FALSE, digits = 4)
  cat(sprintf("Mean SED = %.4g; LSD (5%%) = %.4g\n", x$diagnostics$mean_sed, x$diagnostics$lsd_5pct))
  invisible(x)
}

# ---- Internals --------------------------------------------------------------

# Plot-level matrices of the LV + RC model for observed plots (any order):
# P_R and P_C (1/s and 1/k for plots in the same row or column of a
# replicate), and the LV matrices F_R, F_C, F_RC built from
# U_k = 3 L_k / (k^2 - 1), L_k(i, j) = |i - j|.
rrc_variance_matrices <- function(rep, row, col, k, s) {
  same_rep <- outer(as.integer(rep), as.integer(rep), "==")
  same_row <- same_rep & outer(row, row, "==")
  same_col <- same_rep & outer(col, col, "==")
  dr <- abs(outer(row, row, "-"))
  dc <- abs(outer(col, col, "-"))
  Uc <- 3 * dc / (s^2 - 1)
  Ur <- 3 * dr / (k^2 - 1)
  list(
    eta_r = same_row / s,
    eta_c = same_col / k,
    phi_r = -(same_row * Uc),
    phi_c = -(same_col * Ur),
    phi = -(same_row * Uc + same_col * Ur),
    phi_rc = -(same_rep * Uc * Ur),
    same_rep = same_rep,
    dr = dr,
    dc = dc
  )
}

# REML on error contrasts K'y, K'K = I, K'X = 0, so the variance matrix only
# has to be positive definite on contrasts (Eq. 1 of the paper). sigma^2 is
# profiled out. Returns fitted parameters and GLS estimates of the fixed
# effects; replicate-mean directions (in the column space of X) are added to
# V for GLS, which does not change the estimates.
rrc_reml_engine <- function(y, x, mats, rep, row, col) {
  n <- length(y)
  qx <- qr(x)
  p <- qx$rank
  Q <- qr.Q(qx, complete = TRUE)
  K <- Q[, (p + 1L):n, drop = FALSE]
  m <- ncol(K)
  Ky <- drop(crossprod(K, y))
  KMK <- lapply(mats[c("eta_r", "eta_c", "phi_r", "phi_c", "phi", "phi_rc")], function(M) crossprod(K, M %*% K))
  rep_ind <- stats::model.matrix(~ 0 + factor(rep))

  ar_matrix <- function(rho_r, rho_c) mats$same_rep * rho_r^mats$dr * rho_c^mats$dc

  build_V <- function(par, spec, contrast = TRUE) {
    lin <- intersect(names(par), names(KMK))
    if (contrast) {
      V <- diag(m)
      for (nm in lin) V <- V + par[[nm]] * KMK[[nm]]
      if (spec$ar) V <- V + par[["psi"]] * crossprod(K, ar_matrix(par[["rho_r"]], par[["rho_c"]]) %*% K)
    } else {
      V <- diag(n)
      for (nm in lin) V <- V + par[[nm]] * mats[[nm]]
      if (spec$ar) V <- V + par[["psi"]] * ar_matrix(par[["rho_r"]], par[["rho_c"]])
    }
    V
  }

  neg2ll <- function(par, spec) {
    V <- build_V(par, spec)
    ch <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(ch)) return(1e12)
    z <- backsolve(ch, Ky, transpose = TRUE)
    quad <- sum(z^2)
    2 * sum(log(diag(ch))) + m * log(quad / m) + m * (1 + log(2 * pi))
  }

  fit <- function(spec) {
    names_par <- c(spec$terms, if (spec$ar) c("psi", "rho_r", "rho_c"))
    lower <- ifelse(names_par %in% c("rho_r", "rho_c"), 0.01, 0)
    upper <- ifelse(names_par %in% c("rho_r", "rho_c"), 0.99, 1e4)
    starts <- list(
      ifelse(names_par %in% c("rho_r", "rho_c"), 0.5, ifelse(grepl("^phi", names_par), 0.5, 1)),
      ifelse(names_par %in% c("rho_r", "rho_c"), 0.8, ifelse(grepl("^phi", names_par), 0.1, 5))
    )
    best <- NULL
    for (st in starts) {
      names(st) <- names_par
      obj <- function(theta) neg2ll(stats::setNames(theta, names_par), spec)
      if (obj(st) >= 1e12) st[grepl("^phi", names_par)] <- 0
      o <- tryCatch(
        stats::optim(st, obj, method = "L-BFGS-B", lower = lower, upper = upper),
        error = function(e) NULL
      )
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    par <- stats::setNames(best$par, names_par)
    # sigma^2 and GLS fixed effects.
    Vk <- build_V(par, spec)
    ch <- chol(Vk)
    z <- backsolve(ch, Ky, transpose = TRUE)
    sigma2 <- sum(z^2) / m
    Vfull <- build_V(par, spec, contrast = FALSE) + n * tcrossprod(rep_ind)
    chf <- tryCatch(chol(Vfull), error = function(e) NULL)
    if (is.null(chf)) {
      e <- eigen((Vfull + t(Vfull)) / 2, symmetric = TRUE)
      Vi <- e$vectors %*% (t(e$vectors) / pmax(e$values, 1e-8))
    } else {
      Vi <- chol2inv(chf)
    }
    xtvx <- crossprod(x, Vi %*% x)
    cov_beta <- sigma2 * solve(xtvx)
    beta <- drop(solve(xtvx, crossprod(x, Vi %*% y)))
    # Plot residuals (conditional on the predicted row, column and spatial
    # effects): e = sigma2 * (sigma2 V)^-1 (y - X b) = V^-1 (y - X b). The
    # rep-mean directions added to V lie in the column space of X, so they
    # do not change this.
    residual <- drop(Vi %*% (y - x %*% beta))
    list(
      parameters = par,
      sigma2 = sigma2,
      beta = beta,
      residual = residual,
      cov_beta = cov_beta,
      deviance = best$value,
      n_par = length(names_par) + 1L,
      converged = best$convergence == 0
    )
  }
  list(fit = fit, neg2ll = neg2ll)
}
