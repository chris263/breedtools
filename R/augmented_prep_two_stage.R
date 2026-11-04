#' Two-stage analysis of a multi-location augmented p-rep trial
#'
#' Analyses each location separately (stage 1) and combines the per-location
#' adjusted means across locations (stage 2), weighting them by their
#' precision. This is the approach suggested by Williams, Piepho and
#' Whitaker (2011, Section 3) for adding spatial error components and for
#' heterogeneous error variances across locations: every location gets its
#' own error variance (and, optionally, its own spatial trend), and the
#' stage-2 weights carry that information into the combined means (Möhring
#' and Piepho, 2009).
#'
#' @section Stage 1:
#' \describe{
#'   \item{`stage1 = "blocks"`}{The incomplete-block model of
#'     [analyze_augmented_prep_design()] at each location (replicates and
#'     blocks fixed, genotypes fixed).}
#'   \item{`stage1 = "spatial"`}{A spatial model fitted with the 'SpATS'
#'     package (Rodríguez-Álvarez et al., 2018): genotypes fixed, a smooth
#'     two-dimensional P-spline surface over the field rows and columns
#'     (`PSANOVA`), replicates fixed and incomplete blocks random. It needs
#'     the `row` and `col` field coordinates (use `field_cols` in
#'     [augmented_prep_design()]).}
#' }
#' The stage-1 genotype means are BLUEs with their full covariance matrix.
#' Each location needs some error degrees of freedom, i.e. enough p-rep
#' entries for its blocks. In spatial fits part of the plot-to-plot noise is
#' taken up by the smooth surface, so the reported residual variance is
#' somewhat smaller than in the block model; the stage-1 covariance matrix
#' accounts for this, and the Hartley test on spatial residual variances is
#' approximate.
#'
#' @section Stage 2:
#' The stage-1 means \eqn{\bar y_{gl}} are modelled as
#' \eqn{\bar y_{gl} = \mu + \lambda_l + \gamma_g + (\gamma\lambda)_{gl} + e_{gl}},
#' with fixed locations and genotypes, random genotype-by-location effects
#' with variance \eqn{\sigma^2_{GL}} (estimated by REML) and stage-1 errors
#' \eqn{e} whose covariance is taken from stage 1:
#' \describe{
#'   \item{`weighting = "full"`}{the full stage-1 covariance matrix of each
#'     location (most accurate);}
#'   \item{`weighting = "diagonal"`}{only the stage-1 variances (weights
#'     \eqn{1 / se^2});}
#'   \item{`weighting = "none"`}{an unweighted analysis of the stage-1 means
#'     with a common error, which ignores heterogeneity (for comparison).}
#' }
#' Genotype adjusted means average the location effects with equal weights.
#'
#' @param data A design data frame (usually `design$design` from
#'   [augmented_prep_design()]), or a design object, with a response column.
#' @param response Name of the numeric response column.
#' @param stage1 Stage-1 model: `"blocks"` or `"spatial"`.
#' @param weighting Stage-2 weighting: `"full"`, `"diagonal"`, or `"none"`.
#' @param genotype_by_location `"random"` (default) to estimate a
#'   genotype-by-location variance in stage 2, or `"none"` to omit it.
#' @param genotype_col,location_col,block_col,replicate_col Column names, as
#'   in [analyze_augmented_prep_design()].
#' @param row_col,col_col Field row and column coordinate columns, used by
#'   `stage1 = "spatial"`.
#' @param spatial_nseg Optional number of P-spline segments for columns and
#'   rows, `c(col, row)`. Defaults to half the number of distinct columns and
#'   rows of each location (at least 3, at most 40), which is much faster
#'   than one segment per row or column and usually as accurate.
#' @param spatial_random_blocks Logical. If `TRUE` (default), incomplete
#'   blocks are included as a random effect in the spatial model.
#'
#' @return An object of class `augmented_prep_two_stage`, a list with
#' \describe{
#'   \item{adjusted_means}{Combined genotype means ranked best first, with
#'     standard errors and the number of locations.}
#'   \item{variance_components}{Stage-2 REML genotype-by-location variance.}
#'   \item{locations}{One row per location: stage-1 model, plots, genotypes,
#'     residual variance and df, mean SED and, for spatial fits, the
#'     effective dimension of the spatial surface.}
#'   \item{heterogeneity}{Hartley's Fmax ratio of the stage-1 residual
#'     variances with a simulated 5% critical value.}
#'   \item{stage1_means}{The stage-1 genotype means and standard errors per
#'     location.}
#'   \item{stage1_fits}{The stage-1 fits (analysis objects or `SpATS`
#'     objects), named by location.}
#'   \item{diagnostics}{Settings, number of stage-2 means, REML
#'     log-likelihood, mean SED and approximate 5% LSD.}
#' }
#'
#' @references
#' Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented p-rep
#' designs. Biometrical Journal 53, 19-27. \doi{10.1002/bimj.201000102}
#'
#' Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in two-stage
#' analysis of plant breeding trials. Crop Science 49, 1977-1988.
#' \doi{10.2135/cropsci2009.02.0083}
#'
#' Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A. and Eilers,
#' P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding
#' experiments with P-splines. Spatial Statistics 23, 52-71.
#' \doi{10.1016/j.spasta.2017.10.003}
#'
#' @seealso [analyze_augmented_prep_design()] for the one-stage analysis.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = sprintf("G%02d", 1:60),
#'   locations = c("A", "B", "C"),
#'   method = "alpha",
#'   block_size = 6,
#'   field_cols = 9,
#'   seed = 1
#' )
#' trial <- design$design
#' set.seed(2)
#' # Location C is twice as noisy as A and B.
#' trial$yield <- 50 + stats::rnorm(nrow(trial), sd = ifelse(trial$location == "C", 4, 2))
#'
#' fit <- analyze_augmented_prep_two_stage(trial, response = "yield")
#' fit
#' fit$locations
#'
#' \donttest{
#' if (requireNamespace("SpATS", quietly = TRUE)) {
#'   spatial <- analyze_augmented_prep_two_stage(trial, "yield", stage1 = "spatial")
#'   spatial$locations
#' }
#' }
#'
#' @export
analyze_augmented_prep_two_stage <- function(
    data,
    response,
    stage1 = c("blocks", "spatial"),
    weighting = c("full", "diagonal", "none"),
    genotype_by_location = c("random", "none"),
    genotype_col = "all_entries",
    location_col = "location",
    block_col = "block",
    replicate_col = "replicate",
    row_col = "row",
    col_col = "col",
    spatial_nseg = NULL,
    spatial_random_blocks = TRUE
) {
  stage1 <- match.arg(stage1)
  weighting <- match.arg(weighting)
  genotype_by_location <- match.arg(genotype_by_location)

  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) {
    data <- data$design
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or an augmented p-rep design object.", call. = FALSE)
  }
  required <- c(response, genotype_col, location_col, block_col)
  if (stage1 == "spatial") required <- c(required, row_col, col_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (stage1 == "spatial" && !requireNamespace("SpATS", quietly = TRUE)) {
    stop("`stage1 = \"spatial\"` needs the 'SpATS' package.", call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data[[response]] <- suppressWarnings(as.numeric(data[[response]]))
  locations <- unique(as.character(data[[location_col]][!is.na(data[[response]])]))
  if (length(locations) == 0) {
    stop("No complete observations are available for analysis.", call. = FALSE)
  }

  # ---- Stage 1 --------------------------------------------------------------
  stage1_results <- lapply(locations, function(loc) {
    sub <- data[as.character(data[[location_col]]) == loc, , drop = FALSE]
    if (stage1 == "blocks") {
      prep_stage1_blocks(sub, response, genotype_col, location_col, block_col, replicate_col)
    } else {
      prep_stage1_spatial(
        sub, response, genotype_col, block_col, replicate_col, row_col, col_col,
        spatial_nseg, spatial_random_blocks, loc
      )
    }
  })
  names(stage1_results) <- locations

  location_table <- dplyr::bind_rows(lapply(locations, function(loc) {
    r <- stage1_results[[loc]]
    tibble::tibble(
      location = loc,
      stage1_model = stage1,
      n_plots = r$n_plots,
      n_genotypes = length(r$means),
      residual_variance = r$residual_variance,
      residual_df = r$residual_df,
      mean_sed = prep_mean_sed(r$vcov),
      spatial_effective_dim = r$spatial_effective_dim
    )
  }))

  stage1_means <- dplyr::bind_rows(lapply(locations, function(loc) {
    r <- stage1_results[[loc]]
    tibble::tibble(
      location = loc,
      genotype = names(r$means),
      adjusted_mean = unname(r$means),
      std_error = sqrt(pmax(diag(r$vcov), 0))
    )
  }))

  heterogeneity <- prep_heterogeneity(location_table)

  # ---- Stage 2 --------------------------------------------------------------
  stage2 <- prep_stage2(stage1_results, weighting, genotype_by_location)

  counts <- table(stage1_means$genotype)
  adjusted <- tibble::tibble(
    genotype = stage2$genotypes,
    n_locations = as.integer(counts[stage2$genotypes]),
    adjusted_mean = stage2$estimate,
    std_error = stage2$std_error
  )
  adjusted <- adjusted[order(-adjusted$adjusted_mean, na.last = TRUE), , drop = FALSE]
  adjusted$rank <- ifelse(
    is.na(adjusted$adjusted_mean),
    NA_integer_,
    rank(-adjusted$adjusted_mean, ties.method = "min", na.last = "keep")
  )
  adjusted <- adjusted[, c("rank", "genotype", "n_locations", "adjusted_mean", "std_error")]

  result <- list(
    input = list(
      response = response,
      stage1 = stage1,
      weighting = weighting,
      genotype_by_location = genotype_by_location,
      genotype_col = genotype_col,
      location_col = location_col,
      block_col = block_col,
      replicate_col = replicate_col,
      row_col = row_col,
      col_col = col_col
    ),
    adjusted_means = adjusted,
    variance_components = stage2$variance_components,
    locations = location_table,
    heterogeneity = heterogeneity,
    stage1_means = stage1_means,
    stage1_fits = lapply(stage1_results, `[[`, "fit"),
    diagnostics = tibble::tibble(
      stage1 = stage1,
      weighting = weighting,
      genotype_by_location = genotype_by_location,
      n_locations = length(locations),
      n_genotypes = length(stage2$genotypes),
      n_stage2_means = stage2$n,
      stage2_df = stage2$df,
      reml_loglik = stage2$loglik,
      mean_sed = stage2$precision$mean_sed,
      lsd_5pct = stage2$precision$lsd_5pct
    )
  )
  class(result) <- c("augmented_prep_two_stage", class(result))
  result
}

#' Print a two-stage p-rep analysis
#'
#' @param x An object returned by `analyze_augmented_prep_two_stage()`.
#' @param n Number of top-ranked genotypes to show.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.augmented_prep_two_stage <- function(x, n = 10, ...) {
  cat("Two-stage augmented p-rep analysis\n")
  cat("Stage 1: ", x$input$stage1, " model per location; stage 2: ",
      x$input$weighting, " weighting, genotype x location ",
      x$input$genotype_by_location, "\n\n", sep = "")
  cat("Locations:\n")
  print(as.data.frame(x$locations[, c("location", "n_plots", "n_genotypes",
                                      "residual_variance", "residual_df", "mean_sed")]),
        row.names = FALSE, digits = 4)
  h <- x$heterogeneity
  if (!is.na(h$fmax)) {
    cat(sprintf(
      "\nHartley Fmax = %.2f (5%% critical value %.2f): %s\n",
      h$fmax, h$critical_5pct,
      if (isTRUE(h$heterogeneous)) "error variances differ between locations" else "no evidence of heterogeneity"
    ))
  }
  if (!is.null(x$variance_components)) {
    cat("\nStage-2 variance components:\n")
    print(as.data.frame(x$variance_components), row.names = FALSE, digits = 4)
  }
  cat(sprintf("\nTop adjusted means (%d genotypes):\n", nrow(x$adjusted_means)))
  print(as.data.frame(utils::head(x$adjusted_means, n)), row.names = FALSE, digits = 4)
  d <- x$diagnostics
  cat(sprintf("Mean SED = %.4g; LSD (5%%) = %.4g\n", d$mean_sed, d$lsd_5pct))
  invisible(x)
}

# ---- Stage-1 fits -----------------------------------------------------------

prep_stage1_blocks <- function(sub, response, genotype_col, location_col, block_col, replicate_col) {
  fit <- analyze_augmented_prep_design(
    sub,
    response = response,
    genotype_col = genotype_col,
    location_col = location_col,
    block_col = block_col,
    replicate_col = replicate_col
  )
  am <- fit$adjusted_means
  am <- am[!is.na(am$adjusted_mean), , drop = FALSE]
  vc <- fit$adjusted_means_vcov[am$genotype, am$genotype, drop = FALSE]
  if (fit$diagnostics$residual_df <= 0) {
    stop(
      "A location has no error degrees of freedom in the stage-1 block model; ",
      "it cannot be weighted. Use more p-rep entries per location.",
      call. = FALSE
    )
  }
  list(
    fit = fit,
    means = stats::setNames(am$adjusted_mean, am$genotype),
    vcov = vc,
    n_plots = fit$diagnostics$n_observations,
    residual_variance = fit$diagnostics$sigma^2,
    residual_df = fit$diagnostics$residual_df,
    spatial_effective_dim = NA_real_
  )
}

prep_stage1_spatial <- function(sub, response, genotype_col, block_col, replicate_col,
                                row_col, col_col, spatial_nseg, random_blocks, loc) {
  x <- data.frame(
    y = sub[[response]],
    geno = factor(as.character(sub[[genotype_col]])),
    row = suppressWarnings(as.numeric(sub[[row_col]])),
    col = suppressWarnings(as.numeric(sub[[col_col]])),
    rep = factor(if (!is.null(replicate_col) && replicate_col %in% names(sub)) sub[[replicate_col]] else 1L),
    blk = factor(paste(
      if (!is.null(replicate_col) && replicate_col %in% names(sub)) sub[[replicate_col]] else 1L,
      sub[[block_col]]
    ))
  )
  if (anyNA(x$row) || anyNA(x$col)) {
    stop(
      "Location `", loc, "` has plots without row/column coordinates; the spatial ",
      "model needs them (use `field_cols` in `augmented_prep_design()`).",
      call. = FALSE
    )
  }
  x <- x[!is.na(x$y), , drop = FALSE]
  x <- droplevels(x)
  n_rows <- length(unique(x$row))
  n_cols <- length(unique(x$col))
  if (n_rows < 4 || n_cols < 4) {
    stop(
      "Location `", loc, "` has fewer than 4 field rows or columns, too few for a ",
      "spatial surface. Use `stage1 = \"blocks\"`.",
      call. = FALSE
    )
  }
  nseg <- if (is.null(spatial_nseg)) {
    pmax(pmin(ceiling(c(n_cols, n_rows) / 2), 40L), 3L)
  } else {
    as.integer(rep_len(spatial_nseg, 2))
  }

  spatial <- stats::as.formula(sprintf("~ PSANOVA(col, row, nseg = c(%d, %d))", nseg[1], nseg[2]))
  environment(spatial) <- asNamespace("SpATS")
  fixed <- if (nlevels(x$rep) > 1) ~rep else NULL
  random <- if (isTRUE(random_blocks) && nlevels(x$blk) > nlevels(x$rep)) ~blk else NULL

  fit <- SpATS::SpATS(
    response = "y",
    genotype = "geno",
    genotype.as.random = FALSE,
    spatial = spatial,
    fixed = fixed,
    random = random,
    data = x,
    control = list(monitoring = 0)
  )
  pred <- stats::predict(fit, which = "geno", predFixed = "marginal", return.vcov.matrix = TRUE)
  vc <- as.matrix(attr(pred, "vcov"))
  genos <- as.character(pred$geno)
  dimnames(vc) <- list(genos, genos)

  eff <- fit$eff.dim
  spatial_terms <- grepl("^f\\(|row|col", names(eff))
  list(
    fit = fit,
    means = stats::setNames(pred$predicted.values, genos),
    vcov = vc,
    n_plots = nrow(x),
    residual_variance = unname(fit$psi[1]),
    residual_df = nrow(x) - sum(eff),
    spatial_effective_dim = sum(eff[spatial_terms])
  )
}

prep_mean_sed <- function(vc) {
  m <- nrow(vc)
  if (is.null(m) || m < 2) {
    return(NA_real_)
  }
  sqrt(max(2 * (m * sum(diag(vc)) - sum(vc)) / (m * (m - 1)), 0))
}

prep_heterogeneity <- function(location_table) {
  v <- location_table$residual_variance
  df <- location_table$residual_df
  ok <- !is.na(v) & v > 0 & !is.na(df) & df > 0
  if (sum(ok) < 2) {
    return(tibble::tibble(fmax = NA_real_, critical_5pct = NA_real_, heterogeneous = NA))
  }
  fmax <- max(v[ok]) / min(v[ok])
  critical <- hartley_critical_sim(df[ok], alpha = 0.05, n_sim = 20000)
  tibble::tibble(fmax = fmax, critical_5pct = critical, heterogeneous = fmax > critical)
}

# ---- Stage 2 ----------------------------------------------------------------

# Weighted REML for the stage-2 model. Each location contributes a block
# V_l = S_l + sigma2_gl * I, with S_l the stage-1 covariance (or its
# diagonal). With S_l = Q D Q', V_l^-1/2 = (D + sigma2_gl)^-1/2 Q', so the
# model is transformed to one with identity covariance and solved by least
# squares for any sigma2_gl; sigma2_gl is found by maximizing the REML
# log-likelihood.
prep_stage2 <- function(stage1_results, weighting, genotype_by_location) {
  locations <- names(stage1_results)
  genotypes <- sort(unique(unlist(lapply(stage1_results, function(r) names(r$means)))))
  stacked <- dplyr::bind_rows(lapply(locations, function(loc) {
    r <- stage1_results[[loc]]
    data.frame(location = loc, genotype = names(r$means), y = unname(r$means), stringsAsFactors = FALSE)
  }))
  stacked$location <- factor(stacked$location, levels = locations)
  stacked$genotype <- factor(stacked$genotype, levels = genotypes)
  multi <- length(locations) > 1
  x <- if (multi) {
    stats::model.matrix(~ location + genotype, stacked)
  } else {
    stats::model.matrix(~genotype, stacked)
  }
  n <- nrow(x)

  # Genotype means: intercept + average location effect + genotype effect.
  L <- matrix(0, length(genotypes), ncol(x), dimnames = list(genotypes, colnames(x)))
  L[, "(Intercept)"] <- 1
  if (multi) {
    loc_cols <- paste0("location", locations[-1])
    L[, loc_cols] <- 1 / length(locations)
  }
  g_cols <- match(paste0("genotype", genotypes), colnames(x))
  ok <- !is.na(g_cols)
  L[cbind(which(ok), g_cols[ok])] <- 1

  rows_of <- split(seq_len(n), stacked$location)
  eig <- lapply(locations, function(loc) {
    vc <- stage1_results[[loc]]$vcov
    if (weighting == "diagonal") vc <- diag(diag(vc), nrow(vc))
    e <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
    list(Q = e$vectors, D = pmax(e$values, 0))
  })
  names(eig) <- locations

  transform <- function(sigma2) {
    xt <- matrix(0, n, ncol(x))
    yt <- numeric(n)
    logdet <- 0
    for (loc in locations) {
      idx <- rows_of[[loc]]
      e <- eig[[loc]]
      w <- 1 / sqrt(e$D + sigma2)
      xt[idx, ] <- w * crossprod(e$Q, x[idx, , drop = FALSE])
      yt[idx] <- w * crossprod(e$Q, stacked$y[idx])
      logdet <- logdet + sum(log(e$D + sigma2))
    }
    list(x = xt, y = yt, logdet = logdet)
  }
  reml <- function(sigma2) {
    t <- transform(sigma2)
    q <- qr(t$x)
    r <- qr.resid(q, t$y)
    -0.5 * (t$logdet + 2 * sum(log(abs(diag(q$qr)[seq_len(q$rank)]))) + sum(r^2))
  }

  if (weighting == "none") {
    q <- qr(x)
    coef <- qr.coef(q, stacked$y)
    df <- n - q$rank
    sigma2 <- if (df > 0) sum(qr.resid(q, stacked$y)^2) / df else NA_real_
    est <- arc_linear_estimates(q, coef, L)
    scale <- sigma2
    gxl <- NA_real_
    loglik <- NA_real_
    vcomp <- tibble::tibble(
      component = "residual (genotype x location + stage-1 error)",
      variance = sigma2
    )
  } else {
    gxl <- 0
    if (genotype_by_location == "random" && multi) {
      upper <- max(stats::var(stacked$y) * 2, 1e-8)
      opt <- stats::optimize(reml, c(0, upper), maximum = TRUE)
      gxl <- if (reml(0) >= opt$objective) 0 else opt$maximum
    }
    t <- transform(gxl)
    q <- qr(t$x)
    coef <- qr.coef(q, t$y)
    est <- arc_linear_estimates(q, coef, L)
    df <- n - q$rank
    scale <- 1
    loglik <- reml(gxl)
    vcomp <- if (genotype_by_location == "random" && multi) {
      tibble::tibble(component = "genotype x location", variance = gxl)
    } else {
      NULL
    }
  }

  estimate <- est$estimate
  std_error <- sqrt(scale * colSums(est$W^2))
  estimate[!est$estimable] <- NA_real_
  std_error[!est$estimable] <- NA_real_
  precision <- arc_mean_precision(
    est$W[, est$estimable, drop = FALSE],
    rep(0, sum(est$estimable)),
    scale,
    max(df, 0)
  )

  list(
    genotypes = genotypes,
    estimate = unname(estimate),
    std_error = unname(std_error),
    variance_components = vcomp,
    n = n,
    df = df,
    loglik = loglik,
    precision = precision
  )
}
