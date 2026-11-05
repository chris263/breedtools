#' Single-stage analysis of several trials
#'
#' Fits all trials (locations, environments) together in one mixed model and
#' gives the genotype adjusted means across trials, to compare with (or use
#' instead of) the two-stage analysis of [analyze_two_stage()]:
#' \deqn{y = \mu + \text{trial} + \text{genotype} + \text{rep(trial)} +
#'   \text{block(rep, trial)} + \text{genotype} \times \text{trial} + e,}
#' with trials and genotypes fixed, and replicates, blocks and the
#' genotype-by-trial interaction random (fitted with 'lme4'). Genotype means
#' average the trial effects with equal weights. Without 'lme4' the model is
#' fitted with all effects fixed and no genotype-by-trial term.
#'
#' With `spatial = TRUE` the model is fitted with 'SpATS': genotypes and
#' trials fixed, replicates within trials fixed, blocks and the
#' genotype-by-trial interaction random, and a smooth P-spline surface over
#' the field rows and columns. SpATS fits one surface per model, so the
#' trials are placed side by side in one coordinate system, separated by
#' empty columns: each trial keeps its own spatial trend and the trends do not
#' run from one trial into another.
#'
#' One error variance is shared by all trials here; when trials differ much
#' in precision, the weighted two-stage analysis is the better choice
#' (Möhring and Piepho, 2009).
#'
#' @param data A data frame with one row per plot, e.g. from
#'   [breedbase_phenotypes()].
#' @param response Name of the numeric response column.
#' @param genotype_col,location_col Genotype and trial (location) columns.
#' @param replicate_col,block_col Optional replicate and block columns; by
#'   default `"rep"`/`"replicate"` and `"block"` when the data have them.
#' @param genotype_by_location `"random"` (default) or `"none"`.
#' @param spatial Logical. Fit the single-stage model with 'SpATS' (needs the
#'   field `row` and `col` of every plot).
#' @param row_col,col_col Field row and column columns (for `spatial = TRUE`).
#' @param spatial_nseg Optional number of P-spline segments, `c(col, row)`.
#'
#' @return An object of class `single_stage_analysis`: `adjusted_means`
#'   (ranked, with standard errors and the number of trials),
#'   `adjusted_means_vcov`, `variance_components`, `plot_results` (fitted
#'   values and residuals), `diagnostics` and the fitted `model`.
#'
#' @references
#' Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in two-stage
#' analysis of plant breeding trials. Crop Science 49, 1977-1988.
#' \doi{10.2135/cropsci2009.02.0083}
#'
#' @examples
#' sites <- c("North", "South")
#' met <- do.call(rbind, lapply(1:2, function(i) {
#'   d <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
#'                                    n_candidates = 5, seed = i)$design
#'   d$location <- sites[i]
#'   d
#' }))
#' set.seed(1)
#' g <- stats::setNames(stats::rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
#' met$yield <- 50 + g[met$all_entries] + stats::rnorm(nrow(met))
#' one <- analyze_single_stage(met, "yield", genotype_col = "all_entries", location_col = "location")
#' one
#'
#' @export
analyze_single_stage <- function(data, response, genotype_col = "accession", location_col = "trial_name",
                                 replicate_col = NULL, block_col = NULL,
                                 genotype_by_location = c("random", "none"), spatial = FALSE,
                                 row_col = "row", col_col = "col", spatial_nseg = NULL) {
  genotype_by_location <- match.arg(genotype_by_location)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  d <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  absent <- setdiff(c(response, genotype_col, location_col), names(d))
  if (length(absent)) stop("Missing required column(s): ", paste(absent, collapse = ", "), call. = FALSE)
  if (is.null(replicate_col)) replicate_col <- intersect(c("rep", "replicate"), names(d))[1]
  if (is.null(block_col) && "block" %in% names(d)) block_col <- "block"
  if (length(replicate_col) && is.na(replicate_col)) replicate_col <- NULL
  d$.y <- suppressWarnings(as.numeric(d[[response]]))
  d <- d[!is.na(d$.y) & !is.na(d[[genotype_col]]) & !is.na(d[[location_col]]), , drop = FALSE]
  if (nrow(d) < 3) stop("Not enough observations.", call. = FALSE)
  d$.g <- factor(as.character(d[[genotype_col]]))
  d$.loc <- factor(as.character(d[[location_col]]), levels = unique(as.character(d[[location_col]])))
  multi <- nlevels(d$.loc) > 1
  has_rep <- !is.null(replicate_col) && replicate_col %in% names(d) && length(unique(d[[replicate_col]])) > 1
  has_blk <- !is.null(block_col) && block_col %in% names(d) && length(unique(d[[block_col]])) > 1
  d$.rep <- factor(paste(d$.loc, if (has_rep) d[[replicate_col]] else 1))
  d$.blk <- factor(paste(d$.rep, if (has_blk) d[[block_col]] else 1))
  if (isTRUE(spatial)) {
    return(ss_spatial(d, response, genotype_col, location_col, replicate_col, block_col, has_rep, has_blk, multi,
                      genotype_by_location, row_col, col_col, spatial_nseg))
  }
  fixed <- if (multi) ".y ~ .loc + .g" else ".y ~ .g"
  use_lme4 <- requireNamespace("lme4", quietly = TRUE)
  random <- c(if (has_rep && nlevels(d$.rep) > nlevels(d$.loc)) "(1 | .rep)",
              if (has_blk) "(1 | .blk)",
              if (multi && genotype_by_location == "random") "(1 | .g:.loc)")

  if (use_lme4 && length(random)) {
    model <- lme4::lmer(stats::as.formula(paste(fixed, "+", paste(random, collapse = " + "))), data = d, REML = TRUE)
    beta <- lme4::fixef(model)
    vb <- as.matrix(stats::vcov(model))
    fitted <- unname(stats::fitted(model))
    vc <- as.data.frame(lme4::VarCorr(model))
    labels <- c(.rep = "replicates within trials", .blk = "blocks within replicates",
                `.g:.loc` = "genotype x trial", Residual = "residual")
    variance_components <- tibble::tibble(component = unname(ifelse(vc$grp %in% names(labels), labels[vc$grp], vc$grp)),
                                          variance = vc$vcov)
    model_label <- "mixed (lme4)"
  } else {
    terms <- c(if (multi) ".loc", if (has_rep && nlevels(d$.rep) > nlevels(d$.loc)) ".rep", if (has_blk) ".blk", ".g")
    model <- stats::lm(stats::as.formula(paste(".y ~", paste(terms, collapse = " + "))), data = d)
    beta <- stats::coef(model)
    beta <- beta[!is.na(beta)]
    vb <- stats::vcov(model)[names(beta), names(beta), drop = FALSE]
    fitted <- unname(stats::fitted(model))
    variance_components <- tibble::tibble(component = "residual", variance = summary(model)$sigma^2)
    model_label <- "fixed effects (lm)"
  }
  # Genotype means: intercept + average trial effect + genotype effect.
  genos <- levels(d$.g)
  L <- matrix(0, length(genos), length(beta), dimnames = list(genos, names(beta)))
  L[, "(Intercept)"] <- 1
  loc_cols <- grep("^\\.loc", names(beta), value = TRUE)
  if (length(loc_cols)) L[, loc_cols] <- 1 / nlevels(d$.loc)
  other <- setdiff(names(beta), c("(Intercept)", loc_cols, grep("^\\.g", names(beta), value = TRUE)))
  if (length(other)) {
    x <- stats::model.matrix(model)
    L[, other] <- matrix(colMeans(x[, other, drop = FALSE]), length(genos), length(other), byrow = TRUE)
  }
  g_cols <- match(paste0(".g", genos), names(beta))
  base <- genos[is.na(g_cols)]
  L[cbind(which(!is.na(g_cols)), g_cols[!is.na(g_cols)])] <- 1
  # A genotype dropped from the fixed effects (not the reference level) is not estimable.
  estimable <- !is.na(g_cols) | genos == levels(d$.g)[1]
  est <- drop(L %*% beta)
  cov <- L %*% vb %*% t(L)
  est[!estimable] <- NA_real_
  ok <- which(estimable)
  vcov_means <- as.matrix(cov[ok, ok, drop = FALSE])
  n_loc <- tapply(d$.loc, d$.g, function(x) length(unique(x)))
  adjusted <- tibble::tibble(genotype = genos, n_locations = as.integer(n_loc[genos]), adjusted_mean = unname(est),
                             std_error = sqrt(pmax(diag(as.matrix(cov)), 0)))
  adjusted$std_error[!estimable] <- NA_real_
  adjusted <- adjusted[order(-adjusted$adjusted_mean, na.last = TRUE), , drop = FALSE]
  adjusted$rank <- ifelse(is.na(adjusted$adjusted_mean), NA_integer_,
                          rank(-adjusted$adjusted_mean, ties.method = "min", na.last = "keep"))
  adjusted <- adjusted[, c("rank", "genotype", "n_locations", "adjusted_mean", "std_error")]
  plot_results <- tibble::as_tibble(d[, setdiff(names(d), c(".y", ".g", ".loc", ".rep", ".blk")), drop = FALSE])
  plot_results$fitted <- fitted
  plot_results$residual <- d$.y - fitted
  mean_sed <- prep_mean_sed(vcov_means)
  result <- list(
    input = list(response = response, genotype_col = genotype_col, location_col = location_col,
                 replicate_col = if (has_rep) replicate_col else NULL, block_col = if (has_blk) block_col else NULL,
                 genotype_by_location = genotype_by_location),
    adjusted_means = adjusted,
    adjusted_means_vcov = vcov_means,
    variance_components = variance_components,
    plot_results = plot_results,
    diagnostics = tibble::tibble(model = model_label, n_observations = nrow(d), n_locations = nlevels(d$.loc),
                                 n_genotypes = length(genos), mean_sed = mean_sed),
    model = model
  )
  class(result) <- c("single_stage_analysis", "list")
  result
}

#' Print a single-stage analysis
#'
#' @param x An object returned by [analyze_single_stage()].
#' @param n Number of top-ranked genotypes to show.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.single_stage_analysis <- function(x, n = 10, ...) {
  d <- x$diagnostics
  cat("Single-stage analysis (", d$model, "): ", d$n_observations, " plots, ", d$n_locations, " trials, ",
      d$n_genotypes, " genotypes\n\n", sep = "")
  cat("Variance components:\n")
  print(as.data.frame(x$variance_components), row.names = FALSE, digits = 4)
  cat(sprintf("\nTop adjusted means (%d genotypes):\n", nrow(x$adjusted_means)))
  print(as.data.frame(utils::head(x$adjusted_means, n)), row.names = FALSE, digits = 4)
  cat(sprintf("Mean SED = %.4g\n", d$mean_sed))
  invisible(x)
}

# Single stage with SpATS: trials side by side in one coordinate system,
# separated by empty columns; trial and replicate(trial) fixed; block and
# genotype x trial random; genotype fixed; one P-spline surface.
ss_spatial <- function(d, response, genotype_col, location_col, replicate_col, block_col, has_rep, has_blk, multi,
                       genotype_by_location, row_col, col_col, spatial_nseg) {
  if (!requireNamespace("SpATS", quietly = TRUE)) stop("`spatial = TRUE` needs the 'SpATS' package.", call. = FALSE)
  absent <- setdiff(c(row_col, col_col), names(d))
  if (length(absent)) stop("The spatial model needs the field coordinates: ", paste(absent, collapse = ", "), call. = FALSE)
  d$.row <- suppressWarnings(as.numeric(d[[row_col]]))
  d$.col0 <- suppressWarnings(as.numeric(d[[col_col]]))
  if (anyNA(d$.row) || anyNA(d$.col0)) stop("Some plots have no row or column.", call. = FALSE)
  # Each trial gets its own block of columns, with a gap as wide as the
  # widest trial between trials.
  width <- max(tapply(d$.col0, d$.loc, function(x) diff(range(x)) + 1))
  start <- (seq_len(nlevels(d$.loc)) - 1) * 2 * width
  d$.col <- d$.col0 - stats::ave(d$.col0, d$.loc, FUN = min) + 1 + start[as.integer(d$.loc)]
  nseg <- if (is.null(spatial_nseg)) {
    c(min(max(ceiling(length(unique(d$.col)) / 2), 3L), 60L), min(max(ceiling(length(unique(d$.row)) / 2), 3L), 40L))
  } else {
    as.integer(rep_len(spatial_nseg, 2))
  }
  spatial <- stats::as.formula(sprintf("~ PSANOVA(.col, .row, nseg = c(%d, %d))", nseg[1], nseg[2]))
  environment(spatial) <- asNamespace("SpATS")
  # Replicates within trials contain the trial effect: one of the two only
  # (SpATS needs a fixed part of full rank).
  fixed_terms <- if (has_rep && nlevels(d$.rep) > nlevels(d$.loc)) ".rep" else if (multi) ".loc" else character(0)
  d$.gxl <- factor(paste(d$.g, d$.loc))
  random_terms <- c(if (has_blk) ".blk", if (multi && genotype_by_location == "random") ".gxl")
  fit <- SpATS::SpATS(
    response = ".y", genotype = ".g", genotype.as.random = FALSE, spatial = spatial,
    fixed = if (length(fixed_terms)) stats::as.formula(paste("~", paste(fixed_terms, collapse = " + "))) else NULL,
    random = if (length(random_terms)) stats::as.formula(paste("~", paste(random_terms, collapse = " + "))) else NULL,
    data = d, control = list(monitoring = 0)
  )
  pred <- stats::predict(fit, which = ".g", predFixed = "marginal", return.vcov.matrix = TRUE)
  vc <- as.matrix(attr(pred, "vcov"))
  genos <- as.character(pred$.g)
  dimnames(vc) <- list(genos, genos)
  n_loc <- tapply(d$.loc, d$.g, function(x) length(unique(x)))
  adjusted <- tibble::tibble(genotype = genos, n_locations = as.integer(n_loc[genos]),
                             adjusted_mean = pred$predicted.values, std_error = pred$standard.errors)
  adjusted <- adjusted[order(-adjusted$adjusted_mean), , drop = FALSE]
  adjusted$rank <- rank(-adjusted$adjusted_mean, ties.method = "min")
  adjusted <- adjusted[, c("rank", "genotype", "n_locations", "adjusted_mean", "std_error")]
  plot_results <- tibble::as_tibble(d[, setdiff(names(d), c(".y", ".g", ".loc", ".rep", ".blk", ".row", ".col0",
                                                         ".col", ".gxl")), drop = FALSE])
  plot_results$fitted <- unname(fit$fitted)
  plot_results$residual <- d$.y - plot_results$fitted
  eff <- fit$eff.dim
  labels <- c(.blk = "blocks within replicates", .gxl = "genotype x trial")
  comps <- names(fit$var.comp)
  variance_components <- tibble::tibble(
    component = c(unname(ifelse(comps %in% names(labels), labels[comps], comps)), "residual"),
    variance = c(unname(fit$var.comp), unname(fit$psi[1]))
  )
  variance_components <- variance_components[!grepl("f\\(", variance_components$component), , drop = FALSE]
  result <- list(
    input = list(response = response, genotype_col = genotype_col, location_col = location_col,
                 replicate_col = if (has_rep) replicate_col else NULL, block_col = if (has_blk) block_col else NULL,
                 genotype_by_location = genotype_by_location, spatial = TRUE),
    adjusted_means = adjusted,
    adjusted_means_vcov = vc,
    variance_components = variance_components,
    plot_results = plot_results,
    diagnostics = tibble::tibble(model = "spatial (SpATS)", n_observations = nrow(d), n_locations = nlevels(d$.loc),
                                 n_genotypes = length(genos), mean_sed = prep_mean_sed(vc),
                                 spatial_effective_dim = sum(eff[grepl("^f\\(", names(eff))])),
    model = fit
  )
  class(result) <- c("single_stage_analysis", "list")
  result
}
