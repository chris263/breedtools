#' Two-stage analysis of a multi-location trial, for any design
#'
#' Analyses each location on its own (stage 1) and combines the genotype
#' means across locations (stage 2), weighting them by their precision (Möhring
#' and Piepho, 2009). Works for every design of the package grown at several
#' locations: augmented row-column, augmented p-rep and resolvable
#' row-column designs. For augmented p-rep designs,
#' [analyze_augmented_prep_two_stage()] gives the same analysis with the
#' p-rep block model of the paper.
#'
#' @section Stage 1:
#' \describe{
#'   \item{`stage1 = "spatial"` (default)}{A 'SpATS' spatial model at every
#'     location (Rodríguez-Álvarez et al., 2018): genotypes fixed, a smooth
#'     two-dimensional P-spline surface over the field rows and columns,
#'     replicates fixed and blocks random when `replicate_col` and
#'     `block_col` are given.}
#'   \item{`stage1 = "blocks"`}{A linear model at every location that needs no
#'     design object: replicates and blocks within replicates fixed (when
#'     `replicate_col` and `block_col` are given), genotypes fixed. For data
#'     read from Breedbase without field coordinates.}
#'   \item{`stage1 = "design"`}{The design's own model at every location:
#'     [analyze_augmented_row_column_design()] for augmented row-column
#'     designs, [analyze_resolvable_row_column()] (row-column model) for
#'     resolvable row-column designs, and the block model of
#'     [analyze_augmented_prep_design()] for p-rep designs. The design is
#'     recognised from its columns.}
#' }
#' The stage-1 genotype means are BLUEs with their full covariance matrix.
#'
#' @section Stage 2:
#' As in [analyze_augmented_prep_two_stage()]: fixed locations and genotypes,
#' random genotype-by-location effects (REML), and stage-1 errors with the
#' stage-1 covariance (`weighting = "full"`), only the stage-1 variances
#' (`"diagonal"`) or a common error (`"none"`, for comparison).
#'
#' @param data A data frame with all locations (one row per plot), with the
#'   design columns and a response column.
#' @param response Name of the numeric response column.
#' @param stage1 `"spatial"`, `"blocks"` or `"design"`.
#' @param weighting `"full"`, `"diagonal"` or `"none"`.
#' @param genotype_by_location `"random"` or `"none"`.
#' @param genotype_col,location_col Genotype and location columns.
#' @param row_col,col_col Field row and column coordinates (spatial stage 1).
#' @param replicate_col,block_col Optional replicate (fixed) and block
#'   (random) columns of the spatial stage 1. By default `"replicate"` and
#'   `"block"` when the data have them.
#' @param spatial_nseg Optional number of P-spline segments, `c(col, row)`.
#'
#' @return An object of class `two_stage_analysis`: a list with
#'   `adjusted_means` (combined, ranked), `variance_components`,
#'   `locations` (one row per location), `heterogeneity` (Hartley's Fmax),
#'   `stage1_means`, `stage1_fits` (per location, with `plot_results`, for
#'   [first_stage_values()] and [compare_first_stage()]) and `diagnostics`.
#'
#' @references
#' Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in two-stage
#' analysis of plant breeding trials. Crop Science 49, 1977-1988.
#' \doi{10.2135/cropsci2009.02.0083}
#'
#' Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A. and Eilers,
#' P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding
#' experiments with P-splines. Spatial Statistics 23, 52-71.
#' \doi{10.1016/j.spasta.2017.10.003}
#'
#' @examples
#' # The same augmented row-column layout, randomized at three locations.
#' sites <- c("North", "Central", "South")
#' trials <- lapply(seq_along(sites), function(i) {
#'   d <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                    12, 12, 4, 4, n_candidates = 5, seed = i)$design
#'   d$location <- sites[i]
#'   d
#' })
#' met <- do.call(rbind, trials)
#' set.seed(1)
#' g <- stats::setNames(stats::rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
#' met$yield <- 50 + g[met$all_entries] + 0.5 * met$row + stats::rnorm(nrow(met))
#' fit <- analyze_two_stage(met, "yield", stage1 = "design")
#' fit
#' \donttest{
#' if (requireNamespace("SpATS", quietly = TRUE)) {
#'   spatial <- analyze_two_stage(met, "yield")      # SpATS at every location
#'   spatial$locations
#' }
#' }
#'
#' @export
analyze_two_stage <- function(data, response, stage1 = c("spatial", "blocks", "design"),
                              weighting = c("full", "diagonal", "none"),
                              genotype_by_location = c("random", "none"),
                              genotype_col = "all_entries", location_col = "location",
                              row_col = "row", col_col = "col",
                              replicate_col = NULL, block_col = NULL, spatial_nseg = NULL) {
  stage1 <- match.arg(stage1)
  weighting <- match.arg(weighting)
  genotype_by_location <- match.arg(genotype_by_location)
  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) data <- data$design
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  required <- c(response, genotype_col, location_col, if (stage1 == "spatial") c(row_col, col_col))
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("Missing required column(s): ", paste(absent, collapse = ", "), call. = FALSE)
  if (stage1 == "spatial" && !requireNamespace("SpATS", quietly = TRUE)) {
    stop("`stage1 = \"spatial\"` needs the 'SpATS' package (or use stage1 = \"design\").", call. = FALSE)
  }
  if (is.null(replicate_col) && "replicate" %in% names(data)) replicate_col <- "replicate"
  if (is.null(block_col) && "block" %in% names(data)) block_col <- "block"
  if ("is_filler" %in% names(data)) data <- data[data$is_filler != 1, , drop = FALSE]
  data[[response]] <- suppressWarnings(as.numeric(data[[response]]))
  locations <- unique(as.character(data[[location_col]][!is.na(data[[response]])]))
  if (!length(locations)) stop("No complete observations are available for analysis.", call. = FALSE)

  results <- lapply(locations, function(loc) {
    sub <- data[as.character(data[[location_col]]) == loc, , drop = FALSE]
    if (stage1 == "spatial") {
      ts_stage1_spatial(sub, response, genotype_col, location_col, row_col, col_col, replicate_col,
                        block_col, spatial_nseg, loc)
    } else if (stage1 == "blocks") {
      ts_stage1_blocks(sub, response, genotype_col, replicate_col, block_col, loc)
    } else {
      ts_stage1_design(sub, response, genotype_col, location_col, loc)
    }
  })
  names(results) <- locations

  location_table <- dplyr::bind_rows(lapply(locations, function(loc) {
    r <- results[[loc]]
    tibble::tibble(location = loc, stage1_model = r$model_label, n_plots = r$n_plots,
                   n_genotypes = length(r$means), residual_variance = r$residual_variance,
                   residual_df = r$residual_df, mean_sed = prep_mean_sed(r$vcov),
                   spatial_effective_dim = r$spatial_effective_dim)
  }))
  stage1_means <- dplyr::bind_rows(lapply(locations, function(loc) {
    r <- results[[loc]]
    tibble::tibble(location = loc, genotype = names(r$means), adjusted_mean = unname(r$means),
                   std_error = sqrt(pmax(diag(r$vcov), 0)))
  }))
  stage2 <- prep_stage2(results, weighting, genotype_by_location)
  counts <- table(stage1_means$genotype)
  adjusted <- tibble::tibble(genotype = stage2$genotypes, n_locations = as.integer(counts[stage2$genotypes]),
                             adjusted_mean = stage2$estimate, std_error = stage2$std_error)
  adjusted <- adjusted[order(-adjusted$adjusted_mean, na.last = TRUE), , drop = FALSE]
  adjusted$rank <- ifelse(is.na(adjusted$adjusted_mean), NA_integer_,
                          rank(-adjusted$adjusted_mean, ties.method = "min", na.last = "keep"))
  adjusted <- adjusted[, c("rank", "genotype", "n_locations", "adjusted_mean", "std_error")]

  result <- list(
    input = list(response = response, stage1 = stage1, weighting = weighting,
                 genotype_by_location = genotype_by_location, genotype_col = genotype_col,
                 location_col = location_col, row_col = row_col, col_col = col_col,
                 replicate_col = replicate_col, block_col = block_col),
    adjusted_means = adjusted,
    variance_components = stage2$variance_components,
    locations = location_table,
    heterogeneity = prep_heterogeneity(location_table),
    stage1_means = stage1_means,
    stage1_fits = lapply(results, `[[`, "fit"),
    diagnostics = tibble::tibble(
      stage1 = stage1, weighting = weighting, genotype_by_location = genotype_by_location,
      n_locations = length(locations), n_genotypes = length(stage2$genotypes), n_stage2_means = stage2$n,
      stage2_df = stage2$df, reml_loglik = stage2$loglik, mean_sed = stage2$precision$mean_sed,
      lsd_5pct = stage2$precision$lsd_5pct
    )
  )
  class(result) <- c("two_stage_analysis", class(result))
  result
}

#' Print a two-stage analysis
#'
#' @param x An object returned by [analyze_two_stage()].
#' @param n Number of top-ranked genotypes to show.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.two_stage_analysis <- function(x, n = 10, ...) {
  cat("Two-stage analysis\n")
  cat("Stage 1: ", x$input$stage1, " model per location; stage 2: ", x$input$weighting,
      " weighting, genotype x location ", x$input$genotype_by_location, "\n\n", sep = "")
  cat("Locations:\n")
  print(as.data.frame(x$locations[, c("location", "stage1_model", "n_plots", "n_genotypes",
                                      "residual_variance", "residual_df", "mean_sed")]),
        row.names = FALSE, digits = 4)
  h <- x$heterogeneity
  if (!is.na(h$fmax)) {
    cat(sprintf("\nHartley Fmax = %.2f (5%% critical value %.2f): %s\n", h$fmax, h$critical_5pct,
                if (isTRUE(h$heterogeneous)) "error variances differ between locations" else "no evidence of heterogeneity"))
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

# SpATS at one location; the fit is returned as a small analysis object with
# plot_results (fitted values and residuals) and adjusted means.
ts_stage1_spatial <- function(sub, response, genotype_col, location_col, row_col, col_col, replicate_col,
                              block_col, spatial_nseg, loc) {
  r <- prep_stage1_spatial(sub, response, genotype_col, block_col, replicate_col, row_col, col_col,
                           spatial_nseg, random_blocks = !is.null(block_col), loc)
  sp <- r$fit
  pr <- sub[!is.na(sub[[response]]), , drop = FALSE]
  if (!"plots" %in% names(pr)) pr$plots <- seq_len(nrow(pr))
  i <- match(pr$plots, sp$data$plots)
  pr$fitted <- unname(sp$fitted[i])
  pr$residual <- pr[[response]] - pr$fitted
  r$fit <- structure(list(
    input = list(response = response, genotype_col = genotype_col, spatial = TRUE),
    plot_results = tibble::as_tibble(pr),
    adjusted_means = tibble::tibble(genotype = names(r$means), adjusted_mean = unname(r$means),
                                    std_error = sqrt(pmax(diag(r$vcov), 0))),
    adjusted_means_vcov = r$vcov,
    diagnostics = tibble::tibble(model_type = "spatial", n_observations = r$n_plots,
                                 residual_variance = r$residual_variance,
                                 spatial_effective_dim = r$spatial_effective_dim),
    model = sp
  ), class = c("stage1_analysis", "list"))
  r$model_label <- "spatial (SpATS)"
  r
}

# A linear model at one location without a design object: replicates and
# blocks within replicates fixed, genotypes fixed; genotype means are
# estimable functions (intercept + average design effect + genotype effect)
# with covariance sigma^2 W'W.
ts_stage1_blocks <- function(sub, response, genotype_col, replicate_col, block_col, loc) {
  d <- sub[!is.na(suppressWarnings(as.numeric(sub[[response]]))) & !is.na(sub[[genotype_col]]), , drop = FALSE]
  d$.y <- as.numeric(d[[response]])
  d$.g <- factor(as.character(d[[genotype_col]]))
  terms <- character(0)
  if (!is.null(replicate_col) && replicate_col %in% names(d) && length(unique(d[[replicate_col]])) > 1) {
    d$.rep <- factor(d[[replicate_col]])
    terms <- c(terms, ".rep")
  }
  if (!is.null(block_col) && block_col %in% names(d) && length(unique(d[[block_col]])) > 1) {
    d$.blk <- factor(paste(if (".rep" %in% terms) d$.rep else 1, d[[block_col]]))
    terms <- c(terms, ".blk")
  }
  form <- stats::as.formula(paste(".y ~", paste(c(terms, ".g"), collapse = " + ")))
  x <- stats::model.matrix(form, d)
  qx <- qr(x)
  coef <- qr.coef(qx, d$.y)
  fitted <- drop(x[, qx$pivot[seq_len(qx$rank)], drop = FALSE] %*% coef[qx$pivot[seq_len(qx$rank)]])
  df <- nrow(d) - qx$rank
  if (df <= 0) stop("Location `", loc, "` has no error degrees of freedom in the block model.", call. = FALSE)
  sigma2 <- sum((d$.y - fitted)^2) / df
  genos <- levels(d$.g)
  g_cols <- match(paste0(".g", genos), colnames(x))
  L <- matrix(colMeans(x), length(genos), ncol(x), byrow = TRUE)
  L[, grepl("^\\.g", colnames(x))] <- 0
  L[cbind(which(!is.na(g_cols)), g_cols[!is.na(g_cols)])] <- 1
  est <- arc_linear_estimates(qx, coef, L)
  ok <- which(est$estimable)
  vc <- sigma2 * crossprod(est$W[, ok, drop = FALSE])
  dimnames(vc) <- list(genos[ok], genos[ok])
  pr <- d[, setdiff(names(d), c(".y", ".g", ".rep", ".blk")), drop = FALSE]
  if (!"plots" %in% names(pr)) pr$plots <- seq_len(nrow(pr))
  pr$fitted <- fitted
  pr$residual <- d$.y - fitted
  means <- stats::setNames(est$estimate[ok], genos[ok])
  fit <- structure(list(
    input = list(response = response, genotype_col = genotype_col, spatial = FALSE),
    plot_results = tibble::as_tibble(pr),
    adjusted_means = tibble::tibble(genotype = names(means), adjusted_mean = unname(means),
                                    std_error = sqrt(pmax(diag(vc), 0))),
    adjusted_means_vcov = vc,
    diagnostics = tibble::tibble(model_type = "blocks", n_observations = nrow(d), residual_variance = sigma2,
                                 residual_df = df)
  ), class = c("stage1_analysis", "list"))
  list(fit = fit, means = means, vcov = vc, n_plots = nrow(d), residual_variance = sigma2, residual_df = df,
       spatial_effective_dim = NA_real_, model_label = "blocks (linear model)")
}

# The design's own model at one location.
ts_stage1_design <- function(sub, response, genotype_col, location_col, loc) {
  cols <- names(sub)
  if (all(c("rowgroup", "colgroup") %in% cols)) {
    fit <- analyze_augmented_row_column_design(sub, response, genotype_col = genotype_col)
    label <- "augmented row-column"
    sigma2 <- fit$diagnostics$sigma^2
    df <- fit$diagnostics$residual_df
  } else if (all(c("row_in_rep", "col_in_rep") %in% cols)) {
    fit <- analyze_resolvable_row_column(sub, response, model = "rc", compare = FALSE, genotype_col = genotype_col)
    label <- "resolvable row-column (RC)"
    vp <- fit$variance_parameters
    sigma2 <- vp$estimate[vp$parameter == "sigma2"]
    df <- fit$diagnostics$intrablock_error_df
  } else if ("block" %in% cols) {
    r <- prep_stage1_blocks(sub, response, genotype_col, location_col, "block", "replicate")
    r$model_label <- "p-rep blocks"
    return(r)
  } else {
    stop("Location `", loc, "`: the design was not recognised from its columns; use stage1 = \"spatial\".",
         call. = FALSE)
  }
  am <- fit$adjusted_means
  am <- am[!is.na(am$adjusted_mean), , drop = FALSE]
  vc <- fit$adjusted_means_vcov
  if (is.null(vc)) stop("Location `", loc, "`: the stage-1 analysis gave no covariance matrix.", call. = FALSE)
  keep <- intersect(as.character(am$genotype), rownames(vc))
  if (is.na(df) || df <= 0) {
    stop("Location `", loc, "` has no error degrees of freedom in the stage-1 model; it cannot be weighted.",
         call. = FALSE)
  }
  list(fit = fit, means = stats::setNames(am$adjusted_mean, am$genotype)[keep], vcov = vc[keep, keep, drop = FALSE],
       n_plots = sum(!is.na(sub[[response]])), residual_variance = sigma2, residual_df = df,
       spatial_effective_dim = NA_real_, model_label = label)
}
