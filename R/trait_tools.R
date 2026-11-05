# Looking at trait data before and after the first stage: plots, outliers,
# and raw versus adjusted values.

#' Find outliers in trait data (IQR rule)
#'
#' Flags values below \eqn{Q_1 - k \cdot IQR} or above
#' \eqn{Q_3 + k \cdot IQR}, where \eqn{Q_1}, \eqn{Q_3} are the quartiles and
#' \eqn{k} is `multiplier` (Tukey's fences; 1.5 by default, 3 for only the
#' extreme values). The rule is applied within each group of `by`, usually
#' each location, because locations have different means.
#'
#' @param data A data frame with one row per plot.
#' @param traits Names of the numeric trait columns.
#' @param multiplier The IQR multiplier \eqn{k}.
#' @param by Optional grouping column(s); by default `"location"` when the
#'   data have it.
#' @param id_cols Columns that identify the plot in the report (those
#'   present in `data`).
#'
#' @return A tibble with one row per outlier: the identifying columns, the
#'   group, `trait`, `value`, the fences `lower` and `upper`, and `side`
#'   (`"low"` or `"high"`).
#'
#' @examples
#' d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
#'                 yield = c(rnorm(18, 50, 2), 80, 10))
#' find_outliers(d, "yield")
#'
#' @export
find_outliers <- function(data, traits, multiplier = 1.5,
                          by = if ("location" %in% names(data)) "location" else NULL,
                          id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries", "accession",
                                      "row", "col")) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.numeric(multiplier) || length(multiplier) != 1 || multiplier <= 0) {
    stop("`multiplier` must be one positive number, e.g. 1.5.", call. = FALSE)
  }
  absent <- setdiff(c(traits, by), names(data))
  if (length(absent)) stop("Column(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  id_cols <- intersect(unique(c(by, id_cols)), names(data))
  group <- if (length(by)) interaction(data[by], drop = TRUE, lex.order = TRUE) else factor(rep("all", nrow(data)))
  rows <- list()
  for (t in traits) {
    x <- suppressWarnings(as.numeric(data[[t]]))
    for (g in levels(group)) {
      i <- which(group == g)
      v <- x[i]
      if (sum(!is.na(v)) < 4) next
      q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      iqr <- q[2] - q[1]
      lower <- q[1] - multiplier * iqr
      upper <- q[2] + multiplier * iqr
      out <- i[!is.na(v) & (v < lower | v > upper)]
      if (!length(out)) next
      part <- tibble::as_tibble(data[out, id_cols, drop = FALSE])
      part$trait <- t
      part$value <- x[out]
      part$lower <- lower
      part$upper <- upper
      part$side <- ifelse(x[out] < lower, "low", "high")
      part$.row <- out
      rows[[length(rows) + 1]] <- part
    }
  }
  res <- if (length(rows)) dplyr::bind_rows(rows) else
    tibble::tibble(trait = character(), value = numeric(), lower = numeric(), upper = numeric(),
                   side = character(), .row = integer())
  attr(res, "multiplier") <- multiplier
  res
}

#' Remove outliers from trait data (IQR rule)
#'
#' Finds outliers with [find_outliers()] and sets those values to `NA`, so
#' the plots stay in the design but the values are left out of the
#' analysis. The plots and values removed are printed and returned with the
#' data.
#'
#' @inheritParams find_outliers
#' @param quiet Logical. If `FALSE` (default), print what was removed.
#'
#' @return `data` with the outlying values set to `NA`. The removed values
#'   are in `attr(result, "outliers")` (as from [find_outliers()]).
#'
#' @examples
#' d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
#'                 yield = c(rnorm(18, 50, 2), 80, 10))
#' clean <- remove_outliers(d, "yield")
#' attr(clean, "outliers")
#'
#' @export
remove_outliers <- function(data, traits, multiplier = 1.5,
                            by = if ("location" %in% names(data)) "location" else NULL,
                            id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries",
                                        "accession", "row", "col"),
                            quiet = FALSE) {
  found <- find_outliers(data, traits, multiplier = multiplier, by = by, id_cols = id_cols)
  for (t in unique(found$trait)) data[[t]][found$.row[found$trait == t]] <- NA
  if (!isTRUE(quiet)) {
    n_values <- sum(vapply(traits, function(t) sum(!is.na(suppressWarnings(as.numeric(data[[t]])))), integer(1))) +
      nrow(found)
    if (!nrow(found)) {
      cat("No outliers (IQR x ", multiplier, ") among ", n_values, " values", if (length(by)) paste0(" (within each ", paste(by, collapse = ", "), ")"), ".\n", sep = "")
    } else {
      cat("Removed ", nrow(found), " of ", n_values, " values as outliers (outside Q1 - ", multiplier,
          " IQR and Q3 + ", multiplier, " IQR", if (length(by)) paste0(", within each ", paste(by, collapse = ", ")),
          "):\n", sep = "")
      show <- as.data.frame(found[, setdiff(names(found), ".row"), drop = FALSE])
      print(show, row.names = FALSE, digits = 4)
    }
  }
  attr(data, "outliers") <- found[, setdiff(names(found), ".row"), drop = FALSE]
  data
}

#' Plot a trait: distribution and field map
#'
#' Draws, for each location (or other group), a histogram of the trait and a
#' map of the field with the value of every plot, so gradients, odd plots
#' and outliers are easy to see. Outliers from [find_outliers()] can be
#' circled on the map.
#'
#' @param data A data frame with one row per plot.
#' @param trait Name of the numeric trait column.
#' @param by Optional grouping column; by default `"location"` when the data
#'   have it.
#' @param type `"both"` (histogram and field map), `"histogram"` or `"field"`.
#' @param row_col,col_col Field row and column columns (for the map).
#' @param outliers Optional table from [find_outliers()] to circle on the map.
#' @param breaks Histogram breaks, as in [graphics::hist()].
#'
#' @return Invisibly, a tibble with a summary per group: plots, values,
#'   missing, mean, sd, min, max and CV (%).
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' d$yield <- 50 + 0.8 * d$row + rnorm(nrow(d), sd = 2)
#' plot_trait(d, "yield")
#'
#' @export
plot_trait <- function(data, trait, by = if ("location" %in% names(data)) "location" else NULL,
                       type = c("both", "histogram", "field"), row_col = "row", col_col = "col",
                       outliers = NULL, breaks = "Sturges") {
  type <- match.arg(type)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  absent <- setdiff(c(trait, by), names(data))
  if (length(absent)) stop("Column(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  has_map <- all(c(row_col, col_col) %in% names(data))
  if (type == "field" && !has_map) stop("The field map needs `", row_col, "` and `", col_col, "` columns.", call. = FALSE)
  if (type == "both" && !has_map) type <- "histogram"
  x <- suppressWarnings(as.numeric(data[[trait]]))
  group <- if (length(by)) factor(as.character(data[[by]]), levels = unique(as.character(data[[by]]))) else
    factor(rep(trait, nrow(data)))
  groups <- levels(group)
  panels <- length(groups) * if (type == "both") 2 else 1
  ncol <- if (type == "both") 2 else min(length(groups), 3)
  old <- graphics::par(mfrow = c(ceiling(panels / ncol), ncol), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  zlim <- range(x, na.rm = TRUE)
  for (g in groups) {
    i <- which(group == g)
    title <- if (length(by)) paste(trait, "-", g) else trait
    if (type %in% c("both", "histogram")) {
      graphics::hist(x[i], breaks = breaks, col = "steelblue", border = "white", main = title,
                     xlab = trait, ylab = "Plots")
      graphics::abline(v = mean(x[i], na.rm = TRUE), lty = 2)
    }
    if (type %in% c("both", "field")) {
      mark <- NULL
      if (!is.null(outliers) && nrow(outliers) && "row" %in% names(outliers)) {
        o <- outliers[outliers$trait == trait, , drop = FALSE]
        if (length(by) && by %in% names(o)) o <- o[as.character(o[[by]]) == g, , drop = FALSE]
        mark <- o
      }
      bb_field_image(data[[row_col]][i], data[[col_col]][i], x[i], main = paste(title, "(field)"),
                     zlim = zlim, mark = mark)
    }
  }
  summary <- dplyr::bind_rows(lapply(groups, function(g) {
    v <- x[group == g]
    tibble::tibble(group = g, plots = length(v), values = sum(!is.na(v)), missing = sum(is.na(v)),
                   mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
                   min = suppressWarnings(min(v, na.rm = TRUE)), max = suppressWarnings(max(v, na.rm = TRUE)),
                   cv_percent = 100 * stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))
  }))
  invisible(summary)
}

# Field map: one cell per plot, row 1 at the top; NA plots grey; optional
# circles on marked plots (columns row and col).
bb_field_image <- function(row, col, value, main, zlim = range(value, na.rm = TRUE), mark = NULL) {
  row <- as.integer(row)
  col <- as.integer(col)
  nr <- max(row, na.rm = TRUE)
  nc <- max(col, na.rm = TRUE)
  m <- matrix(NA_real_, nr, nc)
  m[cbind(row, col)] <- value
  pal <- grDevices::hcl.colors(50, "viridis")
  graphics::plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5), xlab = "Column", ylab = "Row",
                 main = main, xaxs = "i", yaxs = "i")
  graphics::rect(0.5, nr + 0.5, nc + 0.5, 0.5, col = "grey85", border = NA)
  if (all(is.finite(zlim)) && diff(zlim) > 0) {
    graphics::image(seq_len(nc), seq_len(nr), t(m), col = pal, zlim = zlim, add = TRUE)
  } else {
    graphics::image(seq_len(nc), seq_len(nr), t(m), col = pal[25], add = TRUE)
  }
  if (!is.null(mark) && nrow(mark)) {
    graphics::points(mark$col, mark$row, pch = 1, cex = 2.2, lwd = 2, col = "red")
  }
  if (all(is.finite(zlim))) {
    graphics::mtext(sprintf("low %.3g (dark) - high %.3g (yellow)", zlim[1], zlim[2]), side = 3, line = 0.2,
                    cex = 0.7)
  }
}

#' Compare raw and first-stage adjusted values
#'
#' Puts each plot's raw value next to its first-stage value (the genotype's
#' adjusted mean plus the plot's residual, as in [first_stage_values()]),
#' summarises how much the first stage changed (per location: correlation of
#' raw and adjusted plot values, the spread before and after, and how much
#' the ranking of the genotypes changed) and draws, per location, the raw
#' and the adjusted field maps and a raw-versus-adjusted plot.
#'
#' @param fit An analysis (e.g. [analyze_augmented_row_column_design()] with
#'   `spatial = TRUE`) or a two-stage analysis ([analyze_two_stage()],
#'   [analyze_augmented_prep_two_stage()]), whose stage-1 fits are compared.
#' @param plot Logical. Draw the plots.
#' @param locations For two-stage analyses: the locations to draw (default
#'   all).
#'
#' @return Invisibly, a list with `plots` (one row per plot: location, plot,
#'   genotype, row, col, `raw`, `adjusted`, `change`) and `summary` (one row
#'   per location). The summary is also printed.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 50 + 0.8 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' cmp <- compare_first_stage(fit)
#'
#' @export
compare_first_stage <- function(fit, plot = TRUE, locations = NULL) {
  fits <- if (!is.null(fit$stage1_fits) && is.null(fit$plot_results)) fit$stage1_fits else list(fit)
  if (is.null(names(fits)) || length(fits) == 1 && identical(fits[[1]], fit)) names(fits) <- "trial"
  per_plot <- dplyr::bind_rows(lapply(names(fits), function(loc) bb_compare_one(fits[[loc]], loc)))
  summary <- dplyr::bind_rows(lapply(unique(per_plot$location), function(loc) {
    p <- per_plot[per_plot$location == loc & !is.na(per_plot$raw) & !is.na(per_plot$adjusted), , drop = FALSE]
    raw_g <- tapply(p$raw, p$genotype, mean)
    adj_g <- tapply(p$adjusted, p$genotype, mean)
    top <- max(1L, ceiling(0.1 * length(raw_g)))
    tibble::tibble(
      location = loc, plots = nrow(p),
      cor_raw_adjusted = stats::cor(p$raw, p$adjusted),
      sd_raw = stats::sd(p$raw), sd_adjusted = stats::sd(p$adjusted),
      genotype_rank_cor = stats::cor(raw_g, adj_g[names(raw_g)], method = "spearman"),
      top10pct_kept = mean(names(sort(-raw_g))[seq_len(top)] %in% names(sort(-adj_g))[seq_len(top)])
    )
  }))
  cat("Raw versus first-stage adjusted values:\n")
  print(as.data.frame(summary), row.names = FALSE, digits = 3)
  cat("(genotype_rank_cor: Spearman correlation of genotype means; top10pct_kept: share of the raw top 10%",
      "still in the adjusted top 10%)\n")
  if (isTRUE(plot)) {
    draw <- if (is.null(locations)) unique(per_plot$location) else intersect(locations, unique(per_plot$location))
    has_map <- !all(is.na(per_plot$row)) && !all(is.na(per_plot$col))
    old <- graphics::par(mfrow = c(length(draw), if (has_map) 3 else 1), mar = c(4, 4, 3, 1))
    on.exit(graphics::par(old), add = TRUE)
    for (loc in draw) {
      p <- per_plot[per_plot$location == loc, , drop = FALSE]
      label <- if (loc == "trial") "" else paste0(" - ", loc)
      if (has_map) {
        zlim <- range(c(p$raw, p$adjusted), na.rm = TRUE)
        bb_field_image(p$row, p$col, p$raw, main = paste0("Raw", label), zlim = zlim)
        bb_field_image(p$row, p$col, p$adjusted, main = paste0("Adjusted (first stage)", label), zlim = zlim)
      }
      graphics::plot(p$raw, p$adjusted, pch = 19, col = grDevices::adjustcolor("steelblue", 0.6),
                     xlab = "Raw value", ylab = "Adjusted value", main = paste0("Raw vs adjusted", label))
      graphics::abline(0, 1, lty = 2)
    }
  }
  invisible(list(plots = per_plot, summary = summary))
}

bb_compare_one <- function(fit, loc) {
  values <- bb_first_stage_plots(fit, loc)
  if (inherits(fit, "SpATS")) {
    d <- fit$data
    return(tibble::tibble(location = loc, plots = as.integer(d$plots), genotype = as.character(d$geno),
                          row = as.numeric(d$row), col = as.numeric(d$col), raw = d$y,
                          adjusted = values$value, change = values$value - d$y))
  }
  pr <- fit$plot_results
  response <- fit$input$response
  gcol <- fit$input$genotype_col %||% "all_entries"
  rcol <- if ("row" %in% names(pr)) "row" else NA
  ccol <- if ("col" %in% names(pr)) "col" else NA
  raw <- if (!is.null(response) && response %in% names(pr)) suppressWarnings(as.numeric(pr[[response]])) else
    pr$fitted + pr$residual
  # A two-stage fit is grouped by its named stage-1 trials. Breedbase can
  # omit location labels, or several trials can share a physical location.
  location <- rep(loc, nrow(pr))
  location_col <- fit$input$location_col %||% "location"
  if (identical(loc, "trial") && location_col %in% names(pr)) {
    labels <- as.character(pr[[location_col]])
    present <- !is.na(labels) & nzchar(trimws(labels))
    location[present] <- labels[present]
  }
  tibble::tibble(
    location = location,
    plots = as.integer(pr$plots), genotype = as.character(pr[[gcol]]),
    row = if (is.na(rcol)) NA_real_ else as.numeric(pr[[rcol]]),
    col = if (is.na(ccol)) NA_real_ else as.numeric(pr[[ccol]]),
    raw = raw, adjusted = values$value, change = values$value - raw
  )
}
