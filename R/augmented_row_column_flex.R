# ============================================================================
# Flexible augmented row-column designs and mixed/spatial analyses
# Piepho, H.-P. and Williams, E. R. (2016). Augmented row-column designs for a
# small number of checks. Agronomy Journal 108, 2256-2262.
# ============================================================================

# Group sizes and checks per block. Rows (columns) are grouped into row
# (column) groups of `kb` (`sb`) adjacent rows (columns); blocks are the
# intersections and block (i, j) receives `pb[i, j]` check plots.
arc_group_spec <- function(rows_in_field, cols_in_field, rows_per_block, cols_per_block,
                           rows_per_rowgroup, cols_per_colgroup, checks_per_block,
                           n_controls) {
  group_sizes <- function(total, per_block, per_group, total_label, per_block_label, group_label) {
    if (!is.null(per_group)) {
      out <- suppressWarnings(as.integer(per_group))
      if (length(out) < 1 || anyNA(out) || any(out < 1)) {
        stop("`", group_label, "` must contain positive integers.", call. = FALSE)
      }
      if (!is.null(total) && sum(out) != arc_positive_integer(total, total_label)) {
        stop("`", group_label, "` must add up to `", total_label, "`.", call. = FALSE)
      }
      return(out)
    }
    if (is.null(total) || is.null(per_block)) {
      stop(
        "Give `", total_label, "` and `", per_block_label, "`, or `", group_label, "`.",
        call. = FALSE
      )
    }
    total <- arc_positive_integer(total, total_label)
    per_block <- arc_positive_integer(per_block, per_block_label)
    if (total %% per_block != 0) {
      stop("`", total_label, "` must be divisible by `", per_block_label, "`.", call. = FALSE)
    }
    rep(per_block, total %/% per_block)
  }

  kb <- group_sizes(rows_in_field, rows_per_block, rows_per_rowgroup,
                    "rows_in_field", "rows_per_block", "rows_per_rowgroup")
  sb <- group_sizes(cols_in_field, cols_per_block, cols_per_colgroup,
                    "cols_in_field", "cols_per_block", "cols_per_colgroup")
  gk <- length(kb)
  gs <- length(sb)
  n_controls <- arc_positive_integer(n_controls, "n_controls")

  pb <- if (is.null(checks_per_block)) {
    matrix(n_controls, gk, gs)
  } else if (is.matrix(checks_per_block)) {
    if (!identical(dim(checks_per_block), c(gk, gs))) {
      stop("`checks_per_block` must be a ", gk, " x ", gs, " matrix (row groups x column groups).", call. = FALSE)
    }
    checks_per_block
  } else if (length(checks_per_block) == 1) {
    matrix(checks_per_block, gk, gs)
  } else if (length(checks_per_block) == gk * gs) {
    matrix(checks_per_block, gk, gs, byrow = TRUE)
  } else {
    stop(
      "`checks_per_block` must be a single number, a ", gk, " x ", gs,
      " matrix, or a vector of length ", gk * gs, " (by row group).",
      call. = FALSE
    )
  }
  storage.mode(pb) <- "integer"
  if (anyNA(pb) || any(pb < 0)) {
    stop("`checks_per_block` must contain non-negative integers.", call. = FALSE)
  }
  too_many <- pb > outer(kb, sb)
  if (any(too_many)) {
    stop("`checks_per_block` is larger than the number of plots of some blocks.", call. = FALSE)
  }

  list(
    kb = kb, sb = sb, pb = pb, gk = gk, gs = gs,
    k = sum(kb), s = sum(sb), n_controls = n_controls,
    regular = length(unique(kb)) == 1 && length(unique(sb)) == 1 && all(pb == n_controls)
  )
}

# Capacity and feasibility of a group specification.
arc_spec_capacity <- function(spec) {
  vc <- spec$n_controls
  check_plots <- sum(spec$pb)
  total <- spec$k * spec$s
  error_df <- check_plots - vc - (spec$k - 1) - (spec$s - 1) - (spec$gk - 1) * (spec$gs - 1)
  # Piepho and Williams (2016), observation 4: the row-group (column-group)
  # sub-designs need at least 1 + (vc - 1) + (gs - 1) + (kb(i) - 1) check
  # plots per row group (and analogously per column group).
  row_slack <- rowSums(spec$pb) - (vc + spec$gs + spec$kb - 2L)
  col_slack <- colSums(spec$pb) - (vc + spec$gk + spec$sb - 2L)
  tibble::tibble(
    rows_in_field = spec$k,
    cols_in_field = spec$s,
    rows_per_block = if (length(unique(spec$kb)) == 1) spec$kb[1] else NA_integer_,
    cols_per_block = if (length(unique(spec$sb)) == 1) spec$sb[1] else NA_integer_,
    n_block_rows = spec$gk,
    n_block_cols = spec$gs,
    n_blocks = spec$gk * spec$gs,
    n_controls = vc,
    total_plots = total,
    check_plots = check_plots,
    entry_plots = total - check_plots,
    error_df = error_df,
    check_fraction = check_plots / total,
    min_rowgroup_slack = min(row_slack),
    min_colgroup_slack = min(col_slack),
    every_row_has_check = all(rowSums(spec$pb) >= spec$kb),
    every_col_has_check = all(colSums(spec$pb) >= spec$sb)
  )
}

# ---- Flexible check allocation ("optimal" method) ---------------------------

arc_flex_layout <- function(spec) {
  row_group <- rep(seq_len(spec$gk), spec$kb)
  col_group <- rep(seq_len(spec$gs), spec$sb)
  list(row_group = row_group, col_group = col_group)
}

# Model matrix of the check plots for
# intercept + row + column + rowgroup x colgroup + check, and its information
# (rank and log-determinant of X'X over the estimable part).
arc_flex_info <- function(row, col, label, spec, layout) {
  n <- length(row)
  k <- spec$k
  s <- spec$s
  gk <- spec$gk
  gs <- spec$gs
  vc <- spec$n_controls
  p <- 1L + (k - 1L) + (s - 1L) + (gk - 1L) * (gs - 1L) + (vc - 1L)
  x <- matrix(0, n, p)
  x[, 1] <- 1
  idx <- seq_len(n)
  off <- 1L
  r_ok <- row > 1L
  x[cbind(idx[r_ok], off + row[r_ok] - 1L)] <- 1
  off <- off + k - 1L
  c_ok <- col > 1L
  x[cbind(idx[c_ok], off + col[c_ok] - 1L)] <- 1
  off <- off + s - 1L
  rg <- layout$row_group[row]
  cg <- layout$col_group[col]
  b_ok <- rg > 1L & cg > 1L
  if (gk > 1 && gs > 1) {
    x[cbind(idx[b_ok], off + (rg[b_ok] - 2L) * (gs - 1L) + cg[b_ok] - 1L)] <- 1
  }
  off <- off + (gk - 1L) * (gs - 1L)
  l_ok <- label > 1L
  x[cbind(idx[l_ok], off + label[l_ok] - 1L)] <- 1

  q <- qr(x)
  d <- abs(diag(q$qr))[seq_len(q$rank)]
  list(rank = q$rank, p = p, logdet = 2 * sum(log(d)))
}

# Number of pairs of check plots that are neighbours (same column in
# consecutive rows, or same row in consecutive columns), split into pairs
# inside one block and pairs across a block boundary.
arc_adjacent_checks <- function(row, col, k, s, row_group, col_group) {
  occupied <- matrix(FALSE, k, s)
  occupied[cbind(row, col)] <- TRUE
  vert <- occupied[-k, , drop = FALSE] & occupied[-1, , drop = FALSE]
  horiz <- occupied[, -s, drop = FALSE] & occupied[, -1, drop = FALSE]
  same_rg <- row_group[-k] == row_group[-1]
  same_cg <- col_group[-s] == col_group[-1]
  inside <- sum(vert[same_rg, , drop = FALSE]) + sum(horiz[, same_cg, drop = FALSE])
  c(inside = inside, across = sum(vert) + sum(horiz) - inside)
}

# Spread penalties: neighbouring check plots (`inside` / `across` blocks),
# uneven numbers of checks per row and column (`total`), and the same check
# more than once in a field row or column (`repeats`).
arc_flex_spread <- function(row, col, label, spec, layout = arc_flex_layout(spec)) {
  occupied <- matrix(FALSE, spec$k, spec$s)
  occupied[cbind(row, col)] <- TRUE
  adj <- arc_adjacent_checks(row, col, spec$k, spec$s, layout$row_group, layout$col_group)
  same_row <- sum(pmax(table(factor(label), factor(row, levels = seq_len(spec$k))) - 1, 0))
  same_col <- sum(pmax(table(factor(label), factor(col, levels = seq_len(spec$s))) - 1, 0))
  row_counts <- tabulate(row, spec$k)
  col_counts <- tabulate(col, spec$s)

  # Even numbers of check plots: within row group i every row gets
  # floor(T_i / kb_i) or ceiling(T_i / kb_i) check plots (T_i = check plots
  # of the row group), and likewise for columns within column groups.
  out_of_range <- function(counts, group, totals, sizes) {
    lo <- floor(totals / sizes)[group]
    hi <- ceiling(totals / sizes)[group]
    sum(pmax(lo - counts, 0) + pmax(counts - hi, 0))
  }
  unbalanced <- out_of_range(row_counts, layout$row_group, rowSums(spec$pb), spec$kb) +
    out_of_range(col_counts, layout$col_group, colSums(spec$pb), spec$sb)

  # Coverage: a row or column with n check plots holds every check at least
  # floor(n / vc) times (a column with as many check plots as checks holds
  # every check).
  vc <- spec$n_controls
  missing_checks <- function(line, n_lines) {
    counts <- matrix(tabulate((line - 1L) * vc + label, nbins = n_lines * vc), n_lines, vc, byrow = TRUE)
    need <- floor(tabulate(line, n_lines) / vc)
    sum(pmax(need - counts, 0))
  }
  uncovered <- missing_checks(col, spec$s) + missing_checks(row, spec$k)

  list(
    total = sum((row_counts - mean(row_counts))^2) + sum((col_counts - mean(col_counts))^2),
    repeats = same_row + same_col,
    adjacent_inside = adj[["inside"]],
    adjacent_across = adj[["across"]],
    unbalanced = unbalanced,
    uncovered = uncovered
  )
}

arc_flex_score <- function(state, spec, layout) {
  info <- arc_flex_info(state$row, state$col, state$label, spec, layout)
  spread <- arc_flex_spread(state$row, state$col, state$label, spec, layout)
  # Estimability first; then no neighbouring checks inside a block, even
  # numbers of check plots per row and column, and every check in each row or
  # column that has room for it; then no neighbours across blocks; then the
  # D-criterion; then fine balance and repeated checks.
  list(
    value = info$rank * 1e6 - 1e4 * spread$adjacent_inside - 1000 * spread$unbalanced -
      1000 * spread$uncovered - 100 * spread$adjacent_across +
      info$logdet - 1e-3 * spread$total - 0.05 * spread$repeats,
    rank = info$rank, p = info$p, logdet = info$logdet, spread = spread$total, repeats = spread$repeats,
    adjacent_inside = spread$adjacent_inside, adjacent_across = spread$adjacent_across,
    unbalanced = spread$unbalanced, uncovered = spread$uncovered
  )
}

# Random initial allocation. Step 2: check labels for each block, balanced
# over row groups, column groups and blocks. Steps 3-4: cells within blocks,
# spreading checks evenly over rows and columns.
arc_flex_initial <- function(spec, layout) {
  vc <- spec$n_controls
  row_start <- c(0L, cumsum(spec$kb))
  col_start <- c(0L, cumsum(spec$sb))
  total <- integer(vc)
  by_rg <- matrix(0L, spec$gk, vc)
  by_cg <- matrix(0L, spec$gs, vc)
  row_load <- integer(spec$k)
  col_load <- integer(spec$s)
  occupied <- matrix(FALSE, spec$k, spec$s)
  label_row <- matrix(0L, vc, spec$k)
  label_col <- matrix(0L, vc, spec$s)

  blocks <- expand.grid(i = seq_len(spec$gk), j = seq_len(spec$gs))
  blocks <- blocks[sample.int(nrow(blocks)), , drop = FALSE]
  out_row <- out_col <- out_label <- integer(0)

  for (b in seq_len(nrow(blocks))) {
    i <- blocks$i[b]
    j <- blocks$j[b]
    n_here <- spec$pb[i, j]
    if (n_here == 0) next
    here <- integer(vc)
    labels <- integer(n_here)
    for (t in seq_len(n_here)) {
      cost <- 1000 * here + 10 * by_rg[i, ] + 10 * by_cg[j, ] + total + stats::runif(vc)
      lab <- which.min(cost)
      labels[t] <- lab
      here[lab] <- here[lab] + 1L
      total[lab] <- total[lab] + 1L
      by_rg[i, lab] <- by_rg[i, lab] + 1L
      by_cg[j, lab] <- by_cg[j, lab] + 1L
    }
    rows <- (row_start[i] + 1L):row_start[i + 1L]
    cols <- (col_start[j] + 1L):col_start[j + 1L]
    for (lab in labels[sample.int(n_here)]) {
      cells <- expand.grid(r = rows, c = cols)
      cells <- cells[!occupied[cbind(cells$r, cells$c)], , drop = FALSE]
      neighbours <- vapply(seq_len(nrow(cells)), function(t) {
        r <- cells$r[t]
        cc <- cells$c[t]
        sum(occupied[cbind(pmin(pmax(r + c(-1L, 1L, 0L, 0L), 1L), spec$k),
                           pmin(pmax(cc + c(0L, 0L, -1L, 1L), 1L), spec$s))])
      }, numeric(1))
      cost <- 4 * row_load[cells$r] + 4 * col_load[cells$c] +
        20 * label_row[lab, cells$r] + 20 * label_col[lab, cells$c] +
        100 * neighbours + stats::runif(nrow(cells))
      pick <- which.min(cost)
      r <- cells$r[pick]
      cc <- cells$c[pick]
      occupied[r, cc] <- TRUE
      row_load[r] <- row_load[r] + 1L
      col_load[cc] <- col_load[cc] + 1L
      label_row[lab, r] <- label_row[lab, r] + 1L
      label_col[lab, cc] <- label_col[lab, cc] + 1L
      out_row <- c(out_row, r)
      out_col <- c(out_col, cc)
      out_label <- c(out_label, lab)
    }
  }
  list(row = out_row, col = out_col, label = out_label)
}

# Local search on the D-criterion (rank first, then log det X'X, then
# spread): move a check plot to a free cell of its block, or swap the labels
# of two check plots.
arc_flex_improve <- function(state, spec, layout, max_iter = NULL, patience = 1500L) {
  n <- length(state$row)
  if (is.null(max_iter)) max_iter <- min(20000L, 60L * n)
  current <- arc_flex_score(state, spec, layout)
  block_of <- function(r, cc) (layout$row_group[r] - 1L) * spec$gs + layout$col_group[cc]
  occupied <- matrix(FALSE, spec$k, spec$s)
  occupied[cbind(state$row, state$col)] <- TRUE
  # Check allocation of Step 2 kept by every move: in each block every check
  # appears between floor(pb / vc) and ceiling(pb / vc) times, so a block with
  # at least as many check plots as checks holds every check (complete), and
  # an incomplete block holds each check at most once.
  max_copies <- as.vector(t(ceiling(spec$pb / spec$n_controls)))
  min_copies <- as.vector(t(floor(spec$pb / spec$n_controls)))
  stale <- 0L
  for (iter in seq_len(max_iter)) {
    if (stale >= patience) break
    trial <- state
    u <- stats::runif(1)
    if (u < 0.6) {
      t <- sample.int(n, 1L)
      rg <- layout$row_group[state$row[t]]
      cg <- layout$col_group[state$col[t]]
      rows <- which(layout$row_group == rg)
      cols <- which(layout$col_group == cg)
      free <- which(!occupied[rows, cols, drop = FALSE], arr.ind = TRUE)
      if (nrow(free) == 0) {
        stale <- stale + 1L
        next
      }
      f <- free[sample.int(nrow(free), 1L), ]
      trial$row[t] <- rows[f[1]]
      trial$col[t] <- cols[f[2]]
    } else {
      pair <- sample.int(n, 2L)
      if (state$label[pair[1]] == state$label[pair[2]]) {
        stale <- stale + 1L
        next
      }
      trial$label[pair] <- state$label[rev(pair)]
      # Keep the check allocation of Step 2 valid (see max_copies).
      b <- block_of(trial$row[pair], trial$col[pair])
      if (b[1] != b[2]) {
        all_blocks <- block_of(trial$row, trial$col)
        ok <- vapply(1:2, function(t) {
          counts <- tabulate(trial$label[all_blocks == b[t]], nbins = spec$n_controls)
          all(counts <= max_copies[b[t]]) && all(counts >= min_copies[b[t]])
        }, logical(1))
        if (!all(ok)) {
          stale <- stale + 1L
          next
        }
      }
    }
    score <- arc_flex_score(trial, spec, layout)
    if (score$value > current$value + 1e-9) {
      stale <- 0L
    } else {
      stale <- stale + 1L
    }
    if (score$value >= current$value - 1e-9) {
      if (u < 0.6) {
        occupied[state$row[t], state$col[t]] <- FALSE
        occupied[trial$row[t], trial$col[t]] <- TRUE
      }
      state <- trial
      current <- score
    }
  }
  list(state = state, score = current)
}

arc_flex_search <- function(spec, n_candidates, patience_restarts = 3L) {
  layout <- arc_flex_layout(spec)
  best <- NULL
  since <- 0L
  # Keep restarting (up to 20 times) while a placement rule is still broken,
  # even if fewer candidates were requested.
  max_restarts <- 20L
  n_restarts <- max(min(n_candidates, 50L), max_restarts)
  evaluated <- 0L
  for (i in seq_len(n_restarts)) {
    evaluated <- i
    improved <- arc_flex_improve(arc_flex_initial(spec, layout), spec, layout)
    if (is.null(best) || improved$score$value > best$score$value + 1e-9) {
      best <- improved
      since <- 0L
    } else {
      since <- since + 1L
    }
    # Stop once the design is estimable, meets every placement rule
    # (no neighbouring checks in a block, even rows and columns, every check
    # where there is room for it) and further restarts stopped helping.
    clean <- best$score$rank == best$score$p && best$score$adjacent_inside == 0 &&
      best$score$unbalanced == 0 && best$score$uncovered == 0
    if (clean && (since >= patience_restarts || i >= n_candidates)) break
    if (!clean && i >= n_candidates && i >= max_restarts) break
  }
  c(best, list(n_evaluated = evaluated))
}

arc_flex_build_design <- function(state, spec, treatments, controls, plot_type) {
  layout <- arc_flex_layout(spec)
  field <- expand.grid(row = seq_len(spec$k), col = seq_len(spec$s), KEEP.OUT.ATTRS = FALSE)
  field$rowgroup <- layout$row_group[field$row]
  field$colgroup <- layout$col_group[field$col]
  field$block <- (field$rowgroup - 1L) * spec$gs + field$colgroup
  field$trt <- NA_character_
  field$type <- "entry"
  check_idx <- (state$col - 1L) * spec$k + state$row
  control_labels <- sample(controls)
  field$trt[check_idx] <- control_labels[state$label]
  field$type[check_idx] <- "check"
  entry_idx <- which(field$type == "entry")
  field$trt[entry_idx] <- treatments[sample.int(length(treatments))]
  field <- arc_make_plot_order(field, plot_type = plot_type)
  field$rep <- ifelse(field$type == "check", field$block, 1L)
  field$is_control <- ifelse(field$type == "check", 1L, 0L)
  field
}

# Neighbouring check plots of a finished design (any method).
arc_design_adjacent <- function(design, spec) {
  checks <- design[design$type == "check", , drop = FALSE]
  layout <- arc_flex_layout(spec)
  arc_adjacent_checks(checks$row, checks$col, spec$k, spec$s, layout$row_group, layout$col_group)
}

# D-criterion of a finished design (any method), as exp(logdet / p).
arc_design_d_criterion <- function(design, spec) {
  checks <- design[design$type == "check", , drop = FALSE]
  label <- as.integer(factor(checks$all_entries))
  info <- arc_flex_info(checks$row, checks$col, label, spec, arc_flex_layout(spec))
  if (info$rank < info$p) {
    return(0)
  }
  exp(info$logdet / info$p)
}

# ---- Mixed-model and spatial analyses ---------------------------------------

# Predictions K = [Lx, Lz] of fixed and random effects from an lmer fit with
# their exact prediction-error covariance sigma^2 K C^-1 K' from Henderson's
# mixed model equations (spherical parameterization of lme4).
arc_lmer_predict <- function(model, Lx, Lz = NULL) {
  X <- Matrix::Matrix(lme4::getME(model, "X"), sparse = TRUE)
  Z <- lme4::getME(model, "Z")
  Lambda <- lme4::getME(model, "Lambda")
  u <- lme4::getME(model, "u")
  sigma2 <- stats::sigma(model)^2
  ZL <- Z %*% Lambda
  q <- ncol(ZL)
  C <- rbind(
    cbind(Matrix::crossprod(X), Matrix::crossprod(X, ZL)),
    cbind(Matrix::crossprod(ZL, X), Matrix::crossprod(ZL) + Matrix::Diagonal(q))
  )
  if (is.null(Lz)) Lz <- matrix(0, nrow(Lx), q)
  K <- cbind(Lx, as.matrix(Lz %*% Lambda))
  estimate <- drop(Lx %*% lme4::fixef(model) + Lz %*% (Lambda %*% u))
  cov <- sigma2 * as.matrix(K %*% Matrix::solve(C, t(K)))
  list(estimate = as.numeric(estimate), cov = cov)
}

arc_variance_components <- function(model, labels) {
  vc <- as.data.frame(lme4::VarCorr(model))
  tibble::tibble(
    component = ifelse(vc$grp %in% names(labels), labels[vc$grp], vc$grp),
    variance = vc$vcov
  )
}

# Mixed model for an augmented row-column trial: random blocking factors
# (rows within row groups, columns within column groups, and blocks) for
# recovery of inter-block information, and/or random entries with fixed
# checks (Piepho, Williams and Fleck, 2006).
arc_rc_mixed <- function(data, random_blocks, random_entries, error_df) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Random blocks or random entries need the 'lme4' package.", call. = FALSE)
  }
  data <- droplevels(as.data.frame(data))
  data$block_nested <- interaction(data$rowgroup, data$colgroup, drop = TRUE)
  is_check <- data$is_check
  data$check <- factor(ifelse(is_check, as.character(data$genotype), ".entries"))
  data$entry <- factor(ifelse(is_check, ".checks", as.character(data$genotype)))
  data$conv <- as.numeric(!is_check)

  fixed_terms <- c(
    if (nlevels(data$rowgroup) > 1) "rowgroup",
    if (nlevels(data$colgroup) > 1) "colgroup",
    if (!random_blocks) c("row_nested", "col_nested", "block_nested"),
    if (random_entries) "check" else "genotype"
  )
  random_terms <- c(
    if (random_blocks) c("(1 | row_nested)", "(1 | col_nested)", if (nlevels(data$block_nested) > 1) "(1 | block_nested)"),
    if (random_entries) "(0 + conv | entry)"
  )
  fit0 <- arc_lm_drop_aliased_terms(
    stats::as.formula(paste("response ~", paste(fixed_terms, collapse = " + "))),
    data
  )
  x_full <- stats::model.matrix(fit0$terms, fit0$model)
  assign <- attr(x_full, "assign")
  kept <- match(colnames(fit0$x), colnames(x_full))

  mixed_data <- data
  mixed_data$X <- fit0$x
  model <- suppressMessages(lme4::lmer(
    stats::as.formula(paste("response ~ 0 + X +", paste(random_terms, collapse = " + "))),
    data = mixed_data,
    REML = TRUE,
    control = lme4::lmerControl(check.rankX = "ignore", calc.derivs = FALSE)
  ))

  factor_term <- if (random_entries) "check" else "genotype"
  term_index <- match(factor_term, attr(fit0$terms, "term.labels"))
  term_cols <- which(assign == term_index)
  base <- colMeans(x_full)
  base[term_cols] <- 0
  levels_f <- levels(data[[factor_term]])
  Lx_level <- matrix(base, length(levels_f), ncol(x_full), byrow = TRUE)
  hit <- match(colnames(x_full)[term_cols], paste0(factor_term, levels_f))
  Lx_level[cbind(hit, term_cols)] <- 1
  rownames(Lx_level) <- levels_f

  qr_full <- qr(x_full)
  estimable_level <- arc_linear_estimates(qr_full, qr.coef(qr_full, fit0$y), Lx_level)$estimable
  names(estimable_level) <- levels_f

  q <- ncol(lme4::getME(model, "Z"))
  if (!random_entries) {
    genotypes <- levels_f
    Lx <- Lx_level[, kept, drop = FALSE]
    Lz <- NULL
    estimable <- unname(estimable_level)
  } else {
    checks <- setdiff(levels_f, ".entries")
    entries <- setdiff(levels(data$entry), ".checks")
    genotypes <- c(checks, entries)
    Lx <- rbind(
      Lx_level[checks, kept, drop = FALSE],
      Lx_level[rep(".entries", length(entries)), kept, drop = FALSE]
    )
    # Columns of Z for the random entry effects.
    zt <- lme4::getME(model, "Ztlist")
    offsets <- c(0L, cumsum(vapply(zt, nrow, integer(1))))
    term <- grep("^entry\\.", names(zt))
    entry_levels <- levels(lme4::getME(model, "flist")$entry)
    Lz <- matrix(0, length(genotypes), q)
    Lz[cbind(length(checks) + seq_along(entries), offsets[term] + match(entries, entry_levels))] <- 1
    estimable <- unname(c(estimable_level[checks], rep(estimable_level[".entries"], length(entries))))
  }

  pred <- arc_lmer_predict(model, Lx, Lz)
  estimate <- pred$estimate
  std_error <- sqrt(pmax(diag(pred$cov), 0))
  estimate[!estimable] <- NA_real_
  std_error[!estimable] <- NA_real_
  ok <- which(estimable)
  vcov <- pred$cov[ok, ok, drop = FALSE]
  dimnames(vcov) <- list(genotypes[ok], genotypes[ok])
  mean_sed <- prep_mean_sed(vcov)

  list(
    model = model,
    # Plot-level fitted values (fixed effects + predicted random effects).
    fitted = unname(stats::fitted(model)),
    adjusted_means = tibble::tibble(
      genotype = genotypes,
      adjusted_mean = estimate,
      std_error = std_error,
      estimable = estimable
    ),
    vcov = vcov,
    precision = list(
      mean_sed = mean_sed,
      lsd_5pct = if (!is.na(mean_sed) && error_df > 0) stats::qt(0.975, error_df) * mean_sed else NA_real_
    ),
    variance_components = arc_variance_components(model, c(
      row_nested = "rows within row groups",
      col_nested = "columns within column groups",
      block_nested = "blocks (row group x column group)",
      entry = "entries",
      Residual = "residual"
    ))
  )
}

# Spatial analysis with SpATS: genotypes fixed, a 2-D P-spline surface over
# rows and columns, row and column groups fixed, and rows within row groups,
# columns within column groups and blocks random.
arc_rc_spatial <- function(data, spatial_nseg, error_df) {
  if (!requireNamespace("SpATS", quietly = TRUE)) {
    stop("`spatial = TRUE` needs the 'SpATS' package.", call. = FALSE)
  }
  x <- droplevels(as.data.frame(data))
  x$rown <- as.numeric(as.character(x$row))
  x$coln <- as.numeric(as.character(x$col))
  if (anyNA(x$rown) || anyNA(x$coln)) {
    stop("The spatial model needs numeric field row and column coordinates.", call. = FALSE)
  }
  x$block_nested <- interaction(x$rowgroup, x$colgroup, drop = TRUE)
  x$geno <- x$genotype
  nseg <- if (is.null(spatial_nseg)) {
    pmax(pmin(ceiling(c(length(unique(x$coln)), length(unique(x$rown))) / 2), 40L), 3L)
  } else {
    as.integer(rep_len(spatial_nseg, 2))
  }
  spatial <- stats::as.formula(sprintf("~ PSANOVA(coln, rown, nseg = c(%d, %d))", nseg[1], nseg[2]))
  environment(spatial) <- asNamespace("SpATS")
  fixed_terms <- c(if (nlevels(x$rowgroup) > 1) "rowgroup", if (nlevels(x$colgroup) > 1) "colgroup")
  random_terms <- c("row_nested", "col_nested", if (nlevels(x$block_nested) > 1) "block_nested")
  fit <- SpATS::SpATS(
    response = "response",
    genotype = "geno",
    genotype.as.random = FALSE,
    spatial = spatial,
    fixed = if (length(fixed_terms)) stats::as.formula(paste("~", paste(fixed_terms, collapse = " + "))) else NULL,
    random = stats::as.formula(paste("~", paste(random_terms, collapse = " + "))),
    data = x,
    control = list(monitoring = 0)
  )
  pred <- stats::predict(fit, which = "geno", predFixed = "marginal", return.vcov.matrix = TRUE)
  vc <- as.matrix(attr(pred, "vcov"))
  genos <- as.character(pred$geno)
  dimnames(vc) <- list(genos, genos)
  mean_sed <- prep_mean_sed(vc)
  list(
    model = fit,
    # Plot-level fitted values (genotype + spatial surface + random effects).
    fitted = unname(fit$fitted),
    adjusted_means = tibble::tibble(
      genotype = genos,
      adjusted_mean = pred$predicted.values,
      std_error = pred$standard.errors,
      estimable = TRUE
    ),
    vcov = vc,
    precision = list(
      mean_sed = mean_sed,
      lsd_5pct = if (!is.na(mean_sed) && error_df > 0) stats::qt(0.975, error_df) * mean_sed else NA_real_
    ),
    variance_components = tibble::tibble(component = "residual", variance = unname(fit$psi[1])),
    spatial_effective_dim = sum(fit$eff.dim[grepl("^f\\(", names(fit$eff.dim))])
  )
}

#' List the analyzable layouts of an augmented row-column field
#'
#' Tries every combination of block size and number of checks for a field
#' and returns the layouts that can be analysed, so a layout can be chosen
#' without writing a loop over [augmented_row_column_capacity()]. Each block
#' gets one plot of every check. Suggested blocks have at least as many
#' rows and at least as many columns as there are checks, leaving room for
#' each check to occupy a distinct row and column.
#'
#' There is a trade-off: more check plots mean fewer entries but more error
#' degrees of freedom, i.e. more reliable F-tests and spatial adjustment. Aim
#' for at least about 10 error df when the field allows it (Clarke and
#' Stefanova, 2011, recommend about 50 for large trials).
#'
#' @param rows_in_field,cols_in_field Field size.
#' @param rows_per_block,cols_per_block Block sizes to try. By default every
#'   divisor of the field size from 2 up to the field size.
#' @param n_controls Numbers of checks to try.
#' @param min_error_df Smallest error df to keep.
#' @param n_entries Optional number of entries to plant. Only layouts with at
#'   least that many entry plots are kept, and `spare_plots` shows how many
#'   entry plots would be left over.
#'
#' @return A tibble with one row per analyzable layout, sorted by entry plots
#'   (most first) and then error df: block sizes, number of checks, blocks,
#'   entry and check plots, `check_fraction`, `error_df`, `method` (the
#'   construction that [augmented_row_column_design()] uses for it: `"latin"`
#'   or `"optimal"`) and, when `n_entries` is given, `spare_plots`.
#'
#' @seealso [augmented_row_column_capacity()], [augmented_row_column_design()]
#'
#' @examples
#' # Suggested analyzable layouts of a 12 x 12 field with 2 to 5 checks:
#' augmented_row_column_layouts(12, 12, n_controls = 2:5)
#'
#' # Layouts that fit 100 entries with at least 10 error df:
#' augmented_row_column_layouts(12, 12, n_controls = 2:6, min_error_df = 10, n_entries = 100)
#'
#' @export
augmented_row_column_layouts <- function(rows_in_field, cols_in_field,
                                         rows_per_block = NULL, cols_per_block = NULL,
                                         n_controls = 2:6, min_error_df = 1,
                                         n_entries = NULL) {
  rows_in_field <- arc_positive_integer(rows_in_field, "rows_in_field")
  cols_in_field <- arc_positive_integer(cols_in_field, "cols_in_field")
  divisors <- function(n) {
    d <- seq_len(n)
    d[n %% d == 0 & d >= 2]
  }
  if (is.null(rows_per_block)) rows_per_block <- divisors(rows_in_field)
  if (is.null(cols_per_block)) cols_per_block <- divisors(cols_in_field)
  grid <- expand.grid(rows_per_block = as.integer(rows_per_block),
                      cols_per_block = as.integer(cols_per_block),
                      n_controls = as.integer(n_controls))
  grid <- grid[rows_in_field %% grid$rows_per_block == 0 & cols_in_field %% grid$cols_per_block == 0 &
                 grid$rows_per_block >= grid$n_controls & grid$cols_per_block >= grid$n_controls &
                 grid$n_controls >= 1 & grid$n_controls <= grid$rows_per_block * grid$cols_per_block,
               , drop = FALSE]

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    cap <- augmented_row_column_capacity(rows_in_field, cols_in_field, g$rows_per_block,
                                         g$cols_per_block, g$n_controls)
    if (cap$entry_plots < 1 || cap$error_df < min_error_df ||
        !cap$every_row_has_check || !cap$every_col_has_check) {
      return(NULL)
    }
    latin <- cap$n_block_cols <= g$rows_per_block && cap$n_block_rows <= g$cols_per_block &&
      g$n_controls <= g$rows_per_block && g$n_controls <= g$cols_per_block
    cap$method <- if (latin) "latin" else "optimal"
    cap
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    message("No analyzable layout: try more checks, other block sizes, or a smaller `min_error_df`.")
    return(out)
  }
  out <- out[, c("rows_per_block", "cols_per_block", "n_controls", "entry_plots", "error_df",
                 "check_plots", "n_blocks", "check_fraction", "method")]
  if (!is.null(n_entries)) {
    n_entries <- arc_positive_integer(n_entries, "n_entries")
    out <- out[out$entry_plots >= n_entries, , drop = FALSE]
    out$spare_plots <- out$entry_plots - n_entries
  }
  out <- out[order(-out$entry_plots, -out$error_df), , drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
