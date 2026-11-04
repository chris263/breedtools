# ============================================================================
# Algebraic constructions of augmented p-rep designs
# Williams, Piepho and Whitaker (2011), Biometrical Journal 53, 19-27.
#
#   method = "alpha"  : partial alpha-arrays (Sections 2.1 and 2.2)
#   method = "alpha2" : augmented alpha(2)-arrays (Section 2.4)
#
# Notation (paper): v entries, l locations, m p-rep entries per block,
# s blocks per replicate, k final block size.
# ============================================================================

#' Feasible parameters for augmented p-rep constructions
#'
#' Lists the parameter sets that `augmented_prep_design()` can build for a
#' given number of entries and locations, in the spirit of Table 3 of
#' Williams, Piepho and Whitaker (2011). For each number `m` of p-rep entries
#' per block, it gives the number of blocks per replicate `s`, the final block
#' size(s), the number of entries duplicated at each location, and how many
#' "dummy" entries are dropped when the number of entries is not a multiple of
#' the array size (dropping them makes some blocks one plot smaller).
#'
#' @param n_treatments Number of entries (`v`).
#' @param n_locations Number of locations (`l`).
#' @param max_block_size Largest final block size to list.
#'
#' @return A tibble with one row per method and `m`: `method`, `m`, `s`,
#'   `block_size` (largest final block), `min_block_size`,
#'   `blocks_per_location`, `plots_per_location` (largest),
#'   `duplicated_per_location` (largest), `dummies`, `feasible`, and `note`.
#'
#' @references Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented
#'   p-rep designs. Biometrical Journal 53, 19-27. \doi{10.1002/bimj.201000102}
#'
#' @examples
#' augmented_prep_parameters(n_treatments = 120, n_locations = 6, max_block_size = 10)
#'
#' @export
augmented_prep_parameters <- function(n_treatments, n_locations, max_block_size = 20) {
  v <- arc_positive_integer(n_treatments, "n_treatments")
  l <- arc_positive_integer(n_locations, "n_locations")
  max_block_size <- arc_positive_integer(max_block_size, "max_block_size")

  rows <- list()
  for (m in seq_len(max(1L, 2L * max_block_size))) {
    for (method in c("alpha", "alpha2")) {
      par <- tryCatch(
        prep_algebraic_parameters(method, v, l, m, check = FALSE),
        error = function(e) NULL
      )
      if (is.null(par) || par$k_max > max_block_size) next
      rows[[length(rows) + 1L]] <- tibble::tibble(
        method = method,
        m = m,
        s = par$s,
        block_size = par$k_max,
        min_block_size = par$k_min,
        blocks_per_location = 2L * par$s,
        plots_per_location = v + par$dup_max,
        duplicated_per_location = par$dup_max,
        dummies = par$d,
        feasible = par$feasible,
        note = par$note
      )
    }
  }
  if (length(rows) == 0) {
    return(tibble::tibble(
      method = character(), m = integer(), s = integer(), block_size = integer(),
      min_block_size = integer(), blocks_per_location = integer(),
      plots_per_location = integer(), duplicated_per_location = integer(),
      dummies = integer(), feasible = logical(), note = character()
    ))
  }
  out <- dplyr::bind_rows(rows)
  out[order(out$method, out$m), , drop = FALSE]
}

# Parameters of an algebraic construction. With `check = TRUE` an infeasible
# set raises an error explaining why.
prep_algebraic_parameters <- function(method, v, l, m, check = TRUE) {
  fail <- function(msg) {
    if (check) stop(msg, call. = FALSE)
    NULL
  }
  note <- ""
  feasible <- TRUE
  if (method == "alpha") {
    if (m < 2L) return(fail("`n_prep_per_block` (m) must be at least 2 so each location is connected."))
    rows <- l * m
    s <- as.integer(ceiling(v / rows))
    if (s < 2L) return(fail("Too few entries for this `n_prep_per_block`: fewer than 2 blocks per replicate."))
    d <- rows * s - v
    k_max <- m + as.integer(ceiling((rows - m) / 2))
    k_min <- m + as.integer(floor((rows - m) / 2))
    if (d > 0) k_min <- k_min - as.integer(ceiling(d / s))
    if (rows > s) note <- "k' > s: start array is alpha(0,1,2); deletions are searched to reach (0,1)"
  } else if (method == "alpha2") {
    if (m < 2L || m %% 2L != 0L) return(fail("Method \"alpha2\" needs an even `n_prep_per_block` (m)."))
    rows <- l * m %/% 2L
    s <- as.integer(ceiling(v / (2L * rows)))
    if (s < 2L) return(fail("Too few entries for this `n_prep_per_block`: fewer than 2 blocks per replicate."))
    if (s %% 2L == 0L) {
      feasible <- FALSE
      note <- "s is even: the two-replicate design at each location would be disconnected"
      if (check) {
        stop(
          "Method \"alpha2\" needs an odd number of blocks per replicate, but these ",
          "parameters give s = ", s, ". Use method \"alpha\", or change the number of ",
          "entries or `n_prep_per_block`.",
          call. = FALSE
        )
      }
    }
    d <- 2L * rows * s - v
    k_max <- rows + m %/% 2L
    k_min <- k_max - as.integer(ceiling(d / (2L * s)))
  } else {
    stop("Unknown method.", call. = FALSE)
  }
  dup <- m * s
  list(
    method = method, v = v, l = l, m = m, rows = rows, s = s, d = d,
    k_max = k_max, k_min = max(1L, k_min), dup_max = dup,
    feasible = feasible, note = note
  )
}

# Largest m whose final block size fits `max_block_size`.
prep_choose_m <- function(method, v, l, max_block_size) {
  candidates <- if (method == "alpha2") seq(2L, 2L * max_block_size, by = 2L) else seq(2L, 2L * max_block_size)
  best <- NULL
  for (m in candidates) {
    par <- tryCatch(prep_algebraic_parameters(method, v, l, m, check = FALSE), error = function(e) NULL)
    if (is.null(par) || !par$feasible || par$k_max > max_block_size) next
    best <- m
  }
  if (is.null(best)) {
    stop(
      "No feasible `n_prep_per_block` for method \"", method, "\" with block_size <= ",
      max_block_size, ". See `augmented_prep_parameters()`.",
      call. = FALSE
    )
  }
  best
}

# Counts of entry pairs by concurrence from (key, multiplicity, pairs per key).
prep_concurrence_counts <- function(keys, pairs_per_key) {
  if (length(keys) == 0) {
    return(c(n_pairs_gt1 = 0, excess = 0, max = 0))
  }
  o <- order(keys)
  r <- rle(keys[o])
  w <- pairs_per_key[o][cumsum(r$lengths)]
  n <- r$lengths
  c(
    n_pairs_gt1 = sum(w[n >= 2L]),
    excess = sum((w * (n - 1L))[n >= 2L]),
    max = max(n)
  )
}

prep_algebraic_score <- function(stats, extra_components) {
  extra_components * 1e10 + stats[["n_pairs_gt1"]] * 1e7 + stats[["excess"]] * 1e5 + stats[["max"]] * 1000
}

# ---------------------------------------------------------------------------
# alpha: partial alpha-array (k' x 2l array mod s, with deletions)
# ---------------------------------------------------------------------------

prep_alpha_random_state <- function(par) {
  rows <- par$rows
  l <- par$l
  m <- par$m
  s <- par$s
  a <- matrix(sample.int(s, rows * 2L * l, replace = TRUE) - 1L, rows, 2L * l)
  a[1, ] <- 0L
  a[, 1] <- 0L
  # keep[i, c]: row i is used in column c. Boxed rows of location j are kept
  # in both of its columns; every other row is kept in exactly one of them.
  keep <- matrix(FALSE, rows, 2L * l)
  n_first <- as.integer(ceiling((rows - m) / 2))
  for (j in seq_len(l)) {
    boxed <- ((j - 1L) * m + 1L):(j * m)
    others <- setdiff(seq_len(rows), boxed)
    c1 <- 2L * j - 1L
    c2 <- 2L * j
    keep[boxed, c(c1, c2)] <- TRUE
    if (length(others) > 0) {
      first <- others[sample.int(length(others), n_first)]
      keep[first, c1] <- TRUE
      keep[setdiff(others, first), c2] <- TRUE
    }
  }
  list(a = a, keep = keep)
}

prep_alpha_state_score <- function(state, par) {
  a <- state$a
  keep <- state$keep
  rows <- par$rows
  s <- par$s
  pair_id <- matrix(seq_len(rows * rows), rows, rows)
  ut <- upper.tri(pair_id)
  keys <- lapply(seq_len(ncol(a)), function(c) {
    kept <- keep[, c]
    both <- ut & outer(kept, kept)
    diff <- outer(a[, c], a[, c], "-") %% s
    (pair_id[both] - 1) * s + diff[both]
  })
  keys <- unlist(keys)
  stats <- prep_concurrence_counts(keys, rep(s, length(keys)))
  comps <- prep_alpha_components(state, par)
  list(score = prep_algebraic_score(stats, sum(comps - 1L)), stats = stats, components = comps)
}

# Connected components of the two-replicate design at each location: the
# blocks are linked by b ~ b + delta_i (mod s), delta_i = a[i, c1] - a[i, c2],
# so the number of components is gcd(s, delta_i - delta_1, ...).
prep_alpha_components <- function(state, par) {
  vapply(seq_len(par$l), function(j) {
    boxed <- ((j - 1L) * par$m + 1L):(j * par$m)
    delta <- (state$a[boxed, 2L * j - 1L] - state$a[boxed, 2L * j]) %% par$s
    as.integer(Reduce(prep_gcd, (delta - delta[1]) %% par$s, par$s))
  }, integer(1))
}

prep_gcd <- function(a, b) {
  a <- abs(as.integer(a))
  b <- abs(as.integer(b))
  while (b > 0L) {
    t <- b
    b <- a %% b
    a <- t
  }
  a
}

prep_alpha_move <- function(state, par) {
  if (stats::runif(1) < 0.7 || par$rows == par$m) {
    if (par$rows < 2L || ncol(state$a) < 2L) return(state)
    i <- sample.int(par$rows - 1L, 1L) + 1L
    c <- sample.int(ncol(state$a) - 1L, 1L) + 1L
    state$a[i, c] <- sample.int(par$s, 1L) - 1L
  } else {
    j <- sample.int(par$l, 1L)
    boxed <- ((j - 1L) * par$m + 1L):(j * par$m)
    c1 <- 2L * j - 1L
    c2 <- 2L * j
    in1 <- setdiff(which(state$keep[, c1]), boxed)
    in2 <- setdiff(which(state$keep[, c2]), boxed)
    if (length(in1) == 0 || length(in2) == 0) return(state)
    r1 <- in1[sample.int(length(in1), 1L)]
    r2 <- in2[sample.int(length(in2), 1L)]
    state$keep[r1, c(c1, c2)] <- c(FALSE, TRUE)
    state$keep[r2, c(c1, c2)] <- c(TRUE, FALSE)
  }
  state
}

# Blocks of an alpha state: one data frame of plots (entry index, location,
# column-replicate, block within replicate) for the full array of v' entries.
prep_alpha_plots <- function(state, par) {
  s <- par$s
  out <- vector("list", ncol(state$a))
  for (c in seq_len(ncol(state$a))) {
    rows_kept <- which(state$keep[, c])
    b <- rep(seq_len(s) - 1L, each = length(rows_kept))
    i <- rep(rows_kept, times = s)
    x <- (state$a[i, c] + b) %% s
    out[[c]] <- data.frame(
      index = (i - 1L) * s + x + 1L,
      row = i,
      location = (c + 1L) %/% 2L,
      replicate = 2L - c %% 2L,
      block = b + 1L
    )
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# alpha2: augmented alpha(2)-array (k x l array of 2-ples mod (2, s))
# ---------------------------------------------------------------------------

prep_alpha2_random_state <- function(par) {
  rows <- par$rows
  l <- par$l
  s <- par$s
  au <- matrix(sample.int(2L, rows * l, replace = TRUE) - 1L, rows, l)
  aw <- matrix(sample.int(s, rows * l, replace = TRUE) - 1L, rows, l)
  au[1, ] <- 0L
  aw[1, ] <- 0L
  au[, 1] <- 0L
  aw[, 1] <- 0L
  # Extra 2-ple on the diagonal: rows of location j get a second 2-ple in
  # column j whose first element differs from the existing one.
  bw <- sample.int(s, rows, replace = TRUE) - 1L
  list(au = au, aw = aw, bw = bw)
}

prep_alpha2_rows_of <- function(par, j) {
  half <- par$m %/% 2L
  ((j - 1L) * half + 1L):(j * half)
}

# Virtual rows of column j: the main 2-ple of every row plus the extra 2-ple
# of the rows duplicated at location j.
prep_alpha2_virtual <- function(state, par, j) {
  extra <- prep_alpha2_rows_of(par, j)
  data.frame(
    row = c(seq_len(par$rows), extra),
    u = c(state$au[, j], 1L - state$au[extra, j]),
    w = c(state$aw[, j], state$bw[extra])
  )
}

prep_alpha2_state_score <- function(state, par) {
  s <- par$s
  rows <- par$rows
  keys <- list()
  weights <- list()
  for (j in seq_len(par$l)) {
    vr <- prep_alpha2_virtual(state, par, j)
    n_vr <- nrow(vr)
    idx <- which(upper.tri(matrix(0, n_vr, n_vr)), arr.ind = TRUE)
    r1 <- vr$row[idx[, 1]]
    r2 <- vr$row[idx[, 2]]
    du <- (vr$u[idx[, 1]] - vr$u[idx[, 2]]) %% 2L
    dw <- (vr$w[idx[, 1]] - vr$w[idx[, 2]]) %% s
    # Orient each pair so that the lower row comes first.
    flip <- r1 > r2
    dw[flip] <- (-dw[flip]) %% s
    lo <- pmin(r1, r2)
    hi <- pmax(r1, r2)
    same <- lo == hi
    # Within a row, delta and -delta give the same pairs of entries.
    code <- du * s + dw
    code_neg <- du * s + (-dw) %% s
    code[same] <- pmin(code[same], code_neg[same])
    keys[[j]] <- ((lo - 1) * rows + (hi - 1)) * (2 * s) + code
    weights[[j]] <- ifelse(same & dw == 0L, s, 2L * s)
  }
  stats <- prep_concurrence_counts(unlist(keys), unlist(weights))
  comps <- prep_alpha2_components(state, par)
  list(score = prep_algebraic_score(stats, sum(comps - 1L)), stats = stats, components = comps)
}

# Components of the two-replicate design at each location: blocks of
# Z2 x Zs linked by b ~ b + delta_i; the component of 0 is the subgroup
# generated by the delta_i.
prep_alpha2_components <- function(state, par) {
  s <- par$s
  vapply(seq_len(par$l), function(j) {
    extra <- prep_alpha2_rows_of(par, j)
    du <- rep(1L, length(extra))
    dw <- (state$aw[extra, j] - state$bw[extra]) %% s
    reach <- 0L
    repeat {
      u <- reach %/% s
      w <- reach %% s
      nxt <- unique(c(reach, unlist(lapply(seq_along(du), function(t) {
        ((u + du[t]) %% 2L) * s + (w + dw[t]) %% s
      }))))
      if (length(nxt) == length(reach)) break
      reach <- nxt
    }
    as.integer((2L * s) / length(reach))
  }, integer(1))
}

prep_alpha2_move <- function(state, par) {
  r <- stats::runif(1)
  if (r < 0.25) {
    state$bw[sample.int(par$rows, 1L)] <- sample.int(par$s, 1L) - 1L
  } else if (par$rows >= 2L && par$l >= 2L) {
    i <- sample.int(par$rows - 1L, 1L) + 1L
    j <- sample.int(par$l - 1L, 1L) + 1L
    if (r < 0.4) {
      state$au[i, j] <- 1L - state$au[i, j]
    } else {
      state$aw[i, j] <- sample.int(par$s, 1L) - 1L
    }
  } else {
    state$bw[sample.int(par$rows, 1L)] <- sample.int(par$s, 1L) - 1L
  }
  state
}

prep_alpha2_plots <- function(state, par) {
  s <- par$s
  out <- vector("list", par$l)
  for (j in seq_len(par$l)) {
    vr <- prep_alpha2_virtual(state, par, j)
    blocks <- expand.grid(w = seq_len(s) - 1L, u = 0:1)
    bi <- rep(seq_len(nrow(blocks)), each = nrow(vr))
    vi <- rep(seq_len(nrow(vr)), times = nrow(blocks))
    u <- (vr$u[vi] + blocks$u[bi]) %% 2L
    w <- (vr$w[vi] + blocks$w[bi]) %% s
    out[[j]] <- data.frame(
      index = (vr$row[vi] - 1L) * 2L * s + u * s + w + 1L,
      row = vr$row[vi],
      location = j,
      replicate = blocks$u[bi] + 1L,
      block = blocks$w[bi] + 1L
    )
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Shared driver
# ---------------------------------------------------------------------------

# Local search over arrays from random starts. Each start is improved by
# random single-element moves (ties accepted, so the search can drift across
# plateaus) until it reaches the best attainable score or stops improving.
# Returns up to `n_keep` distinct optimal states (or the best state found).
prep_algebraic_search <- function(method, par, n_starts, n_keep, max_moves = NULL, patience = 25L) {
  random_state <- if (method == "alpha") prep_alpha_random_state else prep_alpha2_random_state
  score_state <- if (method == "alpha") prep_alpha_state_score else prep_alpha2_state_score
  move <- if (method == "alpha") prep_alpha_move else prep_alpha2_move
  bound <- 1000
  if (is.null(max_moves)) max_moves <- max(2000L, 40L * par$rows * par$l)

  found <- list()
  best <- NULL
  n_evaluated <- 0L
  since_progress <- 0L
  for (start in seq_len(n_starts)) {
    n_evaluated <- start
    state <- random_state(par)
    current <- score_state(state, par)
    stale <- 0L
    moves <- 0L
    while (current$score > bound && moves < max_moves && stale < 600L) {
      moves <- moves + 1L
      trial <- move(state, par)
      trial_score <- score_state(trial, par)
      if (trial_score$score <= current$score) {
        if (trial_score$score < current$score) stale <- 0L else stale <- stale + 1L
        state <- trial
        current <- trial_score
      } else {
        stale <- stale + 1L
      }
    }
    since_progress <- since_progress + 1L
    if (is.null(best) || current$score < best$score) {
      best <- list(state = state, score = current$score)
      since_progress <- 0L
    }
    if (current$score <= bound) {
      found[[length(found) + 1L]] <- state
      since_progress <- 0L
      if (length(found) >= n_keep) break
    }
    # Stop when further random starts no longer help.
    if (since_progress >= patience) break
  }
  optimal <- length(found) > 0
  if (!optimal) found <- list(best$state)
  list(states = found, n_evaluated = n_evaluated, optimal = optimal)
}

# Turn an array state into a candidate for `prep_build_candidate_design()`:
# drop dummy entries, randomize entry labels, block order within replicates
# and plot order within blocks.
prep_algebraic_candidate <- function(method, state, par) {
  plots <- if (method == "alpha") prep_alpha_plots(state, par) else prep_alpha2_plots(state, par)
  v_full <- par$rows * (if (method == "alpha") par$s else 2L * par$s)

  # Dummies: spread round-robin over the locations (each array row is the
  # p-rep row of one location) and, within a row, choose the element that
  # keeps block sizes even without disconnecting that location.
  dummies <- integer(0)
  if (par$d > 0) {
    rows_per_loc <- if (method == "alpha") par$m else par$m %/% 2L
    per_row <- if (method == "alpha") par$s else 2L * par$s
    plots$gblock <- paste(plots$location, plots$replicate, plots$block)
    deficit <- stats::setNames(rep(0, length(unique(plots$gblock))), unique(plots$gblock))
    row_loc <- (seq_len(par$rows) - 1L) %/% rows_per_loc + 1L
    row_order <- unlist(lapply(seq_len(rows_per_loc), function(r) {
      rev((seq_len(par$l) - 1L) * rows_per_loc + (rows_per_loc - r + 1L))
    }))
    for (t in seq_len(par$d)) {
      i <- row_order[(t - 1L) %% par$rows + 1L]
      j <- row_loc[i]
      options <- setdiff((i - 1L) * per_row + seq_len(per_row), dummies)
      loc_plots <- plots[plots$location == j & !plots$index %in% dummies, , drop = FALSE]
      loc_block <- (loc_plots$replicate - 1L) * par$s + loc_plots$block
      base_comps <- prep_block_components(loc_plots$index, loc_block, 2L * par$s)
      cost <- vapply(options, function(idx) {
        keep <- loc_plots$index != idx
        comps <- prep_block_components(loc_plots$index[keep], loc_block[keep], 2L * par$s)
        1e6 * (comps > base_comps) + sum(deficit[plots$gblock[plots$index == idx]])
      }, numeric(1))
      pick <- options[which.min(cost + stats::runif(length(cost)) * 1e-3)]
      dummies <- c(dummies, pick)
      hit <- plots$gblock[plots$index == pick]
      deficit[hit] <- deficit[hit] + 1
    }
    plots <- plots[!plots$index %in% dummies, , drop = FALSE]
    plots$gblock <- NULL
  }

  real <- setdiff(seq_len(v_full), dummies)
  label <- integer(v_full)
  label[real] <- sample.int(length(real))
  plots$entry <- label[plots$index]

  # Location of the p-rep replication of each row.
  dup_loc <- if (method == "alpha") {
    (plots$row - 1L) %/% par$m + 1L
  } else {
    (plots$row - 1L) %/% (par$m %/% 2L) + 1L
  }

  groups <- vector("list", par$l)
  allocations <- vector("list", par$l)
  block_sizes <- vector("list", par$l)
  block_rep <- vector("list", par$l)
  for (j in seq_len(par$l)) {
    p <- plots[plots$location == j, , drop = FALSE]
    p$is_prep <- dup_loc[plots$location == j] == j
    s <- par$s
    # Randomize block order within each replicate.
    perm <- unlist(lapply(1:2, function(r) sample.int(s)))
    p$block_id <- (p$replicate - 1L) * s + perm[(p$replicate - 1L) * s + p$block]
    p <- p[sample.int(nrow(p)), , drop = FALSE]
    p$occurrence <- ifelse(p$is_prep, p$replicate, 1L)
    groups[[j]] <- sort(unique(p$entry[p$is_prep]))
    allocations[[j]] <- list(
      entry = p$entry,
      occurrence = as.integer(p$occurrence),
      block = as.integer(p$block_id),
      is_prep = p$is_prep,
      extra_components = prep_block_components(p$entry, p$block_id, 2L * s) - 1L
    )
    block_sizes[[j]] <- tabulate(p$block_id, nbins = 2L * s)
    block_rep[[j]] <- rep(1:2, each = s)
  }
  list(
    candidate = list(groups = groups, allocations = allocations),
    block_sizes = block_sizes,
    block_rep = block_rep
  )
}

# Average efficiency factor from integer entry and block ids.
prep_efficiency_ids <- function(entry, block, v, max_entries = 2000L) {
  if (v < 2L || v > max_entries) {
    return(NA_real_)
  }
  block <- as.integer(factor(block))
  incidence <- matrix(0, nrow = v, ncol = max(block))
  incidence[cbind(entry, block)] <- 1
  r <- rowSums(incidence)
  k <- colSums(incidence)
  c_mat <- diag(r) - (incidence %*% (t(incidence) / k))
  a_mat <- c_mat / sqrt(outer(r, r))
  e <- eigen(a_mat, symmetric = TRUE, only.values = TRUE)$values
  e <- sort(e, decreasing = TRUE)[seq_len(v - 1L)]
  if (min(e) < 1e-8) {
    return(0)
  }
  (v - 1) / sum(1 / e)
}

prep_candidate_efficiency <- function(candidate, v) {
  entry <- unlist(lapply(candidate$allocations, `[[`, "entry"))
  block <- unlist(lapply(seq_along(candidate$allocations), function(l) {
    paste(l, candidate$allocations[[l]]$block)
  }))
  prep_efficiency_ids(entry, block, v)
}

# Build an algebraic design: search arrays, then keep the candidate with the
# highest average efficiency factor among the optimal ones (second stage of
# Williams, Piepho and Whitaker, 2011, Section 2.3).
prep_algebraic_design <- function(method, inputs, n_prep_per_block, n_efficiency_candidates) {
  v <- length(inputs$treatments)
  l <- length(inputs$locations)
  max_block_size <- min(inputs$block_size)
  m <- if (is.null(n_prep_per_block)) {
    prep_choose_m(method, v, l, max_block_size)
  } else {
    arc_positive_integer(n_prep_per_block, "n_prep_per_block")
  }
  par <- prep_algebraic_parameters(method, v, l, m, check = TRUE)
  if (par$k_max > max_block_size) {
    stop(
      "With `n_prep_per_block` = ", m, " the blocks have ", par$k_max,
      " plots, more than `block_size` (", max_block_size, ").",
      call. = FALSE
    )
  }

  search <- prep_algebraic_search(
    method, par,
    n_starts = inputs$n_candidates,
    n_keep = max(1L, n_efficiency_candidates)
  )

  best <- NULL
  for (state in search$states) {
    built <- prep_algebraic_candidate(method, state, par)
    eff <- if (length(search$states) > 1L) prep_candidate_efficiency(built$candidate, v) else NA_real_
    if (is.null(best) || (!is.na(eff) && (is.na(best$eff) || eff > best$eff))) {
      best <- c(built, list(eff = eff))
    }
  }

  if (!search$optimal) {
    warning(
      "Method \"", method, "\" did not find a (0,1)-design with m = ", m,
      " (array of ", par$rows, " rows for s = ", par$s, " blocks per replicate); ",
      "some pairs of entries share two blocks. Try a smaller `n_prep_per_block` ",
      "or method \"search\".",
      call. = FALSE
    )
  }
  setup <- list(
    block_sizes = best$block_sizes,
    block_rep = best$block_rep
  )
  list(
    candidate = best$candidate,
    setup = setup,
    parameters = par,
    n_evaluated = search$n_evaluated
  )
}
