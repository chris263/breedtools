#' Suggest trial layouts from the number of entries and available field space
#'
#' Lists parameter options for the package's three design families, starting
#' with the information the breeder knows. Rows, columns and check counts
#' can be left unspecified. No randomization or optimization is performed.
#'
#' @param n_entries Number of test entries, excluding checks.
#' @param design Design families to consider: `"augmented_row_column"` (augmented
#'   row-column), `"prep"` (partially replicated), and `"resolvable"`
#'   (fully replicated row-column). Defaults to all three.
#' @param n_checks Distinct check genotypes to consider, a non-negative
#'   integer or vector. `NULL` tries 2 through 6 for augmented row-column designs and
#'   zero for the other families. p-rep uses no separate checks; it is
#'   excluded when only positive check counts are requested. Resolvable
#'   designs replicate checks just like the other genotypes.
#' @param rows,cols Known field dimensions, per location, including all
#'   replicates. Each can be a positive integer, a vector of dimensions to
#'   consider, or `NULL` to search automatically. Supplied dimensions are
#'   exact candidate sizes, not upper bounds. For p-rep, the last row may be
#'   partially occupied, but wholly empty final rows are not proposed.
#' @param n_locations Number of locations. Augmented row-column and resolvable layouts
#'   describe one field to randomize independently at each location. A
#'   p-rep option describes the combined multi-location design.
#' @param n_reps Full replication counts to try for resolvable designs
#'   (at least 2). Other design families use their own replication scheme.
#' @param block_size Maximum incomplete block sizes to try for p-rep.
#' @param prep_prop Proportions of entries to duplicate at each p-rep
#'   location. Counts use `round(n_entries * prep_prop)`, as in the generator.
#'   An entry may be duplicated at only one location, so the sum of these
#'   counts over locations must not exceed `n_entries`. `NULL` tries 15%,
#'   20% and 25%, and, for several locations, a balanced share of entries
#'   per location rounded down to respect this budget.
#' @param max_plots Largest field area (rows times columns) per location.
#'   Defaults to four times the number of genotypes (increased when more
#'   full replicates are requested), or the largest supplied field area
#'   when both dimensions are given.
#' @param max_dimension Largest automatically searched row or column
#'   dimension. Defaults to twice the square root of `max_plots`, with a
#'   minimum of 20. Explicit dimensions are not restricted by this bound.
#'   When neither dimension is supplied, fields with an aspect ratio above
#'   3 are excluded. Supply a dimension explicitly to consider longer fields.
#' @param min_error_df Minimum nominal error degrees of freedom. These are
#'   count-based values for the design model, assuming its effects are
#'   estimable, not degrees of freedom of a later spatial model.
#' @param n_options Maximum options returned per design family. The
#'   shortlist includes different check, full-replication or p-rep counts
#'   before adding further layouts with the same counts.
#' @param n_candidates Search candidates to pass to the augmented row-column
#'   and p-rep generators. Defaults to 200. This is a user-selected search
#'   effort, not a dimension inferred from entry counts.
#'
#' @details
#' This is a parameter planner. Candidates pass plot-count, block-size and
#' nominal error-df checks. Augmented row-column candidates also pass the check
#' coverage and row/column group conditions used by the generator. p-rep
#' candidates have enough repeated plots to connect their blocks and a
#' pair-count budget compatible with (0,1) concurrence. These necessary
#' conditions do not certify a randomized design's connectivity or
#' efficiency: generate the chosen option, inspect its summary and run
#' [check_design()]. The search is bounded and not exhaustive.
#'
#' Augmented row-column options use equal row and column groups. Every block
#' has at least `n_checks` rows and at least `n_checks` columns, leaving room
#' for one copy of each check in a distinct row and column. For example,
#' four checks require blocks of at least 4 rows by 4 columns. Every block
#' contains all checks at least once; extra positions are assigned to extra check
#' replication, so the number of test-entry positions is exactly
#' `n_entries`. Flexible allocations use `method = "optimal"`.
#'
#' Resolvable options may contain filler positions, which are reported in
#' `unused_plots`. Include both entries and checks in the generator's
#' `treatments`. Its native generator treats them all as treatments; if
#' checks need to be marked in Breedbase, set `is_control` in both the
#' generated `design` and `breedbase_design` tables before saving.
#'
#' Within each resource-count group, smaller field areas, fewer unused
#' positions and less elongated fields come first; Latin augmented row-column
#' constructions are preferred over searched allocations at the same area.
#' This ordering does not estimate statistical precision or choose the
#' biologically best design.
#'
#' @return A tibble with option id, design family, entry/check counts,
#'   field dimensions, location and replication counts, active plots per
#'   location and in total, field area, unused positions, check plots,
#'   repeated entries, blocks, nominal `error_df`, method and notes.
#'   Field dimensions are named `rows_in_field` and `cols_in_field`, matching
#'   [augmented_row_column_design()]. Block and replicate dimensions are
#'   included where applicable. `n_candidates` reports the chosen search
#'   effort. `checks_per_block` is a list-column with each augmented
#'   row-column option's check-allocation matrix, or `NULL` for the default
#'   complete check blocks. These columns support an explicit generator call.
#'   `generator` names the function to call; `arguments` is a list-column of
#'   its layout arguments. Add your own `treatments`, `controls` (augmented row-column)
#'   and `seed` before calling it with `do.call()`. p-rep arguments include
#'   placeholder location names, which can be replaced with your own names.
#'   Attribute `diagnostics` reports the number of shortlisted options and
#'   the reason for any excluded family. No matches return a typed empty
#'   tibble and a message suggesting which constraints to relax.
#'
#' @seealso [augmented_row_column_layouts()], [augmented_prep_parameters()],
#'   [augmented_row_column_design()], [augmented_prep_design()],
#'   [resolvable_row_column_design()]
#' @examples
#' # Only the number of entries is known:
#' options <- design_options(100)
#' options[, c("option", "design", "n_checks", "rows_in_field", "cols_in_field", "plots_per_location")]
#'
#' # A known field, four checks and exactly 108 test entries:
#' options <- design_options(108, design = "augmented_row_column", n_checks = 4,
#'                           rows = 12, cols = 12, n_candidates = 5)
#' selected <- options[1, ]
#' trial <- augmented_row_column_design(
#'   treatments = paste0("LINE_", seq_len(108)),
#'   controls = paste0("CHECK_", seq_len(selected$n_checks)),
#'   rows_in_field = selected$rows_in_field,
#'   cols_in_field = selected$cols_in_field,
#'   rows_per_block = selected$rows_per_block,
#'   cols_per_block = selected$cols_per_block,
#'   checks_per_block = selected$checks_per_block[[1]],
#'   method = selected$method,
#'   n_candidates = selected$n_candidates,
#'   seed = 2026
#' )
#'
#' # Fully replicated options with one known dimension:
#' design_options(42, design = "resolvable", cols = 7, n_reps = 3)
#' # Partial replication at three locations:
#' design_options(100, design = "prep", n_locations = 3)
#' @export
design_options <- function(n_entries, design = c("augmented_row_column", "prep", "resolvable"),
                           n_checks = NULL, rows = NULL, cols = NULL,
                           n_locations = 1L, n_reps = 2:3,
                           block_size = c(5L, 8L, 10L, 12L, 16L, 20L),
                           prep_prop = NULL,
                           max_plots = NULL, max_dimension = NULL,
                           min_error_df = 1L, n_options = 6L, n_candidates = 200L) {
  integers <- function(x, name, minimum = 1L, scalar = FALSE) {
    if (!is.numeric(x) || !length(x) || (scalar && length(x) != 1L) ||
        anyNA(x) || any(!is.finite(x) | x < minimum | x > .Machine$integer.max | x != floor(x))) {
      stop("`", name, "` must contain ", if (scalar) "one ", "integer",
           if (!scalar) "s", " >= ", minimum, ".", call. = FALSE)
    }
    sort(unique(as.integer(x)))
  }
  n_entries <- integers(n_entries, "n_entries", scalar = TRUE)
  if ("augmented" %in% design) {
    stop("Use `design = \"augmented_row_column\"` for augmented row-column designs.", call. = FALSE)
  }
  design <- match.arg(design, c("augmented_row_column", "prep", "resolvable"), several.ok = TRUE)
  if (!is.null(n_checks)) n_checks <- integers(n_checks, "n_checks", minimum = 0L)
  if (!is.null(rows)) rows <- integers(rows, "rows")
  if (!is.null(cols)) cols <- integers(cols, "cols")
  n_locations <- integers(n_locations, "n_locations", scalar = TRUE)
  n_reps <- integers(n_reps, "n_reps", minimum = 2L)
  block_size <- integers(block_size, "block_size", minimum = 2L)
  min_error_df <- integers(min_error_df, "min_error_df", scalar = TRUE)
  n_options <- integers(n_options, "n_options", scalar = TRUE)
  n_candidates <- integers(n_candidates, "n_candidates", scalar = TRUE)
  if (is.null(prep_prop)) {
    prep_prop <- c(0.15, 0.20, 0.25)
    if (n_locations > 1L && n_entries >= n_locations) {
      prep_prop <- unique(c(prep_prop, floor(n_entries / n_locations) / n_entries))
    }
  }
  if (!is.numeric(prep_prop) || !length(prep_prop) || anyNA(prep_prop) ||
      any(!is.finite(prep_prop) | prep_prop <= 0 | prep_prop > 1)) {
    stop("`prep_prop` must contain proportions greater than 0 and no larger than 1.", call. = FALSE)
  }
  max_checks <- if (is.null(n_checks)) 6L else max(n_checks)
  if (is.null(max_plots)) {
    multiplier <- if ("resolvable" %in% design) max(4, as.double(n_reps) + 1) else 4
    max_plots <- if (!is.null(rows) && !is.null(cols)) as.double(max(rows)) * max(cols) else
      multiplier * (as.double(n_entries) + max_checks)
  }
  max_plots <- integers(max_plots, "max_plots", scalar = TRUE)
  if (is.null(max_dimension)) max_dimension <- max(20, ceiling(2 * sqrt(max_plots)))
  max_dimension <- integers(max_dimension, "max_dimension", scalar = TRUE)
  fields <- do_option_fields(rows, cols, n_entries, max_plots, max_dimension)
  pieces <- list()
  diagnostics <- list()
  for (kind in unique(design)) {
    checks <- if (is.null(n_checks)) if (kind == "augmented_row_column") 2:6 else 0L else n_checks
    why <- "No layout within the supplied dimensions, field-area budget and error-df limit."
    choices <- switch(kind,
      augmented_row_column = if (any(checks > 0L)) {
        why <- paste0("No augmented row-column layout meets the supplied limits with at least ",
                       "n_checks rows and n_checks columns in every block.")
        do_option_augmented(n_entries, checks[checks > 0L], fields, n_locations, min_error_df, n_options)
      } else {
        why <- "Augmented row-column designs require at least one distinct check."
        NULL
      },
      prep = if (0L %in% checks && n_entries >= 2L) {
        do_option_prep(n_entries, fields, n_locations, block_size, prep_prop, min_error_df, n_options)
      } else {
        why <- "p-rep needs at least two entries and uses no separate checks; allow n_checks = 0."
        NULL
      },
      resolvable = do_option_resolvable(n_entries, checks, fields, n_locations, n_reps, min_error_df, n_options)
    )
    choices <- dplyr::bind_rows(choices)
    if (nrow(choices)) {
      if (kind %in% c("augmented_row_column", "prep")) {
        choices$n_candidates <- n_candidates
        choices$arguments <- lapply(choices$arguments, function(args) {
          args$n_candidates <- n_candidates
          args
        })
      }
      choices <- choices[order(choices$field_plots, choices$unused_plots,
        choices$method != "latin", abs(log(choices$rows_in_field / choices$cols_in_field)), -choices$error_df), ]
      group <- paste(choices$n_checks, choices$n_reps, choices$n_repeated, sep = "/")
      within_group <- ave(seq_len(nrow(choices)), group, FUN = seq_along)
      choices <- choices[head(order(within_group, seq_len(nrow(choices))), n_options), ]
      pieces[[kind]] <- choices
    }
    diagnostics[[kind]] <- tibble::tibble(design = kind, options = nrow(choices),
                                          reason = if (nrow(choices)) "" else why)
  }
  out <- dplyr::bind_rows(pieces)
  if (!nrow(out)) {
    out <- do_option_row()[FALSE, ]
    message("No design options found. Try other field dimensions, a larger `max_plots` or `max_dimension`, ",
            "more replication, or a smaller `min_error_df`. See attr(result, 'diagnostics').")
  }
  out$option <- seq_len(nrow(out))
  out <- out[, c("option", setdiff(names(out), "option"))]
  attr(out, "diagnostics") <- dplyr::bind_rows(diagnostics)
  out
}

do_option_fields <- function(rows, cols, n_entries, max_plots, max_dimension) {
  both_unknown <- is.null(rows) && is.null(cols)
  if (is.null(rows)) rows <- seq_len(min(max_dimension, max_plots))
  choices <- lapply(rows, function(r) {
    candidates <- cols
    if (is.null(candidates)) {
      lo <- max(1L, ceiling(n_entries / r))
      hi <- min(max_dimension, floor(max_plots / r))
      if (lo > hi) return(NULL)
      candidates <- seq.int(lo, hi)
    }
    keep <- as.double(r) * candidates >= n_entries & as.double(r) * candidates <= max_plots
    if (both_unknown) keep <- keep & pmax(r / candidates, candidates / r) <= 3
    if (!any(keep)) return(NULL)
    data.frame(rows = r, cols = candidates[keep])
  })
  out <- dplyr::bind_rows(choices)
  if (!nrow(out)) return(data.frame(rows = integer(), cols = integer()))
  out[order(as.double(out$rows) * out$cols, abs(log(out$rows / out$cols))), ]
}

do_option_row <- function(design = "", n_entries = 0L, n_checks = 0L, rows = 0L, cols = 0L,
                           n_locations = 1L, n_reps = NA_integer_, plots = 0L, check_plots = 0L,
                           repeated = NA_integer_, blocks = 0L, error_df = 0L, method = "",
                           generator = "", arguments = list(), note = "") {
  tibble::tibble(design = design, n_entries = n_entries, n_checks = n_checks,
    rows_in_field = rows, cols_in_field = cols, n_locations = n_locations, n_reps = n_reps,
    rows_per_block = arguments$rows_per_block %||% NA_integer_,
    cols_per_block = arguments$cols_per_block %||% NA_integer_,
    rows_per_rep = arguments$rows_per_rep %||% NA_integer_,
    cols_per_rep = arguments$cols_per_rep %||% NA_integer_,
    block_size = arguments$block_size %||% NA_integer_,
    n_candidates = NA_integer_, checks_per_block = list(arguments$checks_per_block),
    plots_per_location = plots, total_plots = as.double(plots) * n_locations,
    field_plots = as.double(rows) * cols, unused_plots = as.double(rows) * cols - plots,
    check_plots = check_plots, n_repeated = repeated, n_blocks = blocks,
    error_df = error_df, method = method, generator = generator,
    arguments = list(arguments), note = note)
}

do_option_divisors <- function(n) {
  small <- seq_len(floor(sqrt(n)))
  small <- small[n %% small == 0L]
  sort(unique(c(small, n %/% small)))
}

do_option_augmented <- function(n, checks, fields, locations, min_df, limit) {
  if (!length(checks)) return(NULL)
  choices <- list()
  seen_checks <- integer()
  for (i in seq_len(nrow(fields))) {
    r <- fields$rows[i]; c <- fields$cols[i]
    p <- as.double(r) * c - n
    if (p < r + c - 2 + min(checks) + min_df) next
    rb <- do_option_divisors(r); rb <- rb[rb >= 2L]
    cb <- do_option_divisors(c); cb <- cb[cb >= 2L]
    grid <- expand.grid(rb = rb, cb = cb, checks = checks)
    nr <- r / grid$rb; nc <- c / grid$cb; b <- nr * nc
    df <- p - grid$checks - (r - 1) - (c - 1) - (nr - 1) * (nc - 1)
    grid <- grid[grid$rb >= grid$checks & grid$cb >= grid$checks &
                   df >= min_df & floor(p / b) >= grid$checks &
                   ceiling(p / b) < grid$rb * grid$cb, , drop = FALSE]
    for (j in seq_len(nrow(grid))) {
      g <- grid[j, ]; nr <- r %/% g$rb; nc <- c %/% g$cb; b <- nr * nc
      pb <- matrix(p %/% b, nr, nc)
      extra <- matrix(0L, nr, nc)
      for (k in seq_len(p %% b)) {
        slots <- which(extra == 0L, arr.ind = TRUE)
        score <- rowSums(extra)[slots[, 1]] * (b + 1) + colSums(extra)[slots[, 2]]
        pos <- slots[which.min(score), ]
        extra[pos[1], pos[2]] <- 1L
      }
      pb <- pb + extra
      cap <- augmented_row_column_capacity(r, c, g$rb, g$cb, g$checks, checks_per_block = pb)
      if (!cap$every_row_has_check || !cap$every_col_has_check ||
          cap$min_rowgroup_slack < 0 || cap$min_colgroup_slack < 0) next
      regular <- all(pb == g$checks)
      latin <- regular && nc <= g$rb && nr <= g$cb && g$checks <= min(g$rb, g$cb)
      method <- if (latin) "latin" else "optimal"
      args <- list(rows_in_field = r, cols_in_field = c, rows_per_block = g$rb,
                    cols_per_block = g$cb, method = method)
      if (!regular) args$checks_per_block <- pb
      choices[[length(choices) + 1L]] <- do_option_row("augmented_row_column", n, g$checks, r, c, locations,
        1L, r * c, p, blocks = b, error_df = cap$error_df, method = method,
        generator = "augmented_row_column_design", arguments = args,
        note = if (regular) "Every check once per block." else
          "Extra check replication fills the field exactly; allocation needs optimization.")
      seen_checks <- union(seen_checks, g$checks)
    }
    # Search the smallest fields first; keep a pool larger than the shortlist.
    if (length(choices) >= 10L * limit && all(head(checks, limit) %in% seen_checks)) break
  }
  choices
}

do_option_prep <- function(n, fields, locations, sizes, proportions, min_df, limit) {
  choices <- list()
  for (d in unique(as.integer(round(n * proportions)))) {
    if (d < 1L || as.double(d) * locations > n) next
    plots <- n + d
    for (size in sizes[sizes <= n]) {
      blocks <- ceiling(plots / size)
      df <- d - blocks + 1L
      if (blocks < 2L || df < min_df) next
      sizes_here <- prep_balanced_counts(plots, blocks)
      if (sum(sizes_here * (sizes_here - 1) / 2) * locations > as.double(n) * (n - 1) / 2) next
      fits <- fields[fields$rows == ceiling(plots / fields$cols), , drop = FALSE]
      fits <- head(fits, 3L * limit)
      for (i in seq_len(nrow(fits))) {
        r <- fits$rows[i]; c <- fits$cols[i]
        args <- list(locations = paste0("Location_", seq_len(locations)), method = "search",
          block_size = size, n_duplicate_per_location = d, field_cols = c)
        choices[[length(choices) + 1L]] <- do_option_row("prep", n, 0L, r, c, locations,
          plots = plots, repeated = d, blocks = blocks, error_df = df, method = "search",
          generator = "augmented_prep_design", arguments = args,
          note = "Two plots for selected entries; one for the others. Check connectivity and concurrence after search.")
      }
    }
  }
  choices
}

do_option_resolvable <- function(n, checks, fields, locations, reps, min_df, limit) {
  choices <- list()
  for (check in checks) {
    v <- as.double(n) + check
    if (v < 2) next
    for (rep in reps) {
      for (layout in c("vertical", "horizontal")) {
        k <- if (layout == "vertical") fields$rows / rep else fields$rows
        s <- if (layout == "horizontal") fields$cols / rep else fields$cols
        df <- rep * v - rep * (k + s - 1) - (v - 1)
        keep <- k == floor(k) & s == floor(s) & k >= 2 & s >= 2 & k * s >= v & df >= min_df
        for (i in head(which(keep), 3L * limit)) {
          args <- list(n_reps = rep, rows_per_rep = k[i], cols_per_rep = s[i], replicate_layout = layout)
          choices[[length(choices) + 1L]] <- do_option_row("resolvable", n, check,
            fields$rows[i], fields$cols[i], locations, rep, rep * v, check * rep,
            repeated = v, blocks = rep * k[i], error_df = df[i], method = "search",
            generator = "resolvable_row_column_design", arguments = args,
            note = if (check) "Include checks in treatments; mark is_control in the generated field books." else
              "Every entry once per replicate; unused positions are fillers.")
        }
      }
    }
  }
  choices
}
