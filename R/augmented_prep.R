#' Generate an augmented p-rep design
#'
#' Creates an augmented partially replicated (p-rep) design inspired by the
#' Williams, Piepho, and Whitaker approach. Instead of replicated checks, a
#' subset of test entries is replicated at each location, while all entries
#' appear once at each location and the additional replicated entries are placed
#' into incomplete blocks. Candidate designs are scored to reduce pairwise
#' block concurrence, especially avoiding pairs that occur together more than
#' once.
#'
#' @param treatments Character vector of entry names.
#' @param locations Character vector of location names.
#' @param block_size Maximum incomplete block size.
#' @param prep_prop Optional proportion of entries to duplicate per location.
#'   If `NULL`, entries are split across locations so each entry is duplicated
#'   at approximately one location.
#' @param n_duplicate_per_location Optional number of entries duplicated at
#'   each location. Can be a single integer, a vector matching `locations`, or a
#'   named vector indexed by location.
#' @param n_candidates Number of candidate designs to generate and score.
#' @param seed Optional random seed.
#' @param field_cols Optional number of field columns for row/column layout.
#'   Can be a single integer or a named vector indexed by location.
#' @param serpentine Logical. If `TRUE`, row/column coordinates use serpentine
#'   ordering within each location.
#' @param output_file Optional tab-delimited output file for the complete
#'   design.
#'
#' @return A list with design inputs, summary diagnostics, complete plot-level
#' design, Breedbase-compatible design columns, block summary, replication
#' summary, repeated-entry groups, and pairwise concurrence matrix. In the
#' complete design, `repeated_at` is the location where the entry receives its
#' extra p-rep plot, or `"not_repeated"` when the entry is not selected for
#' extra replication.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2", "Loc3"),
#'   block_size = 5,
#'   n_candidates = 5,
#'   seed = 123,
#'   field_cols = 5
#' )
#'
#' design$summary
#' head(design$design)
#'
#' @export
augmented_prep_design <- function(
    treatments,
    locations,
    block_size = 8,
    prep_prop = NULL,
    n_duplicate_per_location = NULL,
    n_candidates = 1000,
    seed = NULL,
    field_cols = NULL,
    serpentine = TRUE,
    output_file = NULL
) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  inputs <- prep_validate_inputs(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    n_candidates = n_candidates,
    field_cols = field_cols
  )

  best <- NULL
  for (i in seq_len(inputs$n_candidates)) {
    candidate <- prep_make_candidate(
      treatments = inputs$treatments,
      locations = inputs$locations,
      block_size = inputs$block_size,
      n_duplicate_per_location = inputs$n_duplicate_per_location,
      field_cols = inputs$field_cols,
      serpentine = serpentine
    )

    if (is.null(best) || candidate$summary$score < best$summary$score) {
      best <- candidate
    }
  }

  block_summary <- best$design |>
    dplyr::group_by(.data$location, .data$block, .data$block_uid) |>
    dplyr::summarise(n_plots = dplyr::n(), .groups = "drop")

  replication_summary <- best$design |>
    dplyr::group_by(.data$all_entries) |>
    dplyr::summarise(
      total_plots = dplyr::n(),
      locations = paste(unique(.data$location), collapse = ", "),
      repeated_at = paste(unique(.data$repeated_at[.data$repeated_at != "not_repeated"]), collapse = ", "),
      .groups = "drop"
    )
  replication_summary$repeated_at[replication_summary$repeated_at == ""] <- "not_repeated"

  breedbase_design <- best$design |>
    dplyr::select(dplyr::all_of(c("plots", "block", "all_entries", "rep", "is_control")))

  result <- list(
    input = inputs,
    summary = best$summary,
    design = best$design,
    breedbase_design = breedbase_design,
    block_summary = block_summary,
    replication_summary = replication_summary,
    repeated_groups = best$repeated_groups,
    concurrence_matrix = best$concurrence_matrix
  )
  class(result) <- c("augmented_prep_design", class(result))

  if (!is.null(output_file)) {
    utils::write.table(result$design, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  }

  result
}

#' Calculate pairwise concurrence for a p-rep design
#'
#' Counts how many times each pair of entries appears together in the same
#' incomplete block across all locations.
#'
#' @param design A p-rep design data frame from `augmented_prep_design()`.
#' @param treatments Optional treatment names defining matrix order.
#'
#' @return A symmetric integer matrix of pairwise block concurrences.
#'
#' @export
augmented_prep_concurrence <- function(design, treatments = NULL) {
  if (is.null(treatments)) {
    treatments <- sort(unique(as.character(design$all_entries)))
  }

  treatments <- as.character(treatments)
  mat <- matrix(0L, length(treatments), length(treatments), dimnames = list(treatments, treatments))
  blocks <- split(as.character(design$all_entries), interaction(design$location, design$block, drop = TRUE))

  for (x in blocks) {
    x <- unique(x)
    if (length(x) > 1) {
      cmb <- utils::combn(x, 2)
      for (j in seq_len(ncol(cmb))) {
        a <- cmb[1, j]
        b <- cmb[2, j]
        mat[a, b] <- mat[a, b] + 1L
        mat[b, a] <- mat[b, a] + 1L
      }
    }
  }

  mat
}

prep_validate_inputs <- function(treatments, locations, block_size, prep_prop,
                                 n_duplicate_per_location, n_candidates, field_cols) {
  treatments <- as.character(treatments)
  locations <- as.character(locations)
  treatments <- treatments[!is.na(treatments) & treatments != ""]
  locations <- locations[!is.na(locations) & locations != ""]
  block_size <- arc_positive_integer(block_size, "block_size")
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

  if (length(treatments) < 2) {
    stop("At least two treatments are required.", call. = FALSE)
  }
  if (length(locations) < 1) {
    stop("At least one location is required.", call. = FALSE)
  }
  if (anyDuplicated(treatments)) {
    stop("`treatments` must be unique.", call. = FALSE)
  }
  if (anyDuplicated(locations)) {
    stop("`locations` must be unique.", call. = FALSE)
  }
  if (!is.null(prep_prop) && (length(prep_prop) != 1 || is.na(prep_prop) ||
    prep_prop <= 0 || prep_prop > 1)) {
    stop("`prep_prop` must be a single number greater than 0 and no larger than 1.", call. = FALSE)
  }

  n_duplicate_per_location <- prep_normalize_duplicate_counts(
    treatments, locations, prep_prop, n_duplicate_per_location
  )
  field_cols <- prep_normalize_field_cols(field_cols, locations)

  list(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    n_candidates = n_candidates,
    field_cols = field_cols
  )
}

prep_normalize_duplicate_counts <- function(treatments, locations, prep_prop, n_duplicate_per_location) {
  n_treatments <- length(treatments)
  n_locations <- length(locations)

  if (is.null(n_duplicate_per_location)) {
    if (is.null(prep_prop)) {
      out <- prep_balanced_counts(n_treatments, n_locations)
    } else {
      out <- rep(round(n_treatments * prep_prop), n_locations)
    }
  } else if (length(n_duplicate_per_location) == 1) {
    out <- rep(as.integer(n_duplicate_per_location), n_locations)
  } else if (!is.null(names(n_duplicate_per_location))) {
    out <- as.integer(n_duplicate_per_location[locations])
  } else {
    out <- as.integer(n_duplicate_per_location)
  }

  if (length(out) != n_locations || any(is.na(out)) || any(out < 0)) {
    stop("`n_duplicate_per_location` must be non-negative and match locations.", call. = FALSE)
  }
  if (sum(out) > n_treatments) {
    stop("Total duplicated entries cannot exceed number of treatments.", call. = FALSE)
  }

  out
}

prep_normalize_field_cols <- function(field_cols, locations) {
  if (is.null(field_cols)) {
    return(NULL)
  }
  if (length(field_cols) == 1) {
    out <- rep(arc_positive_integer(field_cols, "field_cols"), length(locations))
    names(out) <- locations
    return(out)
  }
  if (is.null(names(field_cols))) {
    stop("If `field_cols` has length greater than 1, it must be named by location.", call. = FALSE)
  }
  out <- as.integer(field_cols[locations])
  if (any(is.na(out)) || any(out < 1)) {
    stop("`field_cols` must contain positive integers for every location.", call. = FALSE)
  }
  names(out) <- locations
  out
}

prep_make_candidate <- function(treatments, locations, block_size, n_duplicate_per_location,
                                field_cols, serpentine) {
  repeat_plan <- prep_assign_groups(treatments, locations, n_duplicate_per_location)
  concurrence <- matrix(0L, length(treatments), length(treatments), dimnames = list(treatments, treatments))
  rows <- list()
  row_counter <- 1L

  for (loc in locations) {
    occurrences <- prep_location_occurrences(treatments, repeat_plan$groups[[loc]], repeat_plan$repeat_map, loc)
    occurrences <- occurrences[sample(seq_len(nrow(occurrences))), , drop = FALSE]
    block_sizes <- prep_balanced_block_sizes(nrow(occurrences), block_size)
    allocation <- prep_allocate_blocks(occurrences, block_sizes, concurrence)
    concurrence <- allocation$concurrence

    for (block in seq_along(allocation$block_rows)) {
      tmp <- occurrences[allocation$block_rows[[block]], , drop = FALSE]
      tmp$location <- loc
      tmp$block <- block
      tmp$block_uid <- paste0(loc, "_B", sprintf("%03d", block))
      tmp$block_size <- block_sizes[block]
      tmp$plot_in_block <- seq_len(nrow(tmp))
      rows[[row_counter]] <- tmp
      row_counter <- row_counter + 1L
    }
  }

  design <- dplyr::bind_rows(rows)
  design <- prep_add_field_coordinates(design, locations, field_cols, serpentine)
  design <- prep_finalize_design(design)
  concurrence <- augmented_prep_concurrence(design, treatments)

  list(
    design = design,
    summary = prep_score_design(design, concurrence),
    concurrence_matrix = concurrence,
    repeated_groups = repeat_plan$groups
  )
}

prep_assign_groups <- function(treatments, locations, n_duplicate_per_location) {
  duplicated_entries <- sample(treatments, sum(n_duplicate_per_location))
  groups <- vector("list", length(locations))
  names(groups) <- locations
  start <- 1L

  for (i in seq_along(locations)) {
    n_here <- n_duplicate_per_location[i]
    if (n_here == 0) {
      groups[[i]] <- character()
    } else {
      groups[[i]] <- duplicated_entries[start:(start + n_here - 1L)]
      start <- start + n_here
    }
  }

  repeat_map <- stats::setNames(rep("not_repeated", length(treatments)), treatments)
  for (loc in locations) {
    repeat_map[groups[[loc]]] <- loc
  }

  list(groups = groups, repeat_map = repeat_map)
}

prep_location_occurrences <- function(treatments, repeated_here, repeat_map, location) {
  out <- data.frame(
    all_entries = treatments,
    occurrence = 1L,
    repeated_at = unname(repeat_map[treatments]),
    is_extra_repeat = 0L,
    stringsAsFactors = FALSE
  )

  if (length(repeated_here) > 0) {
    out <- rbind(
      out,
      data.frame(
        all_entries = repeated_here,
        occurrence = 2L,
        repeated_at = location,
        is_extra_repeat = 1L,
        stringsAsFactors = FALSE
      )
    )
  }

  out$entry_type <- ifelse(
    out$repeated_at == location,
    "p_rep_entry_at_this_location",
    "single_entry_at_this_location"
  )
  out
}

prep_allocate_blocks <- function(occurrences, block_sizes, concurrence) {
  block_contents <- vector("list", length(block_sizes))
  block_rows <- vector("list", length(block_sizes))
  remaining <- block_sizes

  for (i in seq_len(nrow(occurrences))) {
    entry <- occurrences$all_entries[i]
    possible_blocks <- which(remaining > 0)
    penalties <- vapply(possible_blocks, function(block) {
      members <- block_contents[[block]]
      if (is.null(members)) members <- character()
      duplicate_penalty <- if (entry %in% members) 1e8 else 0
      pair_penalty <- if (length(members) > 0) {
        previous <- concurrence[entry, unique(members)]
        sum(previous >= 1) * 10000 + sum(previous^2) * 100
      } else {
        0
      }
      duplicate_penalty + pair_penalty + length(members) / block_sizes[block] + stats::runif(1)
    }, numeric(1))
    chosen <- possible_blocks[which.min(penalties)]
    members <- block_contents[[chosen]]
    if (is.null(members)) members <- character()
    for (member in unique(members)) {
      if (member != entry) {
        concurrence[entry, member] <- concurrence[entry, member] + 1L
        concurrence[member, entry] <- concurrence[member, entry] + 1L
      }
    }
    block_contents[[chosen]] <- c(members, entry)
    block_rows[[chosen]] <- c(block_rows[[chosen]], i)
    remaining[chosen] <- remaining[chosen] - 1L
  }

  list(block_rows = block_rows, concurrence = concurrence)
}

prep_add_field_coordinates <- function(design, locations, field_cols, serpentine) {
  design <- design[order(match(design$location, locations), design$block, design$plot_in_block), ]
  design$row <- NA_integer_
  design$col <- NA_integer_

  for (loc in locations) {
    idx <- which(design$location == loc)
    plot_number <- seq_along(idx)
    if (!is.null(field_cols)) {
      n_cols <- field_cols[[loc]]
      row <- ceiling(plot_number / n_cols)
      col <- ((plot_number - 1L) %% n_cols) + 1L
      if (isTRUE(serpentine)) {
        col <- ifelse(row %% 2 == 0, n_cols - col + 1L, col)
      }
      design$row[idx] <- row
      design$col[idx] <- col
    }
  }

  rownames(design) <- NULL
  design
}

prep_finalize_design <- function(design) {
  design |>
    dplyr::arrange(.data$location, .data$block, .data$plot_in_block) |>
    dplyr::mutate(
      plots = dplyr::row_number(),
      rep = .data$occurrence,
      is_control = 0L
    ) |>
    dplyr::select(dplyr::all_of(c(
      "plots",
      "location",
      "row",
      "col",
      "block",
      "block_uid",
      "block_size",
      "plot_in_block",
      "all_entries",
      "occurrence",
      "entry_type",
      "is_extra_repeat",
      "repeated_at",
      "rep",
      "is_control"
    ))) |>
    tibble::as_tibble()
}

prep_score_design <- function(design, concurrence) {
  upper <- concurrence[upper.tri(concurrence)]
  block_entries <- split(design$all_entries, interaction(design$location, design$block, drop = TRUE))
  same_entry_same_block <- sum(vapply(block_entries, function(x) sum(duplicated(x)), numeric(1)))
  replication_counts <- table(design$all_entries)
  max_conc <- if (length(upper) == 0) 0 else max(upper)
  n_pairs_gt1 <- sum(upper > 1)
  total_excess <- sum(pmax(upper - 1, 0))
  score <- same_entry_same_block * 1e12 +
    n_pairs_gt1 * 1e7 +
    total_excess * 1e5 +
    max_conc * 1000 +
    sum((as.numeric(replication_counts) - mean(as.numeric(replication_counts)))^2)

  tibble::tibble(
    score = score,
    max_pairwise_concurrence = max_conc,
    n_pairs_with_concurrence_gt_1 = n_pairs_gt1,
    total_excess_concurrence = total_excess,
    same_entry_same_block = same_entry_same_block,
    min_replication = min(replication_counts),
    max_replication = max(replication_counts),
    mean_replication = mean(replication_counts),
    n_locations = dplyr::n_distinct(design$location),
    n_entries = dplyr::n_distinct(design$all_entries),
    n_plots = nrow(design),
    n_blocks = dplyr::n_distinct(design$block_uid)
  )
}

prep_balanced_counts <- function(total, n_groups) {
  base <- total %/% n_groups
  rem <- total %% n_groups
  out <- rep(base, n_groups)
  if (rem > 0) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  out
}

prep_balanced_block_sizes <- function(n_plots, max_block_size) {
  if (n_plots <= 0) stop("`n_plots` must be positive.", call. = FALSE)
  if (max_block_size <= 1) stop("`block_size` must be greater than 1.", call. = FALSE)
  prep_balanced_counts(n_plots, ceiling(n_plots / max_block_size))
}
