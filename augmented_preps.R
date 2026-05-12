# ============================================================
# Augmented p-rep design generator
# Search-based practical implementation
# ============================================================

# ------------------------------------------------------------
# Helper: split total as evenly as possible across groups
# ------------------------------------------------------------
balanced_counts <- function(total, n_groups) {
  base <- total %/% n_groups
  rem  <- total %% n_groups

  out <- rep(base, n_groups)
  if (rem > 0) {
    out[seq_len(rem)] <- out[seq_len(rem)] + 1
  }

  out
}


# ------------------------------------------------------------
# Helper: create block sizes differing by at most one plot
# ------------------------------------------------------------
balanced_block_sizes <- function(n_plots, max_block_size) {
  if (n_plots <= 0) stop("n_plots must be positive.")
  if (max_block_size <= 1) stop("max_block_size must be > 1.")

  n_blocks <- ceiling(n_plots / max_block_size)

  base <- n_plots %/% n_blocks
  rem  <- n_plots %% n_blocks

  sizes <- rep(base, n_blocks)
  if (rem > 0) {
    sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1
  }

  sizes
}


# ------------------------------------------------------------
# Assign which entries will be duplicated at each location
# Default: each entry is duplicated at exactly one location
# ------------------------------------------------------------
assign_prep_groups <- function(
  treatments,
  locations,
  prep_prop = NULL,
  n_duplicate_per_location = NULL
) {
  treatments <- as.character(treatments)
  locations  <- as.character(locations)

  v <- length(treatments)
  l <- length(locations)

  if (is.null(n_duplicate_per_location)) {

    # Default p-rep logic:
    # approximately 1/l of entries duplicated at each location,
    # so each treatment is duplicated at one location.
    if (is.null(prep_prop)) {
      ndup <- balanced_counts(v, l)
    } else {
      if (prep_prop <= 0 || prep_prop > 1) {
        stop("prep_prop must be > 0 and <= 1.")
      }

      ndup <- rep(round(v * prep_prop), l)

      if (sum(ndup) > v) {
        stop(
          "The requested prep_prop duplicates more entries than available. ",
          "Reduce prep_prop or use n_duplicate_per_location."
        )
      }
    }

  } else {

    if (length(n_duplicate_per_location) == 1) {
      ndup <- rep(as.integer(n_duplicate_per_location), l)
    } else {
      if (!is.null(names(n_duplicate_per_location))) {
        ndup <- as.integer(n_duplicate_per_location[locations])
      } else {
        ndup <- as.integer(n_duplicate_per_location)
      }
    }

    if (length(ndup) != l) {
      stop("n_duplicate_per_location must be length 1 or length equal to locations.")
    }

    if (any(is.na(ndup)) || any(ndup < 0)) {
      stop("n_duplicate_per_location contains NA or negative values.")
    }

    if (sum(ndup) > v) {
      stop("Total duplicated entries cannot be greater than number of treatments.")
    }
  }

  duplicated_entries <- sample(treatments, sum(ndup))

  groups <- vector("list", l)
  names(groups) <- locations

  start <- 1
  for (i in seq_along(locations)) {
    ni <- ndup[i]

    if (ni == 0) {
      groups[[i]] <- character(0)
    } else {
      groups[[i]] <- duplicated_entries[start:(start + ni - 1)]
      start <- start + ni
    }
  }

  repeat_map <- setNames(rep(NA_character_, v), treatments)

  for (loc in locations) {
    if (length(groups[[loc]]) > 0) {
      repeat_map[groups[[loc]]] <- loc
    }
  }

  list(
    groups = groups,
    repeat_map = repeat_map,
    n_duplicate_per_location = ndup
  )
}


# ------------------------------------------------------------
# Pairwise concurrence matrix
# Counts how many times each pair occurs in the same block
# ------------------------------------------------------------
concurrence_matrix <- function(design, treatments = NULL) {
  if (is.null(treatments)) {
    treatments <- sort(unique(as.character(design$accession)))
  }

  treatments <- as.character(treatments)

  mat <- matrix(
    0L,
    nrow = length(treatments),
    ncol = length(treatments),
    dimnames = list(treatments, treatments)
  )

  blocks <- split(
    as.character(design$accession),
    interaction(design$location, design$block, drop = TRUE)
  )

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


# ------------------------------------------------------------
# Add row/column coordinates for field layout
# ------------------------------------------------------------
add_field_coordinates <- function(design, locations, field_cols = NULL, serpentine = TRUE) {
  design <- design[order(
    match(design$location, locations),
    design$block,
    design$plot_in_block
  ), ]

  design$plot_number <- NA_integer_
  design$row <- NA_integer_
  design$col <- NA_integer_

  for (loc in locations) {
    idx <- which(design$location == loc)

    if (length(idx) == 0) next

    design$plot_number[idx] <- seq_along(idx)

    if (!is.null(field_cols)) {
      if (length(field_cols) == 1) {
        ncol_field <- field_cols
      } else {
        if (is.null(names(field_cols))) {
          stop("If field_cols has length > 1, it must be named by location.")
        }
        ncol_field <- field_cols[[loc]]
      }

      r <- ceiling(design$plot_number[idx] / ncol_field)
      c <- ((design$plot_number[idx] - 1) %% ncol_field) + 1

      if (serpentine) {
        c <- ifelse(r %% 2 == 0, ncol_field - c + 1, c)
      }

      design$row[idx] <- r
      design$col[idx] <- c
    }
  }

  row.names(design) <- NULL
  design
}


# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
make_augmented_prep_design <- function(
  treatments,
  locations,
  block_size = 8,
  prep_prop = NULL,
  n_duplicate_per_location = NULL,
  n_iter = 1000,
  seed = NULL,
  field_cols = NULL,
  serpentine = TRUE,
  verbose = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  treatments <- as.character(treatments)
  locations  <- as.character(locations)

  if (anyDuplicated(treatments)) {
    stop("treatments must be unique.")
  }

  if (anyDuplicated(locations)) {
    stop("locations must be unique.")
  }

  v <- length(treatments)

  best_design <- NULL
  best_score  <- Inf
  best_diag   <- NULL
  best_conc   <- NULL
  best_groups <- NULL

  for (iter in seq_len(n_iter)) {

    prep <- assign_prep_groups(
      treatments = treatments,
      locations = locations,
      prep_prop = prep_prop,
      n_duplicate_per_location = n_duplicate_per_location
    )

    repeat_map <- prep$repeat_map

    conc <- matrix(
      0L,
      nrow = v,
      ncol = v,
      dimnames = list(treatments, treatments)
    )

    all_rows <- list()
    row_counter <- 1L

    # --------------------------------------------------------
    # Build location by location
    # --------------------------------------------------------
    for (loc in locations) {

      repeated_here <- prep$groups[[loc]]

      # Every entry appears once at each location
      occ <- data.frame(
        accession = treatments,
        occurrence = 1L,
        repeated_at = unname(repeat_map[treatments]),
        is_extra_repeat = FALSE,
        stringsAsFactors = FALSE
      )

      # Entries assigned to this location get a second plot
      if (length(repeated_here) > 0) {
        occ_extra <- data.frame(
          accession = repeated_here,
          occurrence = 2L,
          repeated_at = loc,
          is_extra_repeat = TRUE,
          stringsAsFactors = FALSE
        )

        occ <- rbind(occ, occ_extra)
      }

      occ <- occ[sample(seq_len(nrow(occ))), , drop = FALSE]

      block_sizes <- balanced_block_sizes(
        n_plots = nrow(occ),
        max_block_size = block_size
      )

      n_blocks <- length(block_sizes)

      block_contents <- vector("list", n_blocks)
      block_rows     <- vector("list", n_blocks)
      remaining      <- block_sizes

      # ------------------------------------------------------
      # Greedy randomized placement
      # Chooses the block that creates the smallest concurrence penalty
      # ------------------------------------------------------
      for (i in seq_len(nrow(occ))) {

        acc <- occ$accession[i]
        possible_blocks <- which(remaining > 0)

        penalties <- vapply(
          possible_blocks,
          function(b) {
            members <- block_contents[[b]]
            if (is.null(members)) members <- character(0)

            # Very strong penalty if the same accession appears twice in one block
            duplicate_penalty <- if (acc %in% members) 1e8 else 0

            pair_penalty <- 0

            if (length(members) > 0) {
              members_unique <- unique(members)

              # Penalize pairs that have already occurred in another block
              previous_conc <- conc[acc, members_unique]

              pair_penalty <-
                sum(previous_conc >= 1) * 10000 +
                sum(previous_conc^2) * 100
            }

            fill_penalty <- length(members) / block_sizes[b]

            duplicate_penalty +
              pair_penalty +
              fill_penalty +
              runif(1, 0, 1)
          },
          numeric(1)
        )

        chosen_block <- possible_blocks[which.min(penalties)]

        members <- block_contents[[chosen_block]]
        if (is.null(members)) members <- character(0)

        # Update concurrence matrix
        if (length(members) > 0) {
          for (memb in unique(members)) {
            if (memb != acc) {
              conc[acc, memb] <- conc[acc, memb] + 1L
              conc[memb, acc] <- conc[memb, acc] + 1L
            }
          }
        }

        block_contents[[chosen_block]] <- c(members, acc)
        block_rows[[chosen_block]] <- c(block_rows[[chosen_block]], i)
        remaining[chosen_block] <- remaining[chosen_block] - 1L
      }

      # Convert blocks to data frame
      for (b in seq_len(n_blocks)) {
        idx <- block_rows[[b]]

        tmp <- occ[idx, , drop = FALSE]
        tmp$location <- loc
        tmp$block <- b
        tmp$block_size <- block_sizes[b]
        tmp$plot_in_block <- seq_len(nrow(tmp))
        tmp$block_uid <- paste0(loc, "_B", sprintf("%03d", b))

        tmp$entry_type <- ifelse(
          !is.na(tmp$repeated_at) & tmp$repeated_at == loc,
          "p_rep_entry_at_this_location",
          "single_entry_at_this_location"
        )

        all_rows[[row_counter]] <- tmp
        row_counter <- row_counter + 1L
      }
    }

    design <- do.call(rbind, all_rows)

    design <- design[, c(
      "location",
      "block",
      "block_uid",
      "block_size",
      "plot_in_block",
      "accession",
      "occurrence",
      "entry_type",
      "is_extra_repeat",
      "repeated_at"
    )]

    design <- add_field_coordinates(
      design = design,
      locations = locations,
      field_cols = field_cols,
      serpentine = serpentine
    )

    # --------------------------------------------------------
    # Diagnostics and score
    # --------------------------------------------------------
    upper_conc <- conc[upper.tri(conc)]

    n_pairs_gt1 <- sum(upper_conc > 1)
    total_excess <- sum(pmax(upper_conc - 1, 0))
    max_conc <- max(upper_conc)

    same_entry_same_block <- sum(
      vapply(
        split(
          design$accession,
          interaction(design$location, design$block, drop = TRUE)
        ),
        function(x) sum(duplicated(x)),
        numeric(1)
      )
    )

    replication_counts <- table(design$accession)

    # Lower is better
    score <-
      same_entry_same_block * 1e12 +
      n_pairs_gt1 * 1e7 +
      total_excess * 1e5 +
      max_conc * 1000 +
      sum((as.numeric(replication_counts) - mean(replication_counts))^2)

    if (score < best_score) {
      best_score  <- score
      best_design <- design
      best_conc   <- conc
      best_groups <- prep$groups

      best_diag <- data.frame(
        score = score,
        max_pairwise_concurrence = max_conc,
        n_pairs_with_concurrence_gt_1 = n_pairs_gt1,
        total_excess_concurrence = total_excess,
        same_entry_same_block = same_entry_same_block,
        min_replication = min(replication_counts),
        max_replication = max(replication_counts),
        mean_replication = mean(replication_counts),
        stringsAsFactors = FALSE
      )
    }

    if (verbose && iter %% max(1, floor(n_iter / 10)) == 0) {
      message("Iteration ", iter, "/", n_iter, " | best score = ", best_score)
    }
  }

  # Final summaries
  block_summary <- aggregate(
    accession ~ location + block + block_uid,
    data = best_design,
    FUN = length
  )
  names(block_summary)[names(block_summary) == "accession"] <- "n_plots"

  replication_summary <- data.frame(
    accession = names(table(best_design$accession)),
    total_plots = as.integer(table(best_design$accession)),
    repeated_at = unname(best_design$repeated_at[
      match(names(table(best_design$accession)), best_design$accession)
    ]),
    stringsAsFactors = FALSE
  )

  list(
    design = best_design,
    repeated_groups = best_groups,
    diagnostics = best_diag,
    block_summary = block_summary,
    replication_summary = replication_summary,
    concurrence_matrix = best_conc
  )
}
