test_that("augmented_prep_design returns complete p-rep outputs", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(12)),
    locations = c("Loc1", "Loc2", "Loc3"),
    block_size = 5,
    n_candidates = 10,
    seed = 123,
    field_cols = 5
  )

  expect_s3_class(design, "augmented_prep_design")
  expect_equal(nrow(design$design), 48)
  expect_equal(nrow(design$breedbase_design), 48)
  expect_equal(design$summary$n_locations, 3)
  expect_equal(design$summary$n_entries, 12)
  expect_equal(design$summary$same_entry_same_block, 0)
  expect_equal(sort(as.integer(table(design$design$all_entries))), rep(4, 12))
  expect_false(anyNA(design$design$all_entries))
  expect_false(anyNA(design$design$repeated_at))
  expect_true(all(c("plots", "block", "all_entries", "rep", "is_control") %in% names(design$breedbase_design)))
})

test_that("augmented_prep_design supports explicit duplicate counts", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(10)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_duplicate_per_location = c(Loc1 = 2, Loc2 = 3),
    n_candidates = 5,
    seed = 456
  )

  expect_equal(nrow(design$design), 25)
  expect_equal(sum(design$design$is_extra_repeat), 5)
  expect_equal(design$input$n_duplicate_per_location, c(2, 3))
  expect_false(anyNA(design$design$all_entries))
  expect_false(anyNA(design$design$repeated_at))
  expect_true("not_repeated" %in% design$design$repeated_at)
})

test_that("augmented_prep_design supports unequal block sizes and replication levels", {
  expect_warning(design <- augmented_prep_design(
    treatments = paste0("G", seq_len(10)),
    locations = c("Loc1", "Loc2"),
    block_size = c(Loc1 = 4, Loc2 = 5),
    n_blocks_per_location = c(Loc1 = 4, Loc2 = 4),
    n_duplicate_per_location = c(Loc1 = 2, Loc2 = 3),
    replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3),
    n_candidates = 5,
    seed = 404,
    field_cols = c(Loc1 = 6, Loc2 = 8)
  ), "cannot all be connected")

  loc1 <- design$design[design$design$location == "Loc1", ]
  loc2 <- design$design[design$design$location == "Loc2", ]

  expect_equal(nrow(loc1), 12)
  expect_equal(nrow(loc2), 16)
  expect_equal(length(unique(loc1$block)), 4)
  expect_equal(length(unique(loc2$block)), 4)
  expect_lte(max(table(loc1$block)), 4)
  expect_lte(max(table(loc2$block)), 5)
  expect_true(3L %in% loc2$occurrence)
  expect_equal(max(loc1$block_size), 3)
  expect_equal(max(loc2$block_size), 4)
  expect_equal(design$summary$same_entry_same_block, 0)
})


test_that("augmented_prep_design supports rectangular block field layouts", {
  expect_warning(design <- augmented_prep_design(
    treatments = paste0("G", seq_len(10)),
    locations = c("Loc1", "Loc2", "Loc3"),
    block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 5),
    n_blocks_per_location = c(Loc1 = 2, Loc2 = 2, Loc3 = 4),
    n_duplicate_per_location = c(Loc1 = 0, Loc2 = 0, Loc3 = 4),
    n_candidates = 5,
    field_cols = 5,
    rows_per_block = c(Loc3 = 2),
    blocks_per_row = c(Loc3 = 2),
    seed = 707
  ), "cannot all be connected")

  loc3 <- design$design[design$design$location == "Loc3", ]
  block1 <- loc3[loc3$block == 1, ]
  block2 <- loc3[loc3$block == 2, ]
  block3 <- loc3[loc3$block == 3, ]

  expect_equal(max(loc3$row), 4)
  expect_equal(sort(unique(block1$row)), 1:2)
  expect_equal(sort(unique(block2$row)), 1:2)
  expect_equal(sort(unique(block3$row)), 3:4)
  expect_lte(max(block1$col), 2)
  expect_gte(min(block2$col), 3)
  expect_lte(max(block2$col), 4)
})

test_that("augmented_prep_concurrence is symmetric with zero diagonal", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(9)),
    locations = c("Loc1", "Loc2", "Loc3"),
    block_size = 4,
    n_candidates = 5,
    seed = 789
  )

  concurrence <- augmented_prep_concurrence(design$design)

  expect_true(isSymmetric(concurrence))
  expect_true(all(diag(concurrence) == 0))
  expect_equal(dim(concurrence), c(9, 9))
})


test_that("analyze_augmented_prep_design returns model outputs", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(8)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_duplicate_per_location = c(Loc1 = 2, Loc2 = 2),
    n_candidates = 5,
    seed = 505
  )

  design_data <- design$design
  design_data$yield <- 100 +
    as.numeric(factor(design_data$all_entries)) * 0.5 +
    ifelse(design_data$location == "Loc2", 3, 0) +
    as.numeric(factor(design_data$block_uid)) * 0.2 +
    seq_len(nrow(design_data)) * 0.01

  analysis <- analyze_augmented_prep_design(design_data, response = "yield")

  expect_s3_class(analysis, "augmented_prep_analysis")
  expect_s3_class(analysis$model, "lm")
  expect_equal(
    analysis$anova_formatted$Source,
    c(
      "Locations",
      "Blocks, nested within locations",
      "Genotypes",
      "Error",
      "Corrected total"
    )
  )
  expect_equal(analysis$anova_formatted$Df, c(1, 4, 7, 7, 19))
  expect_true(all(analysis$anova_formatted$Significance %in% c("", "*", "**")))
  expect_false(any(is.nan(unlist(analysis$anova_formatted[vapply(
    analysis$anova_formatted,
    is.numeric,
    logical(1)
  )]))))
  expect_equal(nrow(analysis$plot_results), nrow(design_data))
  expect_true(all(c("fitted", "residual") %in% names(analysis$plot_results)))
  expect_equal(nrow(analysis$genotype_summary), 8)
  expect_output(print(analysis), "Augmented p-rep design analysis")
})

test_that("analyze_augmented_prep_design accepts design objects and validates inputs", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(6)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_duplicate_per_location = c(Loc1 = 1, Loc2 = 1),
    n_candidates = 3,
    seed = 606
  )
  design$design$phenotype <- seq_len(nrow(design$design))

  analysis <- analyze_augmented_prep_design(design, response = "phenotype")
  expect_s3_class(analysis, "augmented_prep_analysis")

  expect_error(
    analyze_augmented_prep_design(
      data = design$design[, setdiff(names(design$design), "block")],
      response = "phenotype"
    ),
    "Missing required column"
  )
})

test_that("plot_augmented_prep_design draws design objects and data frames", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(10)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_duplicate_per_location = 2,
    n_candidates = 3,
    seed = 202,
    field_cols = 4
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plotted_object <- plot_augmented_prep_design(
    design,
    label = "none",
    fill = "entry_type",
    legend = FALSE
  )
  plotted_data <- plot_augmented_prep_design(
    design$design,
    locations = "Loc1",
    label = "all_entries",
    fill = "is_extra_repeat",
    show_block_labels = TRUE,
    legend = FALSE
  )
  grDevices::dev.off()

  expect_s3_class(plotted_object, "tbl_df")
  expect_s3_class(plotted_data, "tbl_df")
  expect_true(all(c(".plot_row", ".plot_col", ".fill_value", ".panel_index") %in% names(plotted_object)))
  expect_true(file.exists(plot_file))
})

test_that("plot_augmented_prep_design requires field coordinates", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(8)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_duplicate_per_location = 2,
    n_candidates = 3,
    seed = 303
  )

  expect_error(
    plot_augmented_prep_design(design),
    "field_cols"
  )
})

test_that("augmented_prep_design validates inputs", {
  expect_error(
    augmented_prep_design(
      treatments = c("G1", "G1"),
      locations = "Loc1"
    ),
    "unique"
  )

  expect_error(
    augmented_prep_design(
      treatments = paste0("G", seq_len(3)),
      locations = c("Loc1", "Loc2"),
      n_duplicate_per_location = c(2, 2)
    ),
    "cannot exceed"
  )
})

test_that("augmented_prep_design can write tab-delimited output", {
  path <- tempfile(fileext = ".design")
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(6)),
    locations = c("Loc1", "Loc2"),
    block_size = 4,
    n_candidates = 3,
    seed = 101,
    output_file = path
  )

  loaded <- read_breedbase_file(path)

  expect_equal(nrow(loaded), nrow(design$design))
  expect_equal(names(loaded), names(design$design))
})

test_that("analyze_augmented_prep_design handles a single location", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(40)),
    locations = "Station",
    block_size = 6,
    prep_prop = 0.25,
    n_candidates = 5,
    seed = 11
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_prep_design(design_data, response = "yield")

  ref <- stats::anova(stats::lm(
    yield ~ factor(block) + factor(all_entries),
    data = design_data
  ))
  expect_equal(analysis$anova_formatted$Df, c(0, ref$Df, nrow(design_data) - 1))
  expect_equal(analysis$anova_formatted$`Sum Sq`[2:4], ref$`Sum Sq`, tolerance = 1e-8)
})

test_that("augmented p-rep ANOVA df come from the observed plots", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(60)),
    locations = c("A", "B"),
    block_size = 8,
    prep_prop = 0.3,
    n_candidates = 5,
    seed = 21
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  design_data$yield[c(2, 15, 40, 77, 100)] <- NA_real_

  analysis <- analyze_augmented_prep_design(design_data, response = "yield")

  observed <- design_data[!is.na(design_data$yield), ]
  ref <- stats::anova(stats::lm(
    yield ~ location + interaction(location, block, drop = TRUE) + all_entries,
    data = observed
  ))
  expect_equal(analysis$anova_formatted$Df, c(ref$Df, nrow(observed) - 1))
  expect_equal(analysis$anova_formatted$`Sum Sq`[1:4], ref$`Sum Sq`, tolerance = 1e-8)
  expect_equal(analysis$diagnostics$anova_error_df, ref$Df[4])
  expect_gt(analysis$diagnostics$design_error_df, analysis$diagnostics$anova_error_df)
})

test_that("augmented_prep_design stops early once the optimal score is reached", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(300)),
    locations = c("A", "B", "C"),
    block_size = 12,
    n_candidates = 1000,
    seed = 1
  )

  expect_lt(design$summary$n_candidates_evaluated, 1000)
  expect_equal(design$summary$max_pairwise_concurrence, 1)
  expect_equal(design$summary$n_pairs_with_concurrence_gt_1, 0)
  expect_equal(design$summary$same_entry_same_block, 0)
  expect_identical(
    unname(design$concurrence_matrix),
    unname(augmented_prep_concurrence(design$design, rownames(design$concurrence_matrix)))
  )
})

test_that("augmented p-rep adjusted means match a full lm fit", {
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(60)),
    locations = c("A", "B"),
    block_size = 8,
    prep_prop = 0.3,
    n_candidates = 5,
    seed = 21
  )
  set.seed(4)
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  design_data$yield[c(2, 15, 40)] <- NA_real_

  analysis <- analyze_augmented_prep_design(design_data, response = "yield")

  observed <- design_data[!is.na(design_data$yield), ]
  observed$genotype <- factor(observed$all_entries)
  observed$location <- factor(observed$location)
  observed$block_nested <- interaction(observed$location, observed$block, drop = TRUE)
  fit <- stats::lm(yield ~ location + block_nested + genotype, data = observed)
  ref <- reference_adjusted_means(fit, observed)
  got <- analysis$adjusted_means[match(ref$genotype, analysis$adjusted_means$genotype), ]

  expect_equal(got$adjusted_mean, ref$mean, tolerance = 1e-8)
  expect_equal(got$std_error, ref$se, tolerance = 1e-8)
  expect_equal(analysis$adjusted_means$rank[1], 1L)
  expect_gt(analysis$diagnostics$mean_sed, 0)
})

# --- Williams, Piepho and Whitaker (2011) scenarios --------------------------

prep_location_connected <- function(design) {
  vapply(split(design, design$location), function(x) {
    prep_block_components(x$all_entries, x$block, length(unique(x$block))) == 1L
  }, logical(1))
}

test_that("paper Example 1 (v = 20, l = 2, k = 3) gives a connected (0,1)-design", {
  for (resolvable in c(FALSE, TRUE)) {
    design <- augmented_prep_design(
      treatments = paste0("G", 1:20),
      locations = c("L1", "L2"),
      block_size = 3,
      resolvable = resolvable,
      seed = 1
    )
    d <- design$design
    expect_equal(nrow(d), 60)
    expect_equal(design$summary$n_blocks, 20)
    expect_equal(design$summary$max_pairwise_concurrence, 1)
    expect_equal(design$summary$n_disconnected_locations, 0)
    expect_true(all(prep_location_connected(d)))
    # Each block holds m = 2 p-rep plots, as in Table 2 of the paper.
    prep_plots <- d[d$entry_type == "p_rep_entry_at_this_location", ]
    expect_true(all(table(interaction(prep_plots$location, prep_plots$block, drop = TRUE)) == 2))
    expect_gt(design$summary$efficiency_factor, 0)
    expect_lte(design$summary$efficiency_factor, 1)
  }
})

test_that("resolvable designs put every p-rep entry once in each replicate", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:42),
    locations = c("L1", "L2", "L3"),
    block_size = 4,
    resolvable = TRUE,
    seed = 2
  )
  d <- design$design
  expect_true(design$summary$resolvable)
  expect_equal(design$summary$max_pairwise_concurrence, 1)
  prep_plots <- d[d$entry_type == "p_rep_entry_at_this_location", ]
  per_rep <- table(paste(prep_plots$location, prep_plots$all_entries), prep_plots$replicate)
  expect_true(all(per_rep == 1))
  # Same number of blocks in every replicate, blocks nested in replicates.
  blocks <- unique(d[, c("location", "replicate", "block")])
  expect_true(all(table(blocks$location, blocks$replicate) == 7))
  expect_equal(d$rep, d$replicate)
})

test_that("resolvable designs validate the number of blocks", {
  expect_error(
    augmented_prep_design(
      treatments = paste0("G", 1:20),
      locations = c("L1", "L2"),
      block_size = 4,
      n_blocks_per_location = 7,
      resolvable = TRUE
    ),
    "multiple"
  )
})

test_that("the six-location maize case of the paper reaches a (0,1)-design", {
  design <- augmented_prep_design(
    treatments = sprintf("G%03d", 1:120),
    locations = paste0("L", 1:6),
    block_size = 8,
    n_candidates = 100,
    seed = 1
  )
  expect_equal(design$summary$max_pairwise_concurrence, 1)
  expect_equal(design$summary$n_pairs_with_concurrence_gt_1, 0)
  expect_equal(design$summary$n_disconnected_locations, 0)
  expect_lte(max(design$design$block_size) - min(design$design$block_size), 1)
})

test_that("unequal numbers of p-rep entries per location are supported", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:42),
    locations = c("L1", "L2", "L3"),
    block_size = c(L1 = 5, L2 = 6, L3 = 5),
    n_duplicate_per_location = c(L1 = 12, L2 = 18, L3 = 12),
    seed = 3
  )
  expect_equal(as.integer(table(design$design$location)), 42L + c(12L, 18L, 12L))
  expect_equal(design$summary$n_disconnected_locations, 0)
  expect_equal(design$summary$max_pairwise_concurrence, 1)
})

test_that("single-location p-rep designs are connected and fully estimable", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:60),
    locations = "S",
    block_size = 6,
    prep_prop = 0.2,
    seed = 1
  )
  expect_equal(design$summary$n_disconnected_locations, 0)
  expect_gt(design$summary$efficiency_factor, 0)
  d <- design$design
  set.seed(1)
  d$y <- stats::rnorm(nrow(d))
  analysis <- analyze_augmented_prep_design(d, "y")
  expect_false(anyNA(analysis$adjusted_means$adjusted_mean))
})

test_that("too few p-rep plots to connect the blocks gives a warning", {
  expect_warning(
    design <- augmented_prep_design(
      treatments = paste0("G", 1:60),
      locations = "S",
      block_size = 6,
      prep_prop = 0.05,
      seed = 1
    ),
    "cannot all be connected"
  )
  expect_equal(design$summary$n_disconnected_locations, 1)
  expect_equal(design$summary$efficiency_factor, 0)
})

test_that("analysis of a resolvable design fits replicates within locations", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:42),
    locations = c("A", "B", "C"),
    block_size = 4,
    resolvable = TRUE,
    seed = 3
  )
  d <- design$design
  set.seed(1)
  d$y <- stats::rnorm(nrow(d))
  d$y[c(5, 60)] <- NA
  analysis <- analyze_augmented_prep_design(d, "y")

  expect_equal(
    analysis$anova_formatted$Source,
    c("Locations", "Replicates, nested within locations",
      "Blocks, nested within replicates", "Genotypes", "Error", "Corrected total")
  )
  observed <- d[!is.na(d$y), ]
  ref <- stats::anova(stats::lm(
    y ~ location + interaction(location, replicate) +
      interaction(location, replicate, block) + all_entries,
    data = observed
  ))
  expect_equal(analysis$anova_formatted$Df[1:5], ref$Df)
  expect_equal(analysis$anova_formatted$`Sum Sq`[1:5], ref$`Sum Sq`, tolerance = 1e-8)

  ignored <- analyze_augmented_prep_design(d, "y", replicate_col = NULL)
  expect_false("Replicates, nested within locations" %in% ignored$anova_formatted$Source)
})

test_that("genotype-by-location analysis matches lm and recovers mean effects", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:30),
    locations = c("A", "B", "C"),
    block_size = 6,
    prep_prop = 0.3,
    seed = 5
  )
  d <- design$design
  set.seed(2)
  ge <- stats::setNames(stats::rnorm(30), paste0("G", 1:30))
  gl <- matrix(stats::rnorm(90), 30, 3, dimnames = list(names(ge), c("A", "B", "C")))
  be <- stats::rnorm(100)
  d$y <- 10 + ge[d$all_entries] + gl[cbind(d$all_entries, d$location)] +
    be[as.integer(factor(d$block_uid))] + stats::rnorm(nrow(d), sd = 0.3)

  analysis <- analyze_augmented_prep_design(d, "y", genotype_by_location = TRUE)
  expect_true("Genotype x location" %in% analysis$anova_formatted$Source)
  ref <- stats::anova(stats::lm(
    y ~ location + interaction(location, block) + all_entries + location:all_entries,
    data = d
  ))
  n_terms <- nrow(ref)
  expect_equal(analysis$anova_formatted$Df[seq_len(n_terms)], ref$Df)
  expect_equal(analysis$anova_formatted$`Sum Sq`[seq_len(n_terms)], ref$`Sum Sq`, tolerance = 1e-8)
  expect_equal(
    analysis$diagnostics$design_error_df,
    analysis$anova_df$display_df[analysis$anova_df$term == "Residuals"]
  )

  # Without noise, differences of adjusted means equal the location-weighted
  # genotype + genotype-by-location effects.
  d$z <- d$y - (d$y - (10 + ge[d$all_entries] + gl[cbind(d$all_entries, d$location)] +
    be[as.integer(factor(d$block_uid))]))
  exact <- suppressWarnings(analyze_augmented_prep_design(d, "z", genotype_by_location = TRUE))
  w <- as.numeric(table(d$location)[c("A", "B", "C")]) / nrow(d)
  truth <- ge + as.vector(gl %*% w)
  got <- exact$adjusted_means$adjusted_mean[match(names(ge), exact$adjusted_means$genotype)]
  expect_equal(diff(got), unname(diff(truth)), tolerance = 1e-8)
})

test_that("each location of a p-rep design can be analysed on its own", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:60),
    locations = c("A", "B"),
    block_size = 6,
    prep_prop = 0.3,
    resolvable = TRUE,
    seed = 8
  )
  d <- design$design
  set.seed(3)
  d$y <- stats::rnorm(nrow(d))
  for (loc in c("A", "B")) {
    analysis <- analyze_augmented_prep_design(d[d$location == loc, ], "y")
    expect_false(anyNA(analysis$adjusted_means$adjusted_mean))
    expect_gt(analysis$diagnostics$residual_df, 0)
  }
})

# --- Construction methods ---------------------------------------------------

prep_check_paper_structure <- function(design, m) {
  d <- design$design
  # Every entry is duplicated at exactly one location.
  counts <- table(d$all_entries)
  expect_true(all(counts == length(unique(d$location)) + 1L))
  # Resolvable: each p-rep entry once in each of the two replicates.
  prep_plots <- d[d$entry_type == "p_rep_entry_at_this_location", ]
  per_rep <- table(paste(prep_plots$location, prep_plots$all_entries), prep_plots$replicate)
  expect_true(all(per_rep == 1))
  # At most m p-rep plots in each block.
  expect_lte(max(table(interaction(prep_plots$location, prep_plots$block, drop = TRUE))), m)
  expect_equal(design$summary$max_pairwise_concurrence, 1)
  expect_equal(design$summary$n_disconnected_locations, 0)
  expect_equal(design$summary$same_entry_same_block, 0)
}

test_that("augmented_prep_parameters reproduces the paper's settings", {
  ex1 <- augmented_prep_parameters(20, 2, max_block_size = 3)
  row <- ex1[ex1$method == "alpha" & ex1$m == 2, ]
  expect_equal(row$s, 5L)
  expect_equal(row$block_size, 3L)
  expect_equal(row$dummies, 0L)

  ex3 <- augmented_prep_parameters(42, 3, max_block_size = 4)
  row <- ex3[ex3$method == "alpha2" & ex3$m == 2, ]
  expect_equal(row$s, 7L)
  expect_equal(row$block_size, 4L)
  expect_equal(row$blocks_per_location, 14L)
  expect_true(row$feasible)

  maize <- augmented_prep_parameters(120, 6, max_block_size = 8)
  expect_false(maize$feasible[maize$method == "alpha2" & maize$m == 2])
})

test_that("method alpha builds paper Example 1 as a resolvable (0,1)-design", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:20),
    locations = c("L1", "L2"),
    method = "alpha",
    block_size = 3,
    seed = 1
  )
  expect_equal(design$summary$method, "alpha")
  expect_equal(design$summary$m, 2L)
  expect_equal(design$summary$s, 5L)
  expect_equal(nrow(design$design), 60)
  expect_true(all(design$design$block_size == 3))
  expect_true(all(design$input$n_duplicate_per_location == 10))
  prep_check_paper_structure(design, m = 2)
})

test_that("method alpha removes double concurrences when k' > s (paper Example 2)", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:12),
    locations = c("L1", "L2"),
    method = "alpha",
    block_size = 3,
    seed = 2
  )
  expect_equal(design$summary$s, 3L)
  prep_check_paper_structure(design, m = 2)
})

test_that("method alpha drops dummy entries without disconnecting locations", {
  design <- augmented_prep_design(
    treatments = sprintf("G%03d", 1:150),
    locations = paste0("L", 1:6),
    method = "alpha",
    block_size = 8,
    n_candidates = 20,
    seed = 1
  )
  expect_gt(design$summary$n_dummies, 0)
  expect_equal(design$summary$n_entries, 150)
  expect_lte(max(design$design$block_size) - min(design$design$block_size), 1)
  prep_check_paper_structure(design, m = 2)
})

test_that("method alpha2 builds paper Example 3 (Table 5)", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:42),
    locations = c("L1", "L2", "L3"),
    method = "alpha2",
    block_size = 4,
    seed = 3
  )
  expect_equal(design$summary$s, 7L)
  expect_equal(design$summary$n_blocks, 42)
  expect_true(all(design$design$block_size == 4))
  prep_check_paper_structure(design, m = 2)
  expect_gt(design$summary$efficiency_factor, 0)
})

test_that("algebraic methods validate their parameters", {
  expect_error(
    augmented_prep_design(paste0("G", 1:120), paste0("L", 1:6), method = "alpha2", block_size = 8),
    "No feasible"
  )
  expect_error(
    augmented_prep_design(paste0("G", 1:120), paste0("L", 1:6), method = "alpha2",
                          block_size = 8, n_prep_per_block = 2),
    "odd number of blocks"
  )
  expect_error(
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha", prep_prop = 0.3),
    "cannot be used"
  )
  expect_error(
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha",
                          replications_per_repeated_entry = 3),
    "replications_per_repeated_entry = 2"
  )
  expect_error(
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), n_prep_per_block = 2),
    "only used"
  )
  expect_error(
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha",
                          block_size = 3, n_prep_per_block = 4),
    "more than `block_size`"
  )
})

test_that("algebraic designs can be laid out, plotted and analysed", {
  design <- augmented_prep_design(
    treatments = paste0("G", 1:42),
    locations = c("A", "B", "C"),
    method = "alpha",
    block_size = 4,
    field_cols = 4,
    seed = 5
  )
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plotted <- plot_augmented_prep_design(design, fill = "replicate", label = "replicate", legend = FALSE)
  grDevices::dev.off()
  expect_s3_class(plotted, "tbl_df")

  d <- design$design
  set.seed(6)
  d$y <- stats::rnorm(nrow(d))
  analysis <- analyze_augmented_prep_design(d, "y")
  expect_true("Replicates, nested within locations" %in% analysis$anova_formatted$Source)
  expect_false(anyNA(analysis$adjusted_means$adjusted_mean))
})

test_that("the most efficient of several optimal candidates is returned", {
  one <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha",
                               block_size = 4, n_efficiency_candidates = 1, seed = 9)
  many <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha",
                                block_size = 4, n_efficiency_candidates = 10, seed = 9)
  expect_equal(one$summary$max_pairwise_concurrence, 1)
  expect_gte(many$summary$efficiency_factor, one$summary$efficiency_factor - 1e-12)
})

test_that("random genotype-by-location analysis matches lme4", {
  skip_if_not_installed("lme4")
  design <- augmented_prep_design(paste0("G", 1:60), c("A", "B", "C"), method = "alpha",
                                  block_size = 6, seed = 4)
  d <- design$design
  set.seed(1)
  g <- stats::setNames(stats::rnorm(60, sd = 2), paste0("G", 1:60))
  gxl <- interaction(d$all_entries, d$location, drop = TRUE)
  d$y <- 50 + g[d$all_entries] + stats::rnorm(nlevels(gxl), sd = 1.5)[as.integer(gxl)] +
    stats::rnorm(nrow(d))

  analysis <- analyze_augmented_prep_design(d, "y", genotype_by_location = "random")
  expect_true(inherits(analysis$mixed_model, "lmerMod"))
  expect_equal(analysis$variance_components$component, c("genotype x location", "residual"))
  expect_output(print(analysis), "REML variance components")

  d$gen <- factor(d$all_entries)
  d$bn <- interaction(d$location, d$replicate, d$block, drop = TRUE)
  d$rn <- interaction(d$location, d$replicate, drop = TRUE)
  d$gxl <- gxl
  ref <- lme4::lmer(y ~ location + rn + bn + gen + (1 | gxl), data = d,
                    control = lme4::lmerControl(check.rankX = "silent.drop.cols"))
  fe <- lme4::fixef(ref)
  ref_contrast <- unname(c(0, fe[paste0("gen", levels(d$gen)[-1])]))
  got <- analysis$adjusted_means$adjusted_mean[match(levels(d$gen), analysis$adjusted_means$genotype)]
  expect_equal(got - got[1], ref_contrast, tolerance = 1e-5)
  expect_equal(
    analysis$variance_components$variance,
    as.data.frame(lme4::VarCorr(ref))$vcov,
    tolerance = 1e-5
  )

  expect_error(analyze_augmented_prep_design(d, "y", genotype_by_location = "maybe"), "must be")
})

# --- Two-stage analysis -------------------------------------------------------

two_stage_trial <- function(sd_by_location = c(A = 1, B = 1, C = 1), trend = 0, seed = 1) {
  design <- augmented_prep_design(
    treatments = sprintf("G%02d", 1:60),
    locations = names(sd_by_location),
    method = "alpha",
    block_size = 6,
    field_cols = 9,
    seed = 1
  )
  d <- design$design
  set.seed(seed)
  g <- stats::setNames(stats::rnorm(60, sd = 2), sprintf("G%02d", 1:60))
  gxl <- interaction(d$all_entries, d$location, drop = TRUE)
  d$y <- 50 + g[d$all_entries] + stats::rnorm(nlevels(gxl))[as.integer(gxl)] +
    trend * (sin(d$row / 2) + cos(d$col / 3)) +
    stats::rnorm(nrow(d), sd = sd_by_location[d$location])
  list(data = d, truth = g)
}

test_that("two-stage stage 2 matches brute-force REML and GLS", {
  trial <- two_stage_trial()
  fit <- analyze_augmented_prep_two_stage(trial$data, "y")
  expect_s3_class(fit, "augmented_prep_two_stage")
  expect_output(print(fit), "Two-stage augmented p-rep analysis")

  s1 <- lapply(c("A", "B", "C"), function(loc) {
    a <- analyze_augmented_prep_design(trial$data[trial$data$location == loc, ], "y")
    am <- a$adjusted_means[!is.na(a$adjusted_means$adjusted_mean), ]
    list(loc = loc, means = stats::setNames(am$adjusted_mean, am$genotype),
         vcov = a$adjusted_means_vcov[am$genotype, am$genotype])
  })
  st <- do.call(rbind, lapply(s1, function(r) data.frame(loc = r$loc, g = names(r$means), y = r$means)))
  st$loc <- factor(st$loc)
  st$g <- factor(st$g, levels = sort(unique(st$g)))
  X <- stats::model.matrix(~ loc + g, st)
  S <- matrix(0, nrow(st), nrow(st))
  off <- 0
  for (r in s1) {
    idx <- off + seq_along(r$means)
    S[idx, idx] <- r$vcov
    off <- off + length(idx)
  }
  reml <- function(s2) {
    V <- S + diag(s2, nrow(S))
    Vi <- solve(V)
    XtVX <- t(X) %*% Vi %*% X
    b <- solve(XtVX, t(X) %*% Vi %*% st$y)
    res <- st$y - X %*% b
    -0.5 * (determinant(V)$modulus + determinant(XtVX)$modulus + t(res) %*% Vi %*% res)
  }
  opt <- stats::optimize(reml, c(0, 20), maximum = TRUE)
  expect_equal(fit$variance_components$variance, opt$maximum, tolerance = 1e-3)

  Vi <- solve(S + diag(fit$variance_components$variance, nrow(S)))
  XtVX <- t(X) %*% Vi %*% X
  b <- solve(XtVX, t(X) %*% Vi %*% st$y)
  L <- matrix(0, 60, ncol(X))
  L[, 1] <- 1
  L[, grep("^loc", colnames(X))] <- 1 / 3
  for (j in 2:60) L[j, match(paste0("g", levels(st$g)[j]), colnames(X))] <- 1
  got <- fit$adjusted_means[match(levels(st$g), fit$adjusted_means$genotype), ]
  expect_equal(got$adjusted_mean, drop(L %*% b), tolerance = 1e-6)
  expect_equal(got$std_error, sqrt(diag(L %*% solve(XtVX) %*% t(L))), tolerance = 1e-6)
})

test_that("two-stage weighting options run and detect heterogeneous variances", {
  trial <- two_stage_trial(c(A = 0.5, B = 0.5, C = 4), seed = 3)
  full <- analyze_augmented_prep_two_stage(trial$data, "y")
  diag_fit <- analyze_augmented_prep_two_stage(trial$data, "y", weighting = "diagonal")
  none <- analyze_augmented_prep_two_stage(trial$data, "y", weighting = "none")
  no_gxl <- analyze_augmented_prep_two_stage(trial$data, "y", genotype_by_location = "none")

  expect_true(full$heterogeneity$heterogeneous)
  expect_gt(full$locations$residual_variance[3], full$locations$residual_variance[1])
  for (f in list(full, diag_fit, none, no_gxl)) {
    expect_equal(nrow(f$adjusted_means), 60)
    expect_false(anyNA(f$adjusted_means$adjusted_mean))
    expect_gt(f$diagnostics$mean_sed, 0)
  }
  expect_null(no_gxl$variance_components)
  expect_match(none$variance_components$component, "residual")
})

test_that("two-stage analysis of one location equals the one-stage analysis", {
  trial <- two_stage_trial()
  one_loc <- trial$data[trial$data$location == "A", ]
  two <- analyze_augmented_prep_two_stage(one_loc, "y")
  one <- analyze_augmented_prep_design(one_loc, "y")
  got <- two$adjusted_means[match(one$adjusted_means$genotype, two$adjusted_means$genotype), ]
  expect_equal(got$adjusted_mean, one$adjusted_means$adjusted_mean, tolerance = 1e-8)
  expect_equal(got$std_error, one$adjusted_means$std_error, tolerance = 1e-8)
})

test_that("spatial stage 1 with SpATS improves estimates under a field trend", {
  skip_if_not_installed("SpATS")
  trial <- two_stage_trial(trend = 4, seed = 5)
  blocks <- analyze_augmented_prep_two_stage(trial$data, "y")
  spatial <- analyze_augmented_prep_two_stage(trial$data, "y", stage1 = "spatial")
  expect_true(all(spatial$locations$stage1_model == "spatial"))
  expect_true(all(spatial$locations$spatial_effective_dim > 0))
  expect_true(all(vapply(spatial$stage1_fits, inherits, logical(1), "SpATS")))
  r <- function(f) stats::cor(f$adjusted_means$adjusted_mean, trial$truth[f$adjusted_means$genotype])
  expect_gt(r(spatial), r(blocks))
})

test_that("spatial stage 1 needs field coordinates", {
  skip_if_not_installed("SpATS")
  design <- augmented_prep_design(sprintf("G%02d", 1:60), c("A", "B"), method = "alpha",
                                  block_size = 6, seed = 1)
  d <- design$design
  d$y <- stats::rnorm(nrow(d))
  expect_error(analyze_augmented_prep_two_stage(d, "y", stage1 = "spatial"), "coordinates")
})

test_that("the p-rep map shows repeated entries, replicates and the location order", {
  design <- augmented_prep_design(paste0("G", 1:42), c("L3", "L1", "L2"), method = "alpha2",
                                  block_size = 4, field_cols = 7, seed = 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_augmented_prep_design(design)
  expect_equal(unique(out$location), c("L3", "L1", "L2"))
  expect_setequal(unique(out$.fill_value), c("Repeated here", "Single plot here"))
  expect_equal(prep_plot_fill_colors(out$.fill_value, "entry_type")[["Repeated here"]], "#E69F00")
  reps <- plot_augmented_prep_design(design, locations = "L1", fill = "replicate", caption = FALSE)
  expect_setequal(unique(reps$.fill_value), c("Replicate 1", "Replicate 2"))
})
