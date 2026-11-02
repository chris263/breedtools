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
  design <- augmented_prep_design(
    treatments = paste0("G", seq_len(10)),
    locations = c("Loc1", "Loc2"),
    block_size = c(Loc1 = 4, Loc2 = 5),
    n_blocks_per_location = c(Loc1 = 4, Loc2 = 4),
    n_duplicate_per_location = c(Loc1 = 2, Loc2 = 3),
    replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3),
    n_candidates = 5,
    seed = 404,
    field_cols = c(Loc1 = 6, Loc2 = 8)
  )

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
  design <- augmented_prep_design(
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
  )

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
