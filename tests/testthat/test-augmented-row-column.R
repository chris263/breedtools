test_that("augmented_row_column_capacity calculates plot counts", {
  capacity <- augmented_row_column_capacity(
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
    n_controls = 2
  )

  expect_equal(capacity$n_blocks, 4)
  expect_equal(capacity$total_plots, 36)
  expect_equal(capacity$check_plots, 8)
  expect_equal(capacity$entry_plots, 28)
})

test_that("augmented_row_column_design returns complete and breedbase outputs", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 10,
    seed = 123
  )

  expect_s3_class(design, "augmented_row_column_design")
  expect_equal(nrow(design$design), 144)
  expect_equal(nrow(design$breedbase_design), 144)
  expect_equal(
    names(design$design),
    c(
      "plots",
      "row",
      "col",
      "block",
      "rowgroup",
      "colgroup",
      "all_entries",
      "type",
      "rep",
      "is_control"
    )
  )
  expect_equal(
    names(design$breedbase_design),
    c("plots", "block", "all_entries", "rep", "is_control")
  )
  expect_equal(sum(design$design$is_control), 36)
  expect_equal(sum(design$design$type == "entry"), 108)
  expect_setequal(design$input$controls, paste0("Check", seq_len(4)))
  expect_setequal(design$input$treatments, paste0("G", seq_len(108)))
  expect_setequal(design$design$all_entries[design$design$is_control == 0], paste0("G", seq_len(108)))
  expect_equal(design$summary$error_df, 6)
})

test_that("augmented_row_column_design prevents duplicate check rows and columns within blocks", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 10,
    seed = 321
  )

  checks <- design$design[design$design$is_control == 1, ]

  by_block <- split(checks, checks$block)
  for (block_checks in by_block) {
    expect_equal(anyDuplicated(block_checks$row), 0)
    expect_equal(anyDuplicated(block_checks$col), 0)
    expect_equal(anyDuplicated(block_checks$all_entries), 0)
  }
})

test_that("augmented_row_column_design handles large layouts without recursive search", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(576)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 36,
    cols_in_field = 18,
    rows_per_block = 6,
    cols_per_block = 6,
    n_candidates = 3,
    seed = 123
  )

  checks <- design$design[design$design$is_control == 1, ]

  expect_equal(nrow(design$design), 648)
  expect_equal(nrow(checks), 72)
  expect_setequal(checks$row, seq_len(36))
  expect_setequal(checks$col, seq_len(18))
  expect_false(any(table(checks$all_entries, checks$row) > 1))
  expect_false(any(table(checks$all_entries, checks$col) > 1))

  by_block <- split(checks, checks$block)
  for (block_checks in by_block) {
    expect_equal(anyDuplicated(block_checks$row), 0)
    expect_equal(anyDuplicated(block_checks$col), 0)
    expect_equal(anyDuplicated(block_checks$all_entries), 0)
  }
})


test_that("augmented_row_column_design randomizes check positions across seeds", {
  check_signature <- function(seed) {
    design <- augmented_row_column_design(
      treatments = paste0("G", seq_len(108)),
      controls = paste0("Check", seq_len(4)),
      rows_in_field = 12,
      cols_in_field = 12,
      rows_per_block = 4,
      cols_per_block = 4,
      n_candidates = 1,
      seed = seed
    )
    checks <- design$design[design$design$is_control == 1, ]
    paste(sort(paste(checks$all_entries, checks$row, checks$col, sep = ":")), collapse = "|")
  }

  signatures <- vapply(c(101, 202, 303), check_signature, character(1))

  expect_gt(length(unique(signatures)), 1)
})

test_that("plot_augmented_row_column_design draws design objects and data frames", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 111
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plotted_object <- plot_augmented_row_column_design(
    design,
    label = "none",
    fill = "block",
    legend = FALSE
  )
  plotted_data <- plot_augmented_row_column_design(
    design$design,
    label = "all_entries",
    fill = "type",
    legend = FALSE
  )
  grDevices::dev.off()

  expect_s3_class(plotted_object, "tbl_df")
  expect_s3_class(plotted_data, "tbl_df")
  expect_true(all(c(".plot_row", ".plot_col", ".fill_value") %in% names(plotted_object)))
  expect_true(file.exists(plot_file))
})

test_that("augmented_row_column_design validates treatment count", {
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(27)),
      controls = c("Check1", "Check2"),
      rows_in_field = 6,
      cols_in_field = 6,
      rows_per_block = 3,
      cols_per_block = 3,
      n_candidates = 5
    ),
    "Expected 28 treatments"
  )
})

test_that("analyze_augmented_row_column_design returns model outputs", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$yield <- 100 +
    as.numeric(factor(design_data$rowgroup)) * 2 +
    as.numeric(factor(design_data$colgroup)) * 3 +
    ifelse(design_data$is_control == 1, 5, 0) +
    stats::rnorm(nrow(design_data), sd = 1)

  analysis <- analyze_augmented_row_column_design(
    design_data,
    response = "yield"
  )

  expect_s3_class(analysis, "augmented_row_column_analysis")
  expect_s3_class(analysis$model, "lm")
  expect_false(any(is.na(stats::coef(analysis$model))))
  expect_equal(nrow(analysis$plot_results), nrow(design_data))
  expect_equal(nrow(analysis$genotype_summary), dplyr::n_distinct(design_data$all_entries))
  expect_true(all(c("rowgroup:colgroup", "rowgroup:row", "colgroup:col", "genotype") %in% analysis$anova$term))
  expect_equal(
    analysis$anova_formatted$Source,
    c(
      "Row groups",
      "Rows, nested within row groups",
      "Column groups",
      "Columns, nested within column groups",
      "Row groups x column groups (blocks)",
      "Genotypes",
      "Error",
      "Corrected total"
    )
  )
  expect_equal(analysis$anova_formatted$Df, c(2, 9, 2, 9, 4, 111, 6, 143))
  expect_true(all(analysis$anova_formatted$Significance %in% c("", "*", "**")))
  expect_true(analysis$diagnostics$residual_df > 0)
  expect_equal(analysis$diagnostics$anova_error_df, 6)
})

test_that("augmented row-column analysis print method shows formatted ANOVA", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 456
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_output(print(analysis), "Augmented row-column design analysis")
  expect_output(print(analysis), "Corrected total")
})

test_that("analyze_augmented_row_column_design validates required columns", {
  expect_error(
    analyze_augmented_row_column_design(
      data.frame(rowgroup = 1, colgroup = 1, row = 1, col = 1, all_entries = "G1"),
      response = "yield"
    ),
    "Missing required column"
  )
})

test_that("analyze_augmented_row_column_design flags a perfect fit", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$phenotype <- design_data$row * 2 + design_data$col

  expect_warning(
    analysis <- analyze_augmented_row_column_design(design_data, response = "phenotype"),
    NA
  )

  expect_true(analysis$diagnostics$saturated_or_near_perfect_fit)
  expect_false(any(is.nan(unlist(analysis$anova[vapply(analysis$anova, is.numeric, logical(1))]))))
  expect_false(any(is.nan(unlist(analysis$anova_formatted[vapply(
    analysis$anova_formatted,
    is.numeric,
    logical(1)
  )]))))
  expect_match(analysis$message, "essentially perfect fit")
})

test_that("analyze_augmented_row_column_design reports when no error df remain", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  # Keep only the checks of one block, so no genotype is replicated.
  design_data$yield[design_data$is_control == 1 & design_data$block != 1] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_equal(analysis$diagnostics$anova_error_df, 0)
  expect_match(analysis$message, "no positive error degrees of freedom")
})

test_that("augmented ANOVA df come from the observed plots when values are missing", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 789
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  entry_plots <- which(design_data$is_control == 0)[1:5]
  design_data$yield[entry_plots] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_equal(unique(analysis$anova_df$k), 12)
  expect_equal(unique(analysis$anova_df$s), 12)
  expect_equal(analysis$diagnostics$design_error_df, 6)
  # Missing entries remove genotype df only; error df is unchanged.
  expect_equal(analysis$anova_formatted$Df, c(2, 9, 2, 9, 4, 106, 6, 138))
})

test_that("augmented ANOVA df follow the design formulas", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(320)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 24,
    cols_in_field = 16,
    rows_per_block = 6,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 987
  )
  design_data <- design$design
  design_data$phenotype <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_row_column_design(design_data, response = "phenotype")

  expect_equal(unique(analysis$anova_df$k), 24)
  expect_equal(unique(analysis$anova_df$s), 16)
  expect_equal(unique(analysis$anova_df$gk), 4)
  expect_equal(unique(analysis$anova_df$gs), 4)
  expect_equal(unique(analysis$anova_df$vc), 4)
  expect_equal(unique(analysis$anova_df$ve), 320)
  # rows within groups = k - gk = 20; columns within groups = s - gs = 12
  expected_df <- c(3, 20, 3, 12, 9, 323, 13, 383)
  expect_equal(analysis$anova_formatted$Df, expected_df)
  expect_equal(analysis$anova_df$display_df, expected_df[1:7])
})

test_that("augmented ANOVA matches a full lm() fit", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 42
  )
  d <- design$design
  d$y <- stats::rnorm(nrow(d)) + d$row * 0.3
  d$y[c(3, 50, 77)] <- NA_real_

  analysis <- analyze_augmented_row_column_design(d, response = "y")

  ref_data <- data.frame(
    y = d$y,
    rg = factor(d$rowgroup),
    rn = interaction(d$rowgroup, d$row, drop = TRUE),
    cg = factor(d$colgroup),
    cn = interaction(d$colgroup, d$col, drop = TRUE),
    bk = interaction(d$rowgroup, d$colgroup, drop = TRUE),
    g = factor(d$all_entries)
  )
  ref <- stats::anova(stats::lm(y ~ rg + rn + cg + cn + bk + g, data = ref_data))

  expect_equal(analysis$anova$df, ref$Df)
  expect_equal(analysis$anova$sum_sq, ref$`Sum Sq`, tolerance = 1e-8)
  expect_equal(analysis$anova$p_value[1:6], ref$`Pr(>F)`[1:6], tolerance = 1e-6)
})

test_that("augmented_row_column_design rejects layouts that cannot be analyzed", {
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(28)),
      controls = c("Check1", "Check2"),
      rows_in_field = 6,
      cols_in_field = 6,
      rows_per_block = 3,
      cols_per_block = 3,
      n_candidates = 5
    ),
    "no error degrees of freedom"
  )

  capacity <- augmented_row_column_capacity(10, 10, 5, 5, 2)
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(capacity$entry_plots)),
      controls = c("Check1", "Check2"),
      rows_in_field = 10,
      cols_in_field = 10,
      rows_per_block = 5,
      cols_per_block = 5,
      n_candidates = 5
    ),
    "Not every field row can contain a check"
  )
})

test_that("augmented_row_column_capacity reports error df", {
  expect_equal(augmented_row_column_capacity(12, 12, 4, 4, 4)$error_df, 6)
  expect_equal(augmented_row_column_capacity(6, 6, 3, 3, 2)$error_df, -5)
})

test_that("augmented row-column adjusted means match a full lm fit", {
  design <- augmented_row_column_design(
    treatments = paste0("E", seq_len(108)),
    controls = paste0("C", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 9
  )
  set.seed(3)
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data), 50, 5) + 0.5 * design_data$row
  design_data$yield[c(3, 40, 77)] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  observed <- design_data[!is.na(design_data$yield), ]
  observed$genotype <- factor(observed$all_entries)
  observed$rowgroup <- factor(observed$rowgroup)
  observed$colgroup <- factor(observed$colgroup)
  observed$row_nested <- interaction(observed$rowgroup, observed$row, drop = TRUE)
  observed$col_nested <- interaction(observed$colgroup, observed$col, drop = TRUE)
  fit <- stats::lm(
    yield ~ rowgroup + row_nested + colgroup + col_nested + rowgroup:colgroup + genotype,
    data = observed
  )
  ref <- reference_adjusted_means(fit, observed)
  got <- analysis$adjusted_means[match(ref$genotype, analysis$adjusted_means$genotype), ]

  expect_equal(got$adjusted_mean, ref$mean, tolerance = 1e-8)
  expect_equal(got$std_error, ref$se, tolerance = 1e-8)
  expect_equal(nrow(analysis$adjusted_means), 112)
  expect_true(all(is.na(analysis$adjusted_means$adjusted_mean[analysis$adjusted_means$n_observations == 0])))
  expect_false(is.unsorted(-stats::na.omit(analysis$adjusted_means$adjusted_mean)))
  expect_true(all(c("adjusted_mean", "std_error") %in% names(analysis$genotype_summary)))
  expect_gt(analysis$diagnostics$mean_sed, 0)
  expect_equal(
    analysis$diagnostics$lsd_5pct,
    stats::qt(0.975, analysis$diagnostics$anova_error_df) * analysis$diagnostics$mean_sed
  )
  expect_output(print(analysis), "Top adjusted means")
})

test_that("augmented row-column adjusted means are NA when not estimable", {
  design <- augmented_row_column_design(
    treatments = paste0("E", seq_len(108)),
    controls = paste0("C", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 9
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  design_data$yield[design_data$is_control == 1 & design_data$row == 1] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_true(analysis$diagnostics$rank_deficient)
  expect_true(all(is.na(analysis$adjusted_means$adjusted_mean)))
  expect_match(analysis$message, "not estimable")
})
