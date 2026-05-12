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
    treatments = paste0("G", seq_len(28)),
    controls = c("Check1", "Check2"),
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
    n_candidates = 10,
    seed = 123
  )

  expect_s3_class(design, "augmented_row_column_design")
  expect_equal(nrow(design$design), 36)
  expect_equal(nrow(design$breedbase_design), 36)
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
  expect_equal(sum(design$design$is_control), 8)
  expect_equal(sum(design$design$type == "entry"), 28)
  expect_setequal(design$input$controls, c("Check1", "Check2"))
  expect_setequal(design$input$treatments, paste0("G", seq_len(28)))
})

test_that("augmented_row_column_design prevents duplicate check rows and columns within blocks", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(28)),
    controls = c("Check1", "Check2"),
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
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
  expect_false(any(table(checks$all_entries, checks$row) > 1))
  expect_false(any(table(checks$all_entries, checks$col) > 1))

  by_block <- split(checks, checks$block)
  for (block_checks in by_block) {
    expect_equal(anyDuplicated(block_checks$row), 0)
    expect_equal(anyDuplicated(block_checks$col), 0)
    expect_equal(anyDuplicated(block_checks$all_entries), 0)
  }
})

test_that("plot_augmented_row_column_design draws design objects and data frames", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(28)),
    controls = c("Check1", "Check2"),
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
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
  expect_equal(analysis$anova_formatted$Df, c(2, 8, 2, 8, 4, 111, 8, 143))
  expect_true(all(analysis$anova_formatted$Significance %in% c("", "*", "**")))
  expect_true(analysis$diagnostics$residual_df > 0)
  expect_equal(analysis$diagnostics$anova_error_df, 8)
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

test_that("analyze_augmented_row_column_design captures near-perfect fit warning", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(28)),
    controls = c("Check1", "Check2"),
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$phenotype <- seq_len(nrow(design_data))

  expect_warning(
    analysis <- analyze_augmented_row_column_design(
      design_data,
      response = "phenotype"
    ),
    NA
  )

  expect_true(analysis$diagnostics$saturated_or_near_perfect_fit)
  expect_false(any(is.nan(unlist(analysis$anova[vapply(
    analysis$anova,
    is.numeric,
    logical(1)
  )]))))
  expect_false(any(is.nan(unlist(analysis$anova_formatted[vapply(
    analysis$anova_formatted,
    is.numeric,
    logical(1)
  )]))))
  expect_match(analysis$message, "no positive error degrees of freedom")
})

test_that("augmented ANOVA df use full design layout when phenotype values are missing", {
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
  design_data$yield[1:5] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_equal(unique(analysis$anova_df$k), 12)
  expect_equal(unique(analysis$anova_df$s), 12)
  expect_equal(analysis$anova_formatted$Df, c(2, 8, 2, 8, 4, 111, 8, 143))
})

test_that("augmented ANOVA error df is remaining corrected total df", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(36)),
    controls = paste0("Check", seq_len(3)),
    rows_in_field = 6,
    cols_in_field = 9,
    rows_per_block = 3,
    cols_per_block = 3,
    n_candidates = 5,
    seed = 987
  )
  design_data <- design$design
  design_data$phenotype <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_row_column_design(design_data, response = "phenotype")

  expect_equal(unique(analysis$anova_df$k), 6)
  expect_equal(unique(analysis$anova_df$s), 9)
  expect_equal(unique(analysis$anova_df$gk), 2)
  expect_equal(unique(analysis$anova_df$gs), 3)
  expect_equal(unique(analysis$anova_df$vc), 3)
  expect_equal(unique(analysis$anova_df$ve), 36)
  expect_equal(analysis$anova_formatted$Df, c(1, 3, 2, 5, 2, 38, 2, 53))
})
