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
  expect_equal(nrow(analysis$plot_results), nrow(design_data))
  expect_equal(nrow(analysis$genotype_summary), dplyr::n_distinct(design_data$all_entries))
  expect_true(all(c("rowgroup:colgroup", "rowgroup:row", "colgroup:col", "genotype") %in% analysis$anova$term))
  expect_true(analysis$diagnostics$residual_df > 0)
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
