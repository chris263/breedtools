test_that("entry counts alone produce bounded options for each native design family", {
  options <- design_options(100, n_locations = 3)
  expect_setequal(options$design, c("augmented_row_column", "prep", "resolvable"))
  expect_true(all(table(options$design) <= 6))
  expect_equal(options$option, seq_len(nrow(options)))
  expect_true(all(options$n_entries == 100L))
  expect_equal(options$field_plots, options$rows_in_field * options$cols_in_field)
  expect_equal(options$total_plots, 3 * options$plots_per_location)
  expect_equal(options$unused_plots, options$field_plots - options$plots_per_location)
  expect_true(all(options$field_plots <= 424))
  expect_true(all(options$error_df >= 1))
  expect_true(all(pmax(options$rows_in_field / options$cols_in_field,
                      options$cols_in_field / options$rows_in_field) <= 3))
  expect_setequal(options$n_checks[options$design == "augmented_row_column"], 2:6)
  expect_setequal(options$n_reps[options$design == "resolvable"], 2:3)
  expect_equal(nrow(attr(options, "diagnostics")), 3)
})

test_that("augmented row-column options have exactly the requested entry positions, including prime counts", {
  options <- design_options(108, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12)
  expect_gt(nrow(options), 0)
  expect_identical(options$method[1], "latin")
  expect_equal(options$rows_per_block[1], 4)
  expect_equal(options$cols_per_block[1], 4)
  for (n in c(100, 101)) {
    options <- design_options(n, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12, n_options = 2)
    expect_gt(nrow(options), 0)
    for (i in seq_len(nrow(options))) {
      args <- options$arguments[[i]]
      cap <- do.call(augmented_row_column_capacity,
        c(args[intersect(names(args), c("rows_in_field", "cols_in_field", "rows_per_block",
                                      "cols_per_block", "checks_per_block"))], list(n_controls = 4)))
      expect_equal(cap$entry_plots, n)
      expect_equal(cap$check_plots, 144 - n)
      expect_equal(cap$error_df, options$error_df[i])
      expect_gte(cap$min_rowgroup_slack, 0)
      expect_gte(cap$min_colgroup_slack, 0)
    }
    selected <- options[1, ]
    # The explicit constructor needs every field/block setting, including
    # extra check replication for entry counts that do not fill regular blocks.
    generated <- augmented_row_column_design(
      treatments = paste0("G", seq_len(n)), controls = paste0("C", 1:4),
      rows_in_field = selected$rows_in_field, cols_in_field = selected$cols_in_field,
      rows_per_block = selected$rows_per_block, cols_per_block = selected$cols_per_block,
      checks_per_block = selected$checks_per_block[[1]], method = selected$method,
      n_candidates = 5, seed = 7
    )
    expect_equal(sum(generated$design$is_control == 0), n)
    expect_equal(sum(generated$design$is_control == 1), 144 - n)
    expect_equal(anyDuplicated(generated$design$all_entries[generated$design$is_control == 0]), 0)
  }
})

test_that("p-rep options obey repeated-entry budgets and map to generator arguments", {
  options <- design_options(100, "prep", n_locations = 3, cols_in_field = 10)
  expect_gt(nrow(options), 0)
  expect_true(all(options$n_checks == 0))
  expect_true(all(options$n_repeated * 3 <= 100))
  expect_equal(options$plots_per_location, 100 + options$n_repeated)
  expect_equal(options$error_df, options$n_repeated - options$n_blocks + 1)
  expect_equal(options$rows_in_field, ceiling(options$plots_per_location / 10))
  args <- options$arguments[[1]]
  args$treatments <- paste0("G", 1:100)
  args$locations <- c("North", "Central", "South")
  args$n_candidates <- 5
  args$n_efficiency_candidates <- 1
  args$seed <- 1
  generated <- do.call(options$generator[1], args)
  expect_equal(nrow(generated$design), options$total_plots[1])
  expect_equal(generated$summary$n_disconnected_locations, 0)
  expect_equal(max(generated$design$row), options$rows_in_field[1])
  expect_equal(max(generated$design$col), 10)
  expect_message(empty <- design_options(100, "prep", n_checks = 4), "No design options")
  expect_equal(nrow(empty), 0)
  expect_match(attr(empty, "diagnostics")$reason, "no separate checks")
  # The default can adapt duplication to more locations than 15% would allow.
  many_sites <- design_options(200, "prep", n_locations = 8)
  expect_gt(nrow(many_sites), 0)
  expect_true(all(many_sites$n_repeated * 8 <= 200))
})

test_that("augmented row-column blocks have at least one row and column per check", {
  # Check both orientations: a narrow row group or column group is excluded.
  for (dims in list(c(9, 14), c(14, 9), c(12, 12))) {
    options <- design_options(100, "augmented_row_column", n_checks = 4,
                               rows_in_field = dims[1], cols_in_field = dims[2], n_options = 20)
    expect_gt(nrow(options), 0)
    expect_true(all(options$rows_per_block >= 4L))
    expect_true(all(options$cols_per_block >= 4L))
  }
  # The minimum follows each option's own check count.
  options <- design_options(100, "augmented_row_column", n_checks = 2:6)
  expect_setequal(options$n_checks, 2:6)
  expect_true(all(options$rows_per_block >= options$n_checks))
  expect_true(all(options$cols_per_block >= options$n_checks))
  for (dims in list(c(3, 20), c(20, 3))) {
    expect_message(empty <- design_options(20, "augmented_row_column", n_checks = 4,
      rows_in_field = dims[1], cols_in_field = dims[2]), "No design options")
    expect_equal(nrow(empty), 0)
    expect_match(attr(empty, "diagnostics")$reason, "n_checks rows and n_checks columns")
  }
})

test_that("resolvable options account for all checks, replicates and filler positions", {
  options <- design_options(42, "resolvable", n_reps = 3, cols_in_field = 7)
  expect_equal(options$rows_in_field[1], 18)
  expect_equal(options$error_df[1], 49)
  expect_equal(options$unused_plots[1], 0)
  # A larger fixed field needs fillers, rather than extra test entries.
  options <- design_options(40, "resolvable", n_checks = 2, n_reps = 3, rows_in_field = 21, cols_in_field = 7)
  expect_gt(nrow(options), 0)
  expect_equal(options$plots_per_location, 126)
  expect_equal(options$unused_plots, 21)
  expect_equal(options$check_plots, 6)
  args <- options$arguments[[1]]
  args$treatments <- c(paste0("G", 1:40), "CHECK_A", "CHECK_B")
  args$n_iter <- 1
  args$n_starts <- 1
  args$seed <- 1
  generated <- do.call(options$generator[1], args)
  expect_equal(sum(generated$design$is_filler), 21)
  expect_equal(sum(generated$design$is_filler == 0), 126)
  expect_equal(max(generated$design$row), 21)
  expect_equal(max(generated$design$col), 7)
  more_reps <- design_options(17, "resolvable", n_checks = 0, n_reps = 5)
  expect_gt(nrow(more_reps), 0)
  expect_true(all(more_reps$n_reps == 5))
})

test_that("field constraints and empty results keep a consistent schema", {
  options <- design_options(100, rows_in_field = c(12, 18, 24), cols_in_field = c(8, 10), max_plots = 240)
  expect_true(all(options$rows_in_field %in% c(12, 18, 24)))
  expect_true(all(options$cols_in_field %in% c(8, 10)))
  expect_true(all(options$field_plots <= 240))
  expect_message(empty <- design_options(100, rows_in_field = 4, cols_in_field = 4), "No design options")
  expect_equal(nrow(empty), 0)
  expect_identical(names(empty), names(options))
  expect_true(is.list(empty$arguments))
  expect_true(all(attr(empty, "diagnostics")$options == 0))
  expect_message(design_options(100, max_plots = 50), "No design options")
  no_checks <- design_options(100, n_checks = 0)
  expect_false("augmented_row_column" %in% no_checks$design)
  expect_match(attr(no_checks, "diagnostics")$reason[1], "require at least one")
})

test_that("planning is deterministic and does not change random state", {
  withr::local_preserve_seed()
  set.seed(33)
  before <- .Random.seed
  first <- design_options(100)
  expect_identical(.Random.seed, before)
  expect_identical(design_options(100), first)
  expect_true(all(design_options(100, min_error_df = 10)$error_df >= 10))
})

test_that("invalid counts and design requests have explicit errors", {
  for (n in list(0, -1, 2.5, NA_real_, Inf, "100", c(100, 101))) {
    expect_error(design_options(n), "n_entries")
  }
  expect_error(design_options(100, "unknown"), "arg")
  expect_error(design_options(100, "augmented"), "augmented_row_column")
  expect_error(design_options(100, n_checks = -1), "n_checks")
  expect_error(design_options(100, rows_in_field = 0), "rows_in_field")
  expect_error(design_options(100, cols_in_field = 2.5), "cols_in_field")
  expect_error(design_options(100, n_locations = 0), "n_locations")
  expect_error(design_options(100, n_reps = 1), "n_reps")
  expect_error(design_options(100, prep_prop = NA_real_), "prep_prop")
  expect_error(design_options(100, prep_prop = 1.1), "prep_prop")
  expect_error(design_options(100, n_options = 0), "n_options")
  expect_error(design_options(100, n_candidates = 0), "n_candidates")
})

test_that("augmented row-column options expose all generator settings with the correct name", {
  options <- design_options(108, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12,
                             n_candidates = 37)
  expect_true(all(options$design == "augmented_row_column"))
  expect_true(all(c("rows_in_field", "cols_in_field", "rows_per_block", "cols_per_block",
                    "n_candidates", "checks_per_block", "method") %in% names(options)))
  expect_true(all(options$n_candidates == 37L))
  expect_true(all(vapply(options$arguments, function(a) identical(a$n_candidates, 37L), logical(1))))
  expect_identical(options$rows_in_field[1], 12L)
  expect_identical(options$cols_in_field[1], 12L)
  expect_null(options$checks_per_block[[1]])
})

test_that("n_checks x n_checks blocks are the first option when possible, not a requirement", {
  options <- design_options(36, "augmented_row_column", n_checks = 3, rows_in_field = 9)
  expect_equal(options$rows_per_block[1], 3)
  expect_equal(options$cols_per_block[1], 3)
  expect_true(options$square_check_blocks[1])
  expect_false(any(options$square_check_blocks[-1]))
  # Other block shapes are still offered.
  expect_true(any(options$rows_per_block != 3 | options$cols_per_block != 3))
  # With 4 checks, 12 x 12 has 4 x 4 blocks; 108 entries give the Latin layout.
  latin <- design_options(108, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12)
  expect_true(latin$square_check_blocks[1])
  # When no such layout exists, the usual order is kept.
  none <- design_options(100, "augmented_row_column", n_checks = 4, rows_in_field = 9, cols_in_field = 14)
  expect_false(any(none$square_check_blocks))
})

test_that("rows and cols still work as deprecated names of the field dimensions", {
  expect_warning(old <- design_options(108, "augmented_row_column", n_checks = 4, rows = 12, cols = 12),
                 "deprecated")
  new <- design_options(108, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12)
  expect_identical(old, new)
  expect_error(suppressWarnings(design_options(100, rows = 9, rows_in_field = 9)), "rows_in_field")
})
