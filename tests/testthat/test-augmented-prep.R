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
