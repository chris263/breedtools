test_that("read_breedbase_file reads csv files", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      location = c("Loc1", "Loc2"),
      germplasmName = c("G1", "G2"),
      phenotype = c(10, 12)
    ),
    path,
    row.names = FALSE
  )

  result <- read_breedbase_file(path)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(names(result), c("location", "germplasmName", "phenotype"))
})

test_that("read_breedbase_file reads augmented row-column design files with phenotype columns", {
  path <- tempfile(fileext = ".design")
  design <- data.frame(
    plots = 1:2,
    row = c(1, 1),
    col = c(1, 2),
    block = c(1, 1),
    rowgroup = c(1, 1),
    colgroup = c(1, 1),
    all_entries = c("Check1", "G1"),
    type = c("check", "entry"),
    rep = c(1, 1),
    is_control = c(1, 0),
    yield = c(10.5, 12.3),
    height = c(40, 45)
  )
  utils::write.table(
    design,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  result <- read_breedbase_file(path)

  expect_s3_class(result, "tbl_df")
  expect_equal(names(result), names(design))
  expect_equal(result$all_entries, c("Check1", "G1"))
  expect_equal(result$yield, c(10.5, 12.3))
})

test_that("runner helpers can use already loaded data frames", {
  df <- data.frame(
    location = c("Loc1", "Loc1", "Loc2", "Loc2"),
    germplasmName = c("G1", "G2", "G1", "G2"),
    rep = c(1, 1, 1, 1),
    phenotype = c(10, 12, 15, 14)
  )

  prepared <- prepare_breedbase_pheno(df, trait_col = "phenotype")

  expect_equal(nrow(prepared), 4)
  expect_equal(
    names(prepared),
    c("environment", "genotype", "replicate", "block", "value")
  )
})

test_that("prepare_breedbase_pheno detects all_entries as genotype column", {
  df <- data.frame(
    location = c("Loc1", "Loc1", "Loc2", "Loc2"),
    all_entries = c("G1", "G2", "G1", "G2"),
    rep = c(1, 1, 1, 1),
    phenotype = c(10, 12, 15, 14)
  )

  prepared <- prepare_breedbase_pheno(df, trait_col = "phenotype")

  expect_equal(as.character(prepared$genotype), c("G1", "G2", "G1", "G2"))
  expect_equal(attr(prepared, "detected_columns")$genotype_col, "all_entries")
})
