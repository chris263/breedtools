test_that("guess_column detects column names", {
  df <- data.frame(
    germplasmName = c("G1", "G2"),
    phenotype = c(10, 20)
  )

  result <- guess_column(df, c("germplasmName", "accession_name"))

  expect_equal(result, "germplasmName")
})
