local_test_dir <- function(env = parent.frame()) {
  path <- tempfile("breedtools-test-")
  dir.create(path)
  do.call(on.exit, list(bquote(unlink(.(path), recursive = TRUE)), add = TRUE), envir = env)
  path
}

example_pheno <- function() {
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
}

test_that("bundled example data is available", {
  expect_true(file.exists(example_pheno()))
  expect_true(file.exists(
    system.file("extdata", "breedbase_locations.csv", package = "breedtools")
  ))
  expect_true(file.exists(
    system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools")
  ))
})

test_that("run_breedbase_variance_check writes nothing by default", {
  old_wd <- setwd(local_test_dir())
  on.exit(setwd(old_wd), add = TRUE)

  result <- run_breedbase_variance_check(
    example_pheno(),
    trait_col = "fresh root yield|CO_334:0000013",
    n_sim = 500,
    export = FALSE
  )

  expect_equal(nrow(result$env_anova), 4)
  expect_length(list.files("."), 0)

  expect_error(
    run_breedbase_variance_check(
      example_pheno(),
      trait_col = "fresh root yield|CO_334:0000013",
      n_sim = 500,
      export = TRUE
    ),
    "output_dir"
  )
  expect_length(list.files("."), 0)
})

test_that("run_breedbase_variance_check exports to output_dir", {
  out_dir <- file.path(local_test_dir(), "out")

  run_breedbase_variance_check(
    example_pheno(),
    trait_col = "fresh root yield|CO_334:0000013",
    n_sim = 500,
    export = TRUE,
    output_dir = out_dir
  )

  expect_length(list.files(out_dir, pattern = "\\.csv$"), 3)
})

test_that("run_environment_stratification writes JSON only to output_dir", {
  expect_error(
    run_environment_stratification(
      example_pheno(),
      trait = "fresh root yield",
      export_json = TRUE
    ),
    "output_dir"
  )

  out_dir <- file.path(local_test_dir(), "json")
  result <- run_environment_stratification(
    example_pheno(),
    trait = "fresh root yield",
    export_json = TRUE,
    output_dir = out_dir
  )

  expect_true(file.exists(file.path(out_dir, "environment_stratification.json")))
  expect_equal(nrow(result$group_summary), 2)
})

test_that("read_locations accepts run_environment_stratification results", {
  result <- run_environment_stratification(example_pheno(), trait = "fresh root yield")
  locations <- read_locations(
    system.file("extdata", "breedbase_locations.csv", package = "breedtools"),
    phenotype_data = result
  )

  expect_equal(nrow(locations), 4)
  expect_s3_class(map_environment_groups(result, locations), "leaflet")
})
