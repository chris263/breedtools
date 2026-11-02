test_that("hartley_fmax_test works", {
  env_anova <- tibble::tibble(
    environment = c("Loc1", "Loc2", "Loc3"),
    ms_error = c(10, 12, 15),
    df_error = c(20, 20, 20)
  )

  result <- hartley_fmax_test(
    env_anova,
    alpha = 0.05,
    n_sim = 1000,
    seed = 123
  )

  expect_true("fmax_observed" %in% names(result))
  expect_true("homogeneous" %in% names(result))
  expect_equal(result$k, 3)
})

test_that("hartley_critical_sim does not change the global RNG state", {
  set.seed(42)
  before <- .Random.seed
  hartley_critical_sim(c(10, 12), n_sim = 100, seed = 123)
  expect_identical(.Random.seed, before)

  # Reproducible with the same seed
  expect_identical(
    hartley_critical_sim(c(10, 12), n_sim = 100, seed = 1),
    hartley_critical_sim(c(10, 12), n_sim = 100, seed = 1)
  )
})

test_that("hartley_critical_sim removes .Random.seed if none existed", {
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv))

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  hartley_critical_sim(c(10, 12), n_sim = 100, seed = 123)
  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})
