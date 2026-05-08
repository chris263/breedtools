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
