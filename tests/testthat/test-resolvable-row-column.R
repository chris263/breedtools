# Designs from Tables 1 and 2 of Williams, John and Whitaker (2006).
wjw_table1 <- list(
  rbind(c(19, 34, 27, 3, 15, 12, 2), c(40, 4, 7, 20, 21, 33, 18), c(23, 35, 17, 29, 22, 28, 38),
        c(10, 14, 31, 39, 32, 5, 42), c(13, 6, 16, 36, 11, 37, 41), c(25, 8, 1, 26, 9, 30, 24)),
  rbind(c(38, 34, 26, 33, 7, 10, 11), c(1, 13, 15, 35, 5, 3, 18), c(32, 20, 17, 41, 8, 24, 19),
        c(12, 30, 6, 39, 23, 21, 16), c(4, 31, 42, 25, 2, 37, 29), c(36, 22, 40, 27, 9, 14, 28)),
  rbind(c(20, 26, 23, 37, 32, 27, 35), c(10, 22, 16, 18, 42, 8, 12), c(28, 4, 15, 39, 7, 13, 24),
        c(1, 41, 14, 34, 30, 29, 40), c(2, 5, 36, 17, 25, 21, 11), c(6, 19, 33, 9, 3, 38, 31))
)
wjw_table2 <- list(
  rbind(c(11, 8, 33, 23, 10, 6, 29), c(39, 28, 17, 4, 2, 21, 38), c(16, 32, 41, 37, 42, 9, 13),
        c(1, 7, 19, 5, 3, 36, 18), c(27, 14, 34, 30, 40, 31, 25), c(26, 24, 12, 35, 20, 22, 15)),
  rbind(c(6, 2, 31, 20, 13, 14, 5), c(42, 36, 26, 19, 35, 25, 38), c(24, 11, 10, 9, 40, 1, 4),
        c(30, 22, 33, 18, 16, 28, 27), c(21, 12, 32, 15, 17, 37, 8), c(3, 23, 34, 39, 7, 29, 41)),
  rbind(c(5, 28, 20, 30, 8, 7, 26), c(2, 19, 14, 12, 39, 31, 37), c(34, 11, 35, 27, 6, 15, 40),
        c(18, 4, 41, 3, 25, 10, 22), c(24, 33, 21, 13, 36, 16, 29), c(17, 42, 23, 1, 32, 38, 9))
)
wjw_df <- function(tab) {
  do.call(rbind, lapply(seq_along(tab), function(q) {
    data.frame(replicate = q, row_in_rep = rep(1:6, 7), col_in_rep = rep(1:7, each = 6),
               all_entries = as.character(as.vector(tab[[q]])), stringsAsFactors = FALSE)
  }))
}

test_that("efficiency factors reproduce Example 2 of the paper", {
  eff <- resolvable_row_column_efficiency(wjw_df(wjw_table1))
  expect_equal(eff$efficiency_spatial, 0.9204, tolerance = 1e-4)
  expect_equal(eff$upper_bound_spatial, 0.9268, tolerance = 1e-4)
  expect_equal(eff$efficiency_rc, 0.6941, tolerance = 1e-4)
  nb <- eff$neighbours
  # Treatment 1: row neighbours 8, 13, 26, 41 and column neighbours 2, 16,
  # 28, 32, 38 (nine in all); treatments 19 and 38 have only eight.
  expect_equal(nb$row_neighbours[nb$treatment == "1"], 4)
  expect_equal(nb$col_neighbours[nb$treatment == "1"], 5)
  expect_equal(nb$neighbours[nb$treatment == "1"], 9)
  expect_equal(nb$neighbours[nb$treatment %in% c("19", "38")], c(8, 8))
  expect_equal(sum(nb$reps_on_edge == 3), 9)
})

test_that("the AR criterion reproduces Table 3 of the paper", {
  d1 <- wjw_df(wjw_table1)
  d2 <- wjw_df(wjw_table2)
  re <- function(p) {
    ar <- list(eta_r = p[1], eta_c = p[2], psi = p[3], rho_r = p[4], rho_c = p[5])
    resolvable_row_column_efficiency(d2, ar = ar)$ar_criterion /
      resolvable_row_column_efficiency(d1, ar = ar)$ar_criterion
  }
  expect_equal(re(c(35, 60, 10, 0.5, 0.5)), 0.9567, tolerance = 1e-4)
  expect_equal(re(c(70, 60, 1, 0.9, 0.9)), 1.0434, tolerance = 1e-4)
  expect_equal(re(c(0, 0, 0.5, 0.8, 0.8)), 1.0054, tolerance = 1e-4)
})

test_that("the generator builds a resolvable row-column design", {
  design <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5, n_starts = 1, n_iter = 600, seed = 1)
  d <- design$design
  expect_s3_class(design, "resolvable_row_column_design")
  expect_equal(nrow(d), 60)
  expect_true(all(table(d$all_entries, d$replicate) == 1))
  expect_equal(anyDuplicated(paste(d$row, d$col)), 0)
  s <- design$summary
  expect_true(all(c(s$efficiency_rc, s$efficiency_rows, s$efficiency_cols, s$efficiency_spatial) > 0.5))
  expect_lte(s$efficiency_spatial, s$upper_bound_spatial + 1e-8)
  expect_true(is.na(s$latin_ok))
  expect_equal(names(design$breedbase_design),
               c("plots", "block", "all_entries", "rep", "is_control", "row_number", "col_number"))
  # The search improves on a random design.
  random <- design
  set.seed(9)
  random$design$all_entries <- unlist(lapply(split(d$all_entries, d$replicate), sample))
  expect_gt(s$efficiency_rc, resolvable_row_column_efficiency(random)$efficiency_rc)
})

test_that("latinized and t-latinized designs keep treatments out of repeated long lines", {
  latin <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5, latinize = TRUE,
                                        n_starts = 1, n_iter = 400, seed = 2)
  d <- latin$design
  expect_true(latin$summary$latin_ok)
  expect_true(all(table(d$all_entries, d$long_col) <= 1))
  expect_equal(max(d$row), 12)

  tlatin <- resolvable_row_column_design(paste0("T", 1:20), 4, 4, 5, latinize = TRUE, latin_t = 2,
                                         replicate_layout = "horizontal", n_starts = 1, n_iter = 400, seed = 3)
  dt <- tlatin$design
  expect_true(tlatin$summary$latin_ok)
  expect_equal(max(dt$col), 20)
  grp <- (dt$replicate - 1) %/% 2
  expect_false(any(duplicated(paste(grp, dt$all_entries, dt$row_in_rep))))
  expect_error(resolvable_row_column_design(paste0("T", 1:20), 6, 4, 5, latinize = TRUE), "latin")
})

test_that("fewer treatments than plots gives filler plots", {
  design <- resolvable_row_column_design(paste0("P", 1:19), 2, 4, 5, n_starts = 1, n_iter = 300, seed = 4)
  expect_equal(design$summary$n_fillers, 1)
  expect_equal(sum(design$design$is_filler), 2)
  expect_true(all(table(design$design$all_entries[design$design$is_filler == 0]) == 2))
})

test_that("a start design can be laid out and evaluated", {
  d1 <- wjw_df(wjw_table1)
  d1$all_entries <- paste0("T", d1$all_entries)
  design <- resolvable_row_column_design(paste0("T", 1:42), 3, 6, 7, start_design = d1,
                                         n_iter = 1, spatial = FALSE)
  expect_equal(design$summary$efficiency_rc, 0.6941, tolerance = 1e-4)
  expect_equal(design$summary$efficiency_spatial, 0.9204, tolerance = 1e-4)
  expect_error(
    resolvable_row_column_design(paste0("T", 1:42), 3, 6, 7, start_design = d1[-1, ], n_iter = 1),
    "does not match"
  )
})

rrc_trial <- function() {
  design <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5, latinize = TRUE,
                                         n_starts = 1, n_iter = 300, seed = 5)
  d <- design$design
  set.seed(6)
  g <- stats::setNames(stats::rnorm(20, sd = 2), paste0("T", 1:20))
  d$y <- 50 + g[d$all_entries] + c(0, 2, -1)[d$replicate] +
    stats::rnorm(12)[(d$replicate - 1) * 4 + d$row_in_rep] +
    stats::rnorm(15, sd = 0.6)[(d$replicate - 1) * 5 + d$col_in_rep] + stats::rnorm(nrow(d))
  list(design = design, data = d, truth = g)
}

test_that("the RC model matches lme4", {
  skip_if_not_installed("lme4")
  trial <- rrc_trial()
  d <- trial$data
  fit <- analyze_resolvable_row_column(d, "y", model = "rc", compare = FALSE)
  expect_s3_class(fit, "resolvable_row_column_analysis")
  expect_output(print(fit), "Resolvable row-column analysis")
  ref <- lme4::lmer(y ~ factor(replicate) + all_entries + (1 | replicate:row_in_rep) +
                      (1 | replicate:col_in_rep), data = d)
  vc <- as.data.frame(lme4::VarCorr(ref))
  p <- stats::setNames(fit$variance_parameters$estimate, fit$variance_parameters$parameter)
  expect_equal(p[["row_variance"]], vc$vcov[vc$grp == "replicate:row_in_rep"], tolerance = 1e-3)
  expect_equal(p[["column_variance"]], vc$vcov[vc$grp == "replicate:col_in_rep"], tolerance = 1e-3)
  expect_equal(p[["sigma2"]], vc$vcov[vc$grp == "Residual"], tolerance = 1e-3)
  fe <- lme4::fixef(ref)
  lv <- sort(unique(d$all_entries))
  got <- fit$adjusted_means$adjusted_mean[match(lv, fit$adjusted_means$genotype)]
  expect_equal(got - got[1], unname(c(0, fe[paste0("all_entries", lv[-1])])), tolerance = 1e-4)
})

test_that("LV + RC and AR + RC models are fitted and compared", {
  trial <- rrc_trial()
  d <- trial$data
  fit <- analyze_resolvable_row_column(d, "y", model = "lv")
  cmp <- fit$model_comparison
  expect_equal(cmp$model[1], "RC")
  expect_equal(nrow(cmp), 5)
  expect_true(all(is.finite(cmp$deviance)))
  # Nested LV models cannot fit worse than RC.
  expect_true(all(cmp$deviance[2:4] <= cmp$deviance[1] + 1e-4))
  expect_true(all(c("phi_r", "phi_c") %in% fit$variance_parameters$parameter))

  eq <- analyze_resolvable_row_column(d, "y", model = "lv", equal_phi = TRUE, compare = FALSE)
  expect_true("phi" %in% eq$variance_parameters$parameter)
  ar <- analyze_resolvable_row_column(d, "y", model = "ar", compare = FALSE)
  expect_true(all(c("psi", "rho_r", "rho_c") %in% ar$variance_parameters$parameter))
  long <- analyze_resolvable_row_column(d, "y", model = "rc", long_col = "long_col", compare = FALSE)
  expect_false(anyNA(long$adjusted_means$adjusted_mean))
  expect_true("Rows within replicates" %in% fit$anova$term)
  r <- stats::cor(fit$adjusted_means$adjusted_mean, trial$truth[fit$adjusted_means$genotype])
  expect_gt(r, 0.8)
})

test_that("resolvable row-column designs can be plotted", {
  design <- resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5, n_starts = 1, n_iter = 200, seed = 7)
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  out <- plot_resolvable_row_column_design(design, label = "all_entries", fill = "replicate", legend = FALSE)
  grDevices::dev.off()
  expect_s3_class(out, "tbl_df")
})
