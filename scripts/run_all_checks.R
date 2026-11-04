#!/usr/bin/env Rscript
# =============================================================================
# breedtools: automated checks of every design and analysis (for CI)
# For a step-by-step walkthrough, see scripts/validate_designs.R.
#
# Each section generates designs, analyses simulated trials, prints the key
# results (next to published values where the methods come from a paper),
# and records pass/fail checks against independent references (lm(), lme4,
# brute-force matrix algebra, or the published tables).
#
# Sections (`--only` keys in brackets):
#   1. [rc]       Augmented row-column designs, Latin construction
#                 (Federer, Nair and Raghavarao, 1975): layout sweep,
#                 structure, ANOVA and adjusted means vs lm()
#   2. [rcflex]   Flexible augmented row-column designs (Piepho and Williams,
#                 2016): Table 2, Examples 1-3, unequal groups, incomplete
#                 blocks, extra check replication, D-optimal allocation,
#                 random blocks, random entries (BLUPs), SpATS spatial model
#   3. [prep]     Augmented p-rep designs, search method (Williams, Piepho and
#                 Whitaker, 2011): structure, connectedness, resolvability,
#                 (0,1) concurrence, efficiency factor, analysis vs lm()
#   4. [prepalg]  Augmented p-rep designs, alpha and alpha(2) constructions:
#                 the paper's Examples 1-3, Table 3 parameters, dummies
#   5. [prepan]   p-rep analyses: replicates, fixed and random genotype x
#                 location, per-location analysis, two-stage weighted
#                 analysis (Mohring and Piepho, 2009) with block or SpATS
#                 spatial stage 1
#   6. [rrc]      Resolvable (spatial) row-column designs (Williams, John and
#                 Whitaker, 2006): Tables 1-3, generator, latinized and
#                 t-latinized designs, fillers, RC / LV + RC / AR + RC REML
#   7. [env]      Multi-environment tools: Breedbase files, per-environment
#                 ANOVA, Hartley's Fmax, Lin (1982) environment grouping, maps
#   8. [calib]    Statistical calibration: type I error, power, CI coverage
#   9. [repro]    Reproducibility of every generator
#  10. [timing]   Timing on large trials
#
# Usage (from the package root):
#   Rscript scripts/run_all_checks.R                   # full run
#   Rscript scripts/run_all_checks.R --quick           # faster demo
#   Rscript scripts/run_all_checks.R --only rrc,prepalg
#   Rscript scripts/run_all_checks.R --out results.csv
#
# Exits with status 1 if any check fails, so it can be used in CI.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
out_file <- if ("--out" %in% args) args[which(args == "--out") + 1] else NULL
all_sections <- c("rc", "rcflex", "prep", "prepalg", "prepan", "rrc", "env", "calib", "repro", "timing")
only <- if ("--only" %in% args) strsplit(args[which(args == "--only") + 1], ",")[[1]] else all_sections
unknown <- setdiff(only, all_sections)
if (length(unknown) > 0) stop("Unknown section(s): ", paste(unknown, collapse = ", "))
run <- function(key) key %in% only

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(breedtools))
}
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
has_spats <- requireNamespace("SpATS", quietly = TRUE)

set.seed(20261102)
n_sims <- if (quick) 100 else 400
grDevices::pdf(NULL)
script_start <- Sys.time()

# -----------------------------------------------------------------------------
# Result recorder and helpers
# -----------------------------------------------------------------------------
results <- list()

record <- function(section, scenario, check, pass, detail = "") {
  results[[length(results) + 1L]] <<- data.frame(
    section = section,
    scenario = scenario,
    check = check,
    pass = isTRUE(pass),
    detail = as.character(detail)[1],
    stringsAsFactors = FALSE
  )
  invisible(isTRUE(pass))
}

# Runs `expr`; records a failure (instead of aborting) if it errors.
safely <- function(section, scenario, check, expr) {
  tryCatch(expr, error = function(e) {
    record(section, scenario, check, FALSE, paste("error:", conditionMessage(e)))
    NULL
  })
}

banner <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")
say <- function(...) cat("  ", sprintf(...), "\n", sep = "")
close_to <- function(a, b, tol = 1e-4) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol, scale = 1))

expect_error_msg <- function(sec, scenario, expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  record(sec, scenario, paste("errors with", shQuote(pattern)), !is.na(msg) && grepl(pattern, msg), msg)
}

expect_warning_msg <- function(sec, scenario, expr, pattern) {
  msg <- NA_character_
  withCallingHandlers(force(expr), warning = function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  record(sec, scenario, paste("warns with", shQuote(pattern)), !is.na(msg) && grepl(pattern, msg), msg)
}

# Reference sequential ANOVA from stats::lm(); terms in the given order.
reference_anova <- function(formula, data) {
  tab <- suppressWarnings(stats::anova(stats::lm(formula, data = data)))
  data.frame(df = tab$Df, ss = tab$`Sum Sq`, p = tab$`Pr(>F)`)
}

# Adjusted means from a full lm() fit: genotype effect plus the average design
# effect over the observed plots, with SEs from vcov(). Returns NULL when lm()
# drops aliased columns (then the comparison is skipped).
reference_adjusted <- function(formula, data, genotype = "g") {
  fit <- stats::lm(formula, data)
  if (anyNA(stats::coef(fit))) return(NULL)
  X <- stats::model.matrix(fit)
  gcols <- which(attr(X, "assign") == match(genotype, attr(stats::terms(fit), "term.labels")))
  lv <- levels(data[[genotype]])
  L <- matrix(colMeans(X), nrow = length(lv), ncol = ncol(X), byrow = TRUE)
  L[, gcols] <- 0
  L[cbind(match(colnames(X)[gcols], paste0(genotype, lv)), gcols)] <- 1
  V <- stats::vcov(fit)
  data.frame(genotype = lv, mean = drop(L %*% stats::coef(fit)), se = sqrt(rowSums((L %*% V) * L)))
}

compare_adjusted <- function(sec, scenario, tag, a, ref) {
  if (is.null(ref)) return(invisible(NULL))
  got <- a$adjusted_means[match(ref$genotype, a$adjusted_means$genotype), ]
  record(sec, scenario, paste("adjusted means equal lm(),", tag),
         isTRUE(all.equal(got$adjusted_mean, ref$mean, tolerance = 1e-7)))
  record(sec, scenario, paste("adjusted mean SEs equal lm(),", tag),
         isTRUE(all.equal(got$std_error, ref$se, tolerance = 1e-7)))
}

# Correlation of estimated and true genotype values.
accuracy <- function(means, truth, genotypes = names(truth)) {
  stats::cor(means$adjusted_mean[match(genotypes, means$genotype)], truth[genotypes], use = "complete.obs")
}

# -----------------------------------------------------------------------------
# Row-column helpers (sections 1, 2, 8, 9, 10)
# -----------------------------------------------------------------------------
rc_label <- function(R, C, rb, cb, nc) sprintf("%dx%d field, %dx%d blocks, %d checks", R, C, rb, cb, nc)

rc_make <- function(R, C, rb, cb, nc, seed = 1, n_candidates = 50, ...) {
  dots <- list(...)
  layout_args <- dots[intersect(names(dots), c("checks_per_block", "rows_per_rowgroup", "cols_per_colgroup"))]
  cap <- do.call(augmented_row_column_capacity, c(list(R, C, rb, cb, nc), layout_args))
  augmented_row_column_design(
    treatments = paste0("G", seq_len(cap$entry_plots)),
    controls = paste0("Check", seq_len(nc)),
    rows_in_field = R, cols_in_field = C,
    rows_per_block = rb, cols_per_block = cb,
    n_candidates = n_candidates, seed = seed, ...
  )
}

# Structural checks valid for any augmented row-column design.
rc_check_general <- function(design, scenario, sec = "row-column structure") {
  d <- design$design
  ck <- d[d$is_control == 1, ]
  en <- d[d$is_control == 0, ]
  spec <- design$input
  pb <- spec$checks_per_block
  R <- max(d$row)
  C <- max(d$col)
  ok <- c(
    record(sec, scenario, "plot count = rows x cols", nrow(d) == R * C),
    record(sec, scenario, "every (row, col) used exactly once", !anyDuplicated(paste(d$row, d$col))),
    record(sec, scenario, "plot numbers are 1..n", setequal(d$plots, seq_len(R * C))),
    record(sec, scenario, "row/column groups have the requested sizes",
           identical(as.integer(table(d$rowgroup) / C), as.integer(spec$rows_per_rowgroup)) &&
             identical(as.integer(table(d$colgroup) / R), as.integer(spec$cols_per_colgroup))),
    record(sec, scenario, "check plots per block as requested",
           identical(as.vector(table(factor(ck$block, levels = seq_along(pb)))), as.vector(t(pb)))),
    record(sec, scenario, "no check more than ceiling(pb / vc) times in a block",
           all(table(ck$block, ck$all_entries) <= max(ceiling(pb / length(spec$controls))))),
    record(sec, scenario, "every field row has a check", setequal(unique(ck$row), seq_len(R))),
    record(sec, scenario, "every field column has a check", setequal(unique(ck$col), seq_len(C))),
    record(sec, scenario, "each entry appears exactly once",
           !anyDuplicated(en$all_entries) && setequal(en$all_entries, spec$treatments)),
    record(sec, scenario, "breedbase_design matches design",
           identical(design$breedbase_design$all_entries, d$all_entries))
  )
  o <- d[order(d$plots), ]
  if (identical(spec$plot_type, "serpentine")) {
    ok <- c(ok, record(sec, scenario, "serpentine plots are adjacent",
                       all(abs(diff(o$row)) + abs(diff(o$col)) == 1)))
  }
  all(ok)
}

# Extra properties of the Latin construction.
rc_check_latin <- function(design, scenario, sec = "row-column structure") {
  ck <- design$design[design$design$is_control == 1, ]
  record(sec, scenario, "checks in a block use distinct rows and cols",
         all(vapply(split(ck, ck$block), function(b) !anyDuplicated(b$row) && !anyDuplicated(b$col), logical(1))))
  record(sec, scenario, "same check never twice in a field row", !any(table(ck$all_entries, ck$row) > 1))
  record(sec, scenario, "same check never twice in a field column", !any(table(ck$all_entries, ck$col) > 1))
}

# Analysis checks: ANOVA, df and adjusted means against a full lm() fit.
rc_check_analysis <- function(design, scenario, error_df_expected, sec = "row-column analysis") {
  d <- design$design
  geno_eff <- stats::setNames(stats::rnorm(length(unique(d$all_entries)), sd = 2), unique(d$all_entries))
  d$y <- 50 + 0.4 * d$row - 0.2 * d$col + d$block %% 3 + geno_eff[d$all_entries] + stats::rnorm(nrow(d))
  formula <- y ~ rg + rn + cg + cn + bk + g
  to_ref <- function(x) {
    data.frame(y = x$y, rg = factor(x$rowgroup), rn = interaction(x$rowgroup, x$row, drop = TRUE),
               cg = factor(x$colgroup), cn = interaction(x$colgroup, x$col, drop = TRUE),
               bk = interaction(x$rowgroup, x$colgroup, drop = TRUE), g = factor(x$all_entries))
  }
  for (missing in c(FALSE, TRUE)) {
    x <- d
    tag <- if (missing) "missing plots" else "complete data"
    if (missing) x$y[sample(nrow(x), max(3, round(0.03 * nrow(x))))] <- NA
    a <- safely(sec, scenario, paste("analysis runs,", tag), analyze_augmented_row_column_design(x, "y"))
    if (is.null(a)) next
    ref <- reference_anova(formula, to_ref(x[!is.na(x$y), ]))
    record(sec, scenario, paste("df equal lm(),", tag), isTRUE(all.equal(a$anova$df, ref$df)),
           paste(a$anova$df, collapse = ","))
    record(sec, scenario, paste("sums of squares equal lm(),", tag),
           isTRUE(all.equal(a$anova$sum_sq, ref$ss, tolerance = 1e-8)))
    record(sec, scenario, paste("p-values equal lm(),", tag),
           isTRUE(all.equal(a$anova$p_value[1:6], ref$p[1:6], tolerance = 1e-6)))
    record(sec, scenario, paste("residuals = y - fitted,", tag),
           isTRUE(all.equal(a$plot_results$residual, x$y - a$plot_results$fitted)))
    compare_adjusted(sec, scenario, tag, a, reference_adjusted(formula, to_ref(x[!is.na(x$y), ])))
    if (!missing) {
      record(sec, scenario, "error df = capacity error_df", a$diagnostics$anova_error_df == error_df_expected,
             paste(a$diagnostics$anova_error_df, "vs", error_df_expected))
      record(sec, scenario, "design-formula df = fitted df", isTRUE(all.equal(a$anova_df$display_df, a$anova$df)))
      record(sec, scenario, "no aliased design effects", !a$diagnostics$rank_deficient)
    }
  }
}

# Simulated row-column trial with real row, column and genotype effects and a
# smooth trend, used for the mixed and spatial analyses.
rc_effects_trial <- function(design, seed = 3) {
  x <- design$design
  set.seed(seed)
  g <- stats::setNames(stats::rnorm(length(unique(x$all_entries)), sd = 2), unique(x$all_entries))
  x$rn <- interaction(x$rowgroup, x$row)
  x$cn <- interaction(x$colgroup, x$col)
  x$y <- 20 + g[x$all_entries] + stats::rnorm(nlevels(x$rn))[as.integer(x$rn)] +
    stats::rnorm(nlevels(x$cn))[as.integer(x$cn)] + 3 * sin(x$row / 5) + 2 * cos(x$col / 4) +
    stats::rnorm(nrow(x))
  list(data = x, truth = g, entries = x$all_entries[x$is_control == 0])
}

# -----------------------------------------------------------------------------
# p-rep helpers (sections 3, 4, 5, 8)
# -----------------------------------------------------------------------------
prep_connected <- function(d) {
  vapply(split(d, d$location), function(x) {
    prep_block_components(x$all_entries, x$block, length(unique(x$block))) == 1L
  }, logical(1))
}

# The structure of the paper's constructions: every entry duplicated at exactly
# one location, two replicates per location with every p-rep entry once,
# at most m p-rep entries per block, (0,1) concurrence, connected locations.
prep_check_paper_structure <- function(design, sec, scen, m) {
  d <- design$design
  n_loc <- length(unique(d$location))
  record(sec, scen, "every entry duplicated at exactly one location",
         all(table(d$all_entries) == n_loc + 1L))
  pp <- d[d$entry_type == "p_rep_entry_at_this_location", ]
  record(sec, scen, "resolvable: each p-rep entry once per replicate",
         all(table(paste(pp$location, pp$all_entries), pp$replicate) == 1))
  record(sec, scen, sprintf("at most m = %d p-rep entries per block", m),
         max(table(interaction(pp$location, pp$block, drop = TRUE))) <= m)
  record(sec, scen, "(0,1) concurrence", design$summary$max_pairwise_concurrence == 1)
  record(sec, scen, "every location connected", design$summary$n_disconnected_locations == 0)
  record(sec, scen, "no entry twice in a block", design$summary$same_entry_same_block == 0)
}

# =============================================================================
# 1. Augmented row-column designs: Latin construction
# =============================================================================
if (run("rc")) {
  banner("1. Augmented row-column designs (Latin construction)")
  cat("Every layout in the sweep is either generated (Latin construction when it\n",
      "applies, D-optimal search otherwise) or rejected with a clear error.\n", sep = "")

  dims_rows <- if (quick) c(8, 12, 16, 24) else c(8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 36)
  block_sizes <- if (quick) c(3, 4, 6, 8) else c(3, 4, 5, 6, 8, 10)
  checks <- if (quick) 2:5 else 2:6
  sweep <- expand.grid(R = dims_rows, C = dims_rows, rb = block_sizes, cb = block_sizes, nc = checks)
  sweep <- sweep[sweep$R %% sweep$rb == 0 & sweep$C %% sweep$cb == 0, ]

  n_latin <- 0L
  n_optimal <- 0L
  n_rejected <- 0L
  max_optimal <- if (quick) 6L else 25L
  sweep_time <- system.time({
    for (i in seq_len(nrow(sweep))) {
      s <- sweep[i, ]
      scen <- rc_label(s$R, s$C, s$rb, s$cb, s$nc)
      cap <- augmented_row_column_capacity(s$R, s$C, s$rb, s$cb, s$nc)
      if (cap$entry_plots < 1) next
      nbr <- s$R / s$rb
      nbc <- s$C / s$cb
      feasible <- cap$error_df > 0 && cap$every_row_has_check && cap$every_col_has_check
      latin <- nbc <= s$rb && nbr <= s$cb && s$nc <= s$rb && s$nc <= s$cb
      if (feasible && !latin) {
        if (n_optimal >= max_optimal) next
        n_optimal <- n_optimal + 1L
      }
      design <- tryCatch(rc_make(s$R, s$C, s$rb, s$cb, s$nc, seed = i, n_candidates = if (latin) 20 else 2),
                         error = function(e) conditionMessage(e))
      if (!feasible) {
        n_rejected <- n_rejected + 1L
        record("row-column sweep", scen, "infeasible layout is rejected", is.character(design),
               if (is.character(design)) substr(design, 1, 60) else "design was generated")
        next
      }
      if (is.character(design)) {
        record("row-column sweep", scen, "feasible layout is generated", FALSE, design)
        next
      }
      record("row-column sweep", scen, "feasible layout is generated", TRUE)
      record("row-column sweep", scen, "method chosen as expected",
             design$summary$method == if (latin) "latin" else "optimal", design$summary$method)
      rc_check_general(design, scen)
      if (latin) {
        n_latin <- n_latin + 1L
        rc_check_latin(design, scen)
        if (n_latin %% 5 == 1) rc_check_analysis(design, scen, cap$error_df)
      } else if (n_optimal %% 3 == 1) {
        rc_check_analysis(design, scen, cap$error_df)
      }
    }
  })[["elapsed"]]
  say("%d Latin layouts, %d D-optimal layouts, %d infeasible layouts rejected (%.1fs)",
      n_latin, n_optimal, n_rejected, sweep_time)

  for (cfg in list(c(12, 12, 4, 4, 4), c(36, 18, 6, 6, 4), c(24, 16, 6, 4, 4), c(40, 40, 8, 8, 5))) {
    scen <- do.call(rc_label, as.list(cfg))
    design <- safely("row-column named", scen, "design generated", do.call(rc_make, as.list(cfg)))
    if (is.null(design)) next
    cap <- do.call(augmented_row_column_capacity, as.list(unname(cfg)))
    rc_check_general(design, scen)
    rc_check_latin(design, scen)
    rc_check_analysis(design, scen, cap$error_df)
    say("%-38s entries %4d  checks %3d  error df %3d  D = %.4f", scen, cap$entry_plots, cap$check_plots,
        cap$error_df, design$summary$d_criterion)
  }

  sec <- "row-column validation"
  expect_error_msg(sec, "wrong treatment count",
    augmented_row_column_design(paste0("G", 1:100), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
    "Expected 108 treatments")
  expect_error_msg(sec, "duplicated treatments",
    augmented_row_column_design(c(paste0("G", 1:107), "G1"), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
    "duplicated")
  expect_error_msg(sec, "treatment also a control",
    augmented_row_column_design(c(paste0("G", 1:107), "C1"), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
    "must not contain the same name")
  expect_error_msg(sec, "field not divisible by block", augmented_row_column_capacity(10, 12, 4, 4, 4), "divisible")
  expect_error_msg(sec, "no error df",
    augmented_row_column_design(paste0("G", 1:28), paste0("C", 1:2), 6, 6, 3, 3, n_candidates = 2),
    "no error degrees of freedom")
  expect_error_msg(sec, "rows without a check",
    augmented_row_column_design(paste0("G", 1:92), paste0("C", 1:2), 10, 10, 5, 5, n_candidates = 2),
    "Not every field row can contain a check")

  design <- rc_make(12, 12, 4, 4, 4)
  for (fill in c("type", "block", "rowgroup", "colgroup", "is_control")) {
    p <- safely("row-column plot", fill, "plot draws",
                plot_augmented_row_column_design(design, fill = fill, show_block_labels = TRUE))
    if (!is.null(p)) record("row-column plot", paste("fill =", fill), "plot draws", nrow(p) == 144)
  }
  tmp <- tempfile(fileext = ".txt")
  invisible(rc_make(12, 12, 4, 4, 4, output_file = tmp))
  record("row-column export", "tab-delimited file", "written and readable",
         nrow(read_breedbase_file(tmp)) == 144)
}

# =============================================================================
# 2. Flexible augmented row-column designs (Piepho and Williams, 2016)
# =============================================================================
if (run("rcflex")) {
  banner("2. Augmented row-column designs: Piepho and Williams (2016)")
  sec <- "row-column (Piepho & Williams)"
  pb3 <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))

  cat("Table 2 of the paper (entries / check plots / error df):\n")
  ex1 <- augmented_row_column_capacity(6, 9, 3, 3, n_controls = 3)
  ex2 <- augmented_row_column_capacity(36, 18, 6, 6, n_controls = 4)
  ex3 <- augmented_row_column_capacity(n_controls = 4, rows_per_rowgroup = c(6, 5, 6, 5),
                                       cols_per_colgroup = c(3, 4, 3), checks_per_block = pb3)
  tab2 <- data.frame(example = 1:3, entries = c(ex1$entry_plots, ex2$entry_plots, ex3$entry_plots),
                     check_plots = c(ex1$check_plots, ex2$check_plots, ex3$check_plots),
                     error_df = c(ex1$error_df, ex2$error_df, ex3$error_df),
                     paper_error_df = c(0, 6, 7))
  print(tab2, row.names = FALSE)
  record(sec, "Table 2", "entries 36 / 576 / 173", identical(tab2$entries, c(36L, 576L, 173L)))
  record(sec, "Table 2", "error df 0 / 6 / 7", all(tab2$error_df == tab2$paper_error_df))
  record(sec, "Example 1", "row-group condition 7 <= 9", ex1$min_rowgroup_slack == 2)
  record(sec, "Example 2", "row-group 11 <= 12, column-group 14 <= 24",
         ex2$min_rowgroup_slack == 1 && ex2$min_colgroup_slack == 10)
  expect_error_msg(sec, "Example 1 (no error df)",
    augmented_row_column_design(paste0("E", 1:36), LETTERS[1:3], 6, 9, 3, 3, n_candidates = 2),
    "no error degrees of freedom")

  d3 <- safely(sec, "Example 3", "design generated",
               augmented_row_column_design(paste0("E", 1:173), LETTERS[1:4], rows_per_rowgroup = c(6, 5, 6, 5),
                                           cols_per_colgroup = c(3, 4, 3), checks_per_block = pb3,
                                           n_candidates = 3, seed = 1))
  if (!is.null(d3)) {
    record(sec, "Example 3", "method optimal", d3$summary$method == "optimal")
    rc_check_general(d3, "Example 3 (unequal groups, 3-5 checks per block)", sec)
    rc_check_analysis(d3, "Example 3 (unequal groups, 3-5 checks per block)", 7, sec)
    say("Example 3: 22 x 10 field, %d check plots, %d entries, error df %d, D = %.4f",
        d3$summary$check_plots, d3$summary$entry_plots, d3$summary$error_df, d3$summary$d_criterion)
  }

  inc <- safely(sec, "incomplete blocks", "design generated",
                augmented_row_column_design(paste0("E", 1:279), paste0("C", 1:6), 18, 18, 6, 6,
                                            checks_per_block = 5, n_candidates = 2, seed = 1))
  if (!is.null(inc)) {
    rc_check_general(inc, "incomplete blocks: 6 checks, 5 per block", sec)
    ck <- inc$design[inc$design$is_control == 1, ]
    record(sec, "incomplete blocks", "no check twice in a block", all(table(ck$block, ck$all_entries) <= 1))
  }
  extra <- safely(sec, "extra replication", "design generated",
                  augmented_row_column_design(paste0("E", 1:90), LETTERS[1:4], 12, 12, 4, 4,
                                              checks_per_block = 6, n_candidates = 2, seed = 1))
  if (!is.null(extra)) {
    rc_check_general(extra, "extra replication: 4 checks, 6 per block", sec)
    ck <- extra$design[extra$design$is_control == 1, ]
    record(sec, "extra replication", "every check in every block", all(table(ck$block, ck$all_entries) >= 1))
    rc_check_analysis(extra, "extra replication: 4 checks, 6 per block", extra$summary$error_df, sec)
  }

  latin <- rc_make(36, 18, 6, 6, 4, n_candidates = 100)
  opt <- safely(sec, "Example 2", "optimal design generated",
                rc_make(36, 18, 6, 6, 4, method = "optimal", n_candidates = 3))
  if (!is.null(opt)) {
    say("Example 2 D-criterion: Latin %.4f, D-optimal search %.4f", latin$summary$d_criterion, opt$summary$d_criterion)
    record(sec, "Example 2", "D-optimal search within 3% of Latin on the D-criterion",
           opt$summary$d_criterion >= 0.97 * latin$summary$d_criterion)
    record(sec, "Example 2", "no neighbouring checks (both methods)",
           opt$summary$adjacent_checks_in_blocks + opt$summary$adjacent_checks_across_blocks == 0 &&
             latin$summary$adjacent_checks_in_blocks == 0)
    rc_check_general(opt, "Example 2 (optimal)", sec)
  }
  expect_error_msg(sec, "latin with incomplete blocks",
    augmented_row_column_design(paste0("G", 1:117), paste0("C", 1:4), 12, 12, 4, 4,
                                checks_per_block = 3, method = "latin"), "latin")
  expect_error_msg(sec, "wrong checks_per_block shape",
    augmented_row_column_capacity(12, 12, 4, 4, 4, checks_per_block = matrix(4, 2, 2)), "3 x 3 matrix")

  cat("\nAnalyses of a 36 x 18 trial with row, column and spatial effects\n",
      "(correlation of estimated and true entry effects):\n", sep = "")
  tr <- rc_effects_trial(latin)
  fx <- analyze_augmented_row_column_design(tr$data, "y")
  say("%-34s r = %.3f  mean SED %.3f", "fixed (Eq. 1)", accuracy(fx$adjusted_means, tr$truth, tr$entries),
      fx$diagnostics$mean_sed)
  if (has_lme4) {
    rb <- analyze_augmented_row_column_design(tr$data, "y", random_blocks = TRUE)
    say("%-34s r = %.3f  mean SED %.3f", "random blocks (REML)", accuracy(rb$adjusted_means, tr$truth, tr$entries),
        rb$diagnostics$mean_sed)
    x <- tr$data
    x$gen <- factor(x$all_entries)
    x$bl <- factor(x$block)
    ref <- suppressMessages(lme4::lmer(y ~ factor(rowgroup) + factor(colgroup) + gen + (1 | rn) + (1 | cn) + (1 | bl), data = x))
    fe <- lme4::fixef(ref)
    lv <- levels(x$gen)
    got <- rb$adjusted_means$adjusted_mean[match(lv, rb$adjusted_means$genotype)]
    record(sec, "random blocks", "genotype contrasts equal lme4",
           isTRUE(all.equal(got - got[1], unname(c(0, fe[paste0("gen", lv[-1])])), tolerance = 1e-4)))
    record(sec, "random blocks", "variance components equal lme4",
           isTRUE(all.equal(sort(rb$variance_components$variance), sort(as.data.frame(lme4::VarCorr(ref))$vcov),
                            tolerance = 1e-4)))
    record(sec, "random blocks", "ANOVA is still the fixed-effects table",
           identical(rb$anova_formatted, fx$anova_formatted))
    record(sec, "random blocks", "more accurate than the fixed model",
           accuracy(rb$adjusted_means, tr$truth, tr$entries) > accuracy(fx$adjusted_means, tr$truth, tr$entries))

    re <- analyze_augmented_row_column_design(tr$data, "y", random_blocks = TRUE, entries = "random")
    say("%-34s r = %.3f  entry variance %.2f (true 4)", "random blocks + random entries", accuracy(re$adjusted_means, tr$truth, tr$entries),
        re$variance_components$variance[re$variance_components$component == "entries"])
    x$conv <- as.numeric(x$is_control == 0)
    x$chk <- factor(ifelse(x$is_control == 1, x$all_entries, ".e"))
    x$ent <- factor(ifelse(x$is_control == 1, ".c", x$all_entries))
    ref2 <- suppressMessages(lme4::lmer(y ~ factor(rowgroup) + factor(colgroup) + chk + (1 | rn) + (1 | cn) + (1 | bl) +
                                          (0 + conv | ent), data = x))
    u <- lme4::ranef(ref2)$ent[tr$entries, 1]
    ours <- re$adjusted_means$adjusted_mean[match(tr$entries, re$adjusted_means$genotype)]
    record(sec, "random entries", "BLUPs equal lme4 ranef()",
           isTRUE(all.equal(ours - mean(ours), u - mean(u), tolerance = 1e-4)))
    record(sec, "random entries", "message says entries are BLUPs", grepl("BLUP", re$message))
  } else {
    cat("  (lme4 not installed: mixed-model checks skipped)\n")
  }
  if (has_spats) {
    sp <- analyze_augmented_row_column_design(tr$data, "y", spatial = TRUE)
    say("%-34s r = %.3f  mean SED %.3f  effective dim %.1f", "spatial (SpATS P-splines)",
        accuracy(sp$adjusted_means, tr$truth, tr$entries), sp$diagnostics$mean_sed, sp$diagnostics$spatial_effective_dim)
    record(sec, "spatial", "SpATS fit is more accurate than the fixed model under a trend",
           accuracy(sp$adjusted_means, tr$truth, tr$entries) > accuracy(fx$adjusted_means, tr$truth, tr$entries))
    expect_error_msg(sec, "spatial with random entries",
      analyze_augmented_row_column_design(tr$data, "y", spatial = TRUE, entries = "random"), "only available")
  } else {
    cat("  (SpATS not installed: spatial checks skipped)\n")
  }
}

# =============================================================================
# 3. Augmented p-rep designs: search method
# =============================================================================
if (run("prep")) {
  banner("3. Augmented p-rep designs: search method (Williams, Piepho and Whitaker, 2011)")

  prep_scenarios <- list(
    list(name = "12 entries, 3 locs, block 5",
         args = list(treatments = paste0("G", 1:12), locations = c("Loc1", "Loc2", "Loc3"), block_size = 5, field_cols = 5)),
    list(name = "12 entries, uneven per-location settings",
         args = list(treatments = paste0("G", 1:12), locations = c("Loc1", "Loc2", "Loc3"),
                     block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 6),
                     n_blocks_per_location = c(Loc1 = 4, Loc2 = 4, Loc3 = 4),
                     n_duplicate_per_location = c(Loc1 = 4, Loc2 = 2, Loc3 = 4),
                     replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3, Loc3 = 2), field_cols = 8)),
    list(name = "60 entries, 2 locs, prep_prop 0.3", expect_01 = TRUE,
         args = list(treatments = paste0("G", 1:60), locations = c("A", "B"), block_size = 8, prep_prop = 0.3, field_cols = 10)),
    list(name = "60 entries, 2 locs, resolvable", expect_01 = TRUE,
         args = list(treatments = paste0("G", 1:60), locations = c("A", "B"), block_size = 8, prep_prop = 0.3,
                     resolvable = TRUE, field_cols = 10)),
    list(name = "100 entries, 4 locs, default split",
         args = list(treatments = paste0("G", 1:100), locations = paste0("L", 1:4), block_size = 10, field_cols = 10)),
    list(name = "48 entries, 2 locs, block layout",
         args = list(treatments = paste0("G", 1:48), locations = c("A", "B"), block_size = 6,
                     n_duplicate_per_location = 12, rows_per_block = 2, blocks_per_row = 2)),
    list(name = "120 entries, 6 locs, blocks of 8 (maize)", expect_01 = TRUE,
         args = list(treatments = sprintf("G%03d", 1:120), locations = paste0("L", 1:6), block_size = 8, field_cols = 8)),
    list(name = "200 entries, 1 loc, 25% replicated", expect_01 = TRUE,
         args = list(treatments = paste0("G", 1:200), locations = "Station", block_size = 10, prep_prop = 0.25, field_cols = 20))
  )

  for (sc in prep_scenarios) {
    sec <- "p-rep structure"
    scen <- sc$name
    design <- safely(sec, scen, "design generated",
                     do.call(augmented_prep_design, c(sc$args, list(n_candidates = if (quick) 50 else 200, seed = 7))))
    if (is.null(design)) next
    d <- design$design
    inp <- design$input
    reps <- inp$replications_per_repeated_entry

    record(sec, scen, "every entry at every location",
           all(table(factor(d$all_entries, inp$treatments), factor(d$location, inp$locations)) >= 1))
    rep_ok <- TRUE
    for (loc in inp$locations) {
      x <- d[d$location == loc, ]
      counts <- table(x$all_entries)
      repeated <- unique(x$all_entries[x$repeated_at == loc])
      rep_ok <- rep_ok && all(counts == ifelse(names(counts) %in% repeated, reps[[loc]], 1L))
    }
    record(sec, scen, "replication per location as requested", rep_ok)
    record(sec, scen, "no entry twice in the same block", !anyDuplicated(paste(d$location, d$block, d$all_entries)))
    sizes <- table(paste(d$location, d$block))
    by_loc <- split(as.integer(sizes), sub(" .*", "", names(sizes)))
    record(sec, scen, "block sizes within a location differ by <= 1",
           all(vapply(by_loc, function(s) max(s) - min(s) <= 1, logical(1))))
    pp <- d[d$entry_type == "p_rep_entry_at_this_location", ]
    if (nrow(pp) > 0) {
      per_block <- table(interaction(pp$location, pp$block, drop = TRUE))
      record(sec, scen, "p-rep plots spread evenly over blocks (range <= 2)", max(per_block) - min(per_block) <= 2,
             paste(range(per_block), collapse = "-"))
    }
    record(sec, scen, "every location connected", all(prep_connected(d)))
    record(sec, scen, "no pair of entries shares more than two blocks", design$summary$max_pairwise_concurrence <= 2,
           paste("max concurrence", design$summary$max_pairwise_concurrence))
    if (isTRUE(sc$expect_01)) {
      record(sec, scen, "(0,1) design: no pair shares two blocks", design$summary$max_pairwise_concurrence <= 1,
             paste("max concurrence", design$summary$max_pairwise_concurrence))
    }
    record(sec, scen, "efficiency factor in (0, 1]",
           is.na(design$summary$efficiency_factor) || (design$summary$efficiency_factor > 0 && design$summary$efficiency_factor <= 1))
    if (isTRUE(inp$resolvable)) {
      record(sec, scen, "resolvable: each p-rep entry once per replicate",
             all(table(paste(pp$location, pp$all_entries), pp$replicate) == 1))
    }
    record(sec, scen, "field coordinates unique per location", !anyDuplicated(paste(d$location, d$row, d$col)))
    conc <- augmented_prep_concurrence(d, treatments = rownames(design$concurrence_matrix))
    record(sec, scen, "concurrence matrix consistent", identical(unname(conc), unname(design$concurrence_matrix)))
    say("%-44s plots %4d  max conc %d  E = %.4f  (%d candidates)", scen, design$summary$n_plots,
        design$summary$max_pairwise_concurrence, design$summary$efficiency_factor, design$summary$n_candidates_evaluated)

    # Analysis vs lm()
    sec <- "p-rep analysis"
    geno_eff <- stats::setNames(stats::rnorm(length(inp$treatments), sd = 2), inp$treatments)
    loc_eff <- stats::setNames(seq_along(inp$locations) * 3, inp$locations)
    d$y <- 100 + loc_eff[d$location] + (d$block %% 4) + geno_eff[d$all_entries] + stats::rnorm(nrow(d))
    for (missing in c(FALSE, TRUE)) {
      x <- d
      tag <- if (missing) "missing plots" else "complete data"
      if (missing) x$y[sample(nrow(x), max(2, round(0.03 * nrow(x))))] <- NA
      a <- safely(sec, scen, paste("analysis runs,", tag), analyze_augmented_prep_design(x, "y"))
      if (is.null(a)) next
      xr <- x[!is.na(x$y), ]
      ref_data <- data.frame(y = xr$y, loc = factor(xr$location),
                             rp = interaction(xr$location, xr$replicate, drop = TRUE),
                             blk = interaction(xr$location, xr$replicate, xr$block, drop = TRUE),
                             g = factor(xr$all_entries))
      multi <- length(unique(ref_data$loc)) > 1
      has_rep <- nlevels(ref_data$rp) > nlevels(ref_data$loc)
      ref_formula <- stats::as.formula(paste("y ~", paste(c(if (multi) "loc", if (has_rep) "rp", "blk", "g"), collapse = " + ")))
      ref <- reference_anova(ref_formula, ref_data)
      terms <- c(if (multi) "location", if (has_rep) "replicate_nested", "block_nested", "genotype", "Residuals")
      fitted_df <- a$anova$df[match(terms, a$anova$term)]
      record(sec, scen, paste("fitted df equal lm(),", tag), isTRUE(all.equal(fitted_df, ref$df)),
             paste(fitted_df, collapse = ","))
      ss <- a$anova$sum_sq[match(terms, a$anova$term)]
      record(sec, scen, paste("sums of squares equal lm(),", tag), isTRUE(all.equal(ss, ref$ss, tolerance = 1e-8)))
      compare_adjusted(sec, scen, tag, a, reference_adjusted(ref_formula, ref_data))
    }
    p <- safely("p-rep plot", scen, "plot draws", plot_augmented_prep_design(design))
    if (!is.null(p)) record("p-rep plot", scen, "plot draws", TRUE)
  }

  sec <- "p-rep validation"
  expect_warning_msg(sec, "too few p-rep plots",
    augmented_prep_design(paste0("G", 1:60), "S", block_size = 6, prep_prop = 0.05, seed = 1),
    "cannot all be connected")
  expect_error_msg(sec, "resolvable block count",
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), block_size = 4, n_blocks_per_location = 7, resolvable = TRUE),
    "multiple")
  expect_error_msg(sec, "duplicated entries", augmented_prep_design(c("G1", "G1"), "A"), "unique")
  tmp <- tempfile(fileext = ".txt")
  invisible(augmented_prep_design(paste0("G", 1:24), c("A", "B"), block_size = 6, seed = 1, output_file = tmp))
  record("p-rep export", "tab-delimited file", "written and readable", nrow(read_breedbase_file(tmp)) == 72)
}

# =============================================================================
# 4. Augmented p-rep designs: alpha and alpha(2) constructions
# =============================================================================
if (run("prepalg")) {
  banner("4. Augmented p-rep designs: alpha and alpha(2) constructions")
  sec <- "p-rep alpha constructions"

  cat("Feasible parameter sets for 120 entries in 6 locations (cf. Table 3 of the paper):\n")
  print(as.data.frame(augmented_prep_parameters(120, 6, max_block_size = 10)[, c("method", "m", "s", "block_size",
                                                                             "dummies", "feasible")]), row.names = FALSE)
  ex1 <- augmented_prep_parameters(20, 2, max_block_size = 3)
  r1 <- ex1[ex1$method == "alpha" & ex1$m == 2, ]
  record(sec, "parameters", "Example 1: m = 2 gives s = 5, k = 3", r1$s == 5 && r1$block_size == 3)
  ex3 <- augmented_prep_parameters(42, 3, max_block_size = 4)
  r3 <- ex3[ex3$method == "alpha2" & ex3$m == 2, ]
  record(sec, "parameters", "Example 3: alpha2 m = 2 gives s = 7, k = 4, 14 blocks",
         r3$s == 7 && r3$block_size == 4 && r3$blocks_per_location == 14)

  cases <- list(
    list(name = "Example 1: 20 entries, 2 locs, k = 3 (alpha)", method = "alpha", v = 20, l = 2, k = 3, m = 2),
    list(name = "Example 2: 12 entries, 2 locs, k' > s (alpha)", method = "alpha", v = 12, l = 2, k = 3, m = 2),
    list(name = "Example 3: 42 entries, 3 locs, k = 4 (alpha2)", method = "alpha2", v = 42, l = 3, k = 4, m = 2),
    list(name = "Example 3 with alpha", method = "alpha", v = 42, l = 3, k = 4, m = 2),
    list(name = "150 entries, 6 locs, k = 8, dummies (alpha)", method = "alpha", v = 150, l = 6, k = 8, m = 2),
    list(name = "400 entries, 4 locs, k = 10 (alpha2)", method = "alpha2", v = 400, l = 4, k = 10, m = 4)
  )
  for (cs in cases) {
    design <- safely(sec, cs$name, "design generated",
                     augmented_prep_design(paste0("G", seq_len(cs$v)), paste0("L", seq_len(cs$l)), method = cs$method,
                                           block_size = cs$k, n_candidates = 100, seed = 1, field_cols = cs$k))
    if (is.null(design)) next
    prep_check_paper_structure(design, sec, cs$name, design$summary$m)
    record(sec, cs$name, "block sizes differ by <= 1",
           max(design$design$block_size) - min(design$design$block_size) <= 1)
    say("%-46s m %d  s %2d  blocks %3d  dummies %d  E = %.4f", cs$name, design$summary$m, design$summary$s,
        design$summary$n_blocks, design$summary$n_dummies, design$summary$efficiency_factor)
  }

  record(sec, "Example 3 alpha2", "all blocks have (l + 1) m / 2 = 4 plots", {
    d <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha2", block_size = 4, seed = 2)
    all(d$design$block_size == 4)
  })
  one <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha", block_size = 4,
                               n_efficiency_candidates = 1, seed = 9)
  many <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha", block_size = 4,
                                n_efficiency_candidates = 10, seed = 9)
  say("efficiency of first optimal design %.4f, best of 10 optimal designs %.4f",
      one$summary$efficiency_factor, many$summary$efficiency_factor)
  record(sec, "efficiency selection", "best of several optimal designs is at least as efficient",
         many$summary$efficiency_factor >= one$summary$efficiency_factor - 1e-12)

  expect_error_msg(sec, "alpha2 with even s",
    augmented_prep_design(paste0("G", 1:120), paste0("L", 1:6), method = "alpha2", block_size = 8, n_prep_per_block = 2),
    "odd number of blocks")
  expect_error_msg(sec, "prep_prop with alpha",
    augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha", prep_prop = 0.3), "cannot be used")
  expect_warning_msg(sec, "no (0,1)-design for m = 4",
    augmented_prep_design(sprintf("G%03d", 1:120), paste0("L", 1:4), method = "alpha", block_size = 10,
                          n_candidates = if (quick) 10 else 40, seed = 1),
    "did not find a \\(0,1\\)-design")
}

# =============================================================================
# 5. p-rep analyses
# =============================================================================
if (run("prepan")) {
  banner("5. p-rep analyses: replicates, genotype x location, two-stage, spatial")
  sec <- "p-rep analyses"

  design <- augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha", block_size = 4,
                                  field_cols = 7, seed = 3)
  d <- design$design
  set.seed(1)
  d$y <- stats::rnorm(nrow(d))
  d$y[c(5, 60)] <- NA
  a <- analyze_augmented_prep_design(d, "y")
  o <- d[!is.na(d$y), ]
  ref <- reference_anova(y ~ location + interaction(location, replicate) +
                           interaction(location, replicate, block) + all_entries, o)
  record(sec, "replicates within locations", "df and SS equal lm()",
         isTRUE(all.equal(a$anova_formatted$Df[1:5], ref$df)) &&
           isTRUE(all.equal(a$anova_formatted$`Sum Sq`[1:5], ref$ss, tolerance = 1e-8)))

  # Fixed genotype x location: matches lm() and recovers exact effects.
  set.seed(2)
  ge <- stats::setNames(stats::rnorm(42), paste0("G", 1:42))
  gl <- matrix(stats::rnorm(126), 42, 3, dimnames = list(names(ge), c("A", "B", "C")))
  be <- stats::rnorm(200)
  d$z <- 10 + ge[d$all_entries] + gl[cbind(d$all_entries, d$location)] + be[as.integer(factor(d$block_uid))]
  d$w <- d$z + stats::rnorm(nrow(d), sd = 0.5)
  gxl <- analyze_augmented_prep_design(d, "w", genotype_by_location = "fixed")
  ref <- reference_anova(w ~ location + interaction(location, replicate) + interaction(location, replicate, block) +
                           all_entries + location:all_entries, d)
  record(sec, "fixed genotype x location", "df and SS equal lm()",
         isTRUE(all.equal(gxl$anova_formatted$Df[1:6], ref$df)) &&
           isTRUE(all.equal(gxl$anova_formatted$`Sum Sq`[1:6], ref$ss, tolerance = 1e-8)))
  exact <- suppressWarnings(analyze_augmented_prep_design(d, "z", genotype_by_location = "fixed"))
  w_loc <- as.numeric(table(d$location)[c("A", "B", "C")]) / nrow(d)
  truth <- ge + as.vector(gl %*% w_loc)
  got <- exact$adjusted_means$adjusted_mean[match(names(ge), exact$adjusted_means$genotype)]
  record(sec, "fixed genotype x location", "noise-free data recover the effects exactly",
         isTRUE(all.equal(diff(got), unname(diff(truth)), tolerance = 1e-8)))

  if (has_lme4) {
    rd <- augmented_prep_design(paste0("G", 1:60), c("A", "B", "C"), method = "alpha", block_size = 6, seed = 4)$design
    set.seed(1)
    g <- stats::setNames(stats::rnorm(60, sd = 2), paste0("G", 1:60))
    gx <- interaction(rd$all_entries, rd$location, drop = TRUE)
    rd$y <- 50 + g[rd$all_entries] + stats::rnorm(nlevels(gx), sd = 1.5)[as.integer(gx)] + stats::rnorm(nrow(rd))
    ra <- analyze_augmented_prep_design(rd, "y", genotype_by_location = "random")
    rd$gen <- factor(rd$all_entries)
    rd$bn <- interaction(rd$location, rd$replicate, rd$block, drop = TRUE)
    rd$rn <- interaction(rd$location, rd$replicate, drop = TRUE)
    rd$gx <- gx
    lm4 <- suppressMessages(lme4::lmer(y ~ location + rn + bn + gen + (1 | gx), data = rd,
                                       control = lme4::lmerControl(check.rankX = "silent.drop.cols")))
    fe <- lme4::fixef(lm4)
    lv <- levels(rd$gen)
    got <- ra$adjusted_means$adjusted_mean[match(lv, ra$adjusted_means$genotype)]
    record(sec, "random genotype x location", "genotype contrasts equal lme4",
           isTRUE(all.equal(got - got[1], unname(c(0, fe[paste0("gen", lv[-1])])), tolerance = 1e-5)))
    record(sec, "random genotype x location", "variance components equal lme4",
           isTRUE(all.equal(ra$variance_components$variance, as.data.frame(lme4::VarCorr(lm4))$vcov, tolerance = 1e-5)))
    say("random GxL: REML GxL variance %.2f (true 2.25), residual %.2f (true 1)",
        ra$variance_components$variance[1], ra$variance_components$variance[2])
  }

  # Per-location analyses are estimable because locations are connected.
  pl <- augmented_prep_design(paste0("G", 1:60), c("A", "B"), block_size = 6, prep_prop = 0.3, resolvable = TRUE, seed = 8)$design
  pl$y <- stats::rnorm(nrow(pl))
  ok <- vapply(c("A", "B"), function(loc) {
    fit <- analyze_augmented_prep_design(pl[pl$location == loc, ], "y")
    !anyNA(fit$adjusted_means$adjusted_mean) && fit$diagnostics$residual_df > 0
  }, logical(1))
  record(sec, "per-location analysis", "every location estimable on its own", all(ok))

  # Two-stage analysis: stage 2 equals brute-force REML and GLS.
  tsd <- augmented_prep_design(sprintf("G%02d", 1:60), c("A", "B", "C"), method = "alpha", block_size = 6,
                               field_cols = 9, seed = 1)$design
  set.seed(1)
  g <- stats::setNames(stats::rnorm(60, sd = 2), sprintf("G%02d", 1:60))
  gx <- interaction(tsd$all_entries, tsd$location, drop = TRUE)
  sdl <- c(A = 0.5, B = 0.5, C = 4)
  tsd$y <- 50 + g[tsd$all_entries] + stats::rnorm(nlevels(gx))[as.integer(gx)] +
    stats::rnorm(nrow(tsd), sd = sdl[tsd$location])
  two <- analyze_augmented_prep_two_stage(tsd, "y")
  s1 <- lapply(c("A", "B", "C"), function(loc) {
    f <- analyze_augmented_prep_design(tsd[tsd$location == loc, ], "y")
    am <- f$adjusted_means[!is.na(f$adjusted_means$adjusted_mean), ]
    list(loc = loc, means = stats::setNames(am$adjusted_mean, am$genotype), vcov = f$adjusted_means_vcov[am$genotype, am$genotype])
  })
  st <- do.call(rbind, lapply(s1, function(r) data.frame(loc = r$loc, g = names(r$means), y = r$means)))
  st$loc <- factor(st$loc)
  st$g <- factor(st$g, levels = sort(unique(st$g)))
  X <- stats::model.matrix(~ loc + g, st)
  S <- matrix(0, nrow(st), nrow(st))
  off <- 0
  for (r in s1) {
    idx <- off + seq_along(r$means)
    S[idx, idx] <- r$vcov
    off <- off + length(idx)
  }
  reml <- function(s2) {
    V <- S + diag(s2, nrow(S))
    Vi <- solve(V)
    XtVX <- t(X) %*% Vi %*% X
    b <- solve(XtVX, t(X) %*% Vi %*% st$y)
    res <- st$y - X %*% b
    -0.5 * (determinant(V)$modulus + determinant(XtVX)$modulus + t(res) %*% Vi %*% res)
  }
  opt <- stats::optimize(reml, c(0, 20), maximum = TRUE)
  record(sec, "two-stage", "stage-2 REML variance equals brute force",
         isTRUE(all.equal(two$variance_components$variance, opt$maximum, tolerance = 1e-3)))
  Vi <- solve(S + diag(two$variance_components$variance, nrow(S)))
  XtVX <- t(X) %*% Vi %*% X
  b <- solve(XtVX, t(X) %*% Vi %*% st$y)
  L <- matrix(0, 60, ncol(X))
  L[, 1] <- 1
  L[, grep("^loc", colnames(X))] <- 1 / 3
  for (j in 2:60) L[j, match(paste0("g", levels(st$g)[j]), colnames(X))] <- 1
  got <- two$adjusted_means[match(levels(st$g), two$adjusted_means$genotype), ]
  record(sec, "two-stage", "combined means and SEs equal brute-force GLS",
         isTRUE(all.equal(got$adjusted_mean, drop(L %*% b), tolerance = 1e-6)) &&
           isTRUE(all.equal(got$std_error, sqrt(diag(L %*% solve(XtVX) %*% t(L))), tolerance = 1e-6)))
  record(sec, "two-stage", "Hartley test detects the noisy location", isTRUE(two$heterogeneity$heterogeneous),
         sprintf("Fmax %.2f vs %.2f", two$heterogeneity$fmax, two$heterogeneity$critical_5pct))
  say("two-stage: stage-1 residual variances %s; Fmax %.1f (5%% critical %.1f)",
      paste(sprintf("%.2f", two$locations$residual_variance), collapse = ", "),
      two$heterogeneity$fmax, two$heterogeneity$critical_5pct)
  wts <- c(full = "full", diagonal = "diagonal", none = "none")
  acc <- vapply(wts, function(w) accuracy(analyze_augmented_prep_two_stage(tsd, "y", weighting = w)$adjusted_means, g), numeric(1))
  say("two-stage accuracy by weighting: full %.3f, diagonal %.3f, unweighted %.3f", acc[1], acc[2], acc[3])
  record(sec, "two-stage", "weighted analysis at least as accurate as unweighted", acc[["full"]] >= acc[["none"]] - 0.01)
  one_loc <- tsd[tsd$location == "A", ]
  t1 <- analyze_augmented_prep_two_stage(one_loc, "y")
  o1 <- analyze_augmented_prep_design(one_loc, "y")
  record(sec, "two-stage", "one location equals the one-stage analysis",
         isTRUE(all.equal(t1$adjusted_means$adjusted_mean[match(o1$adjusted_means$genotype, t1$adjusted_means$genotype)],
                          o1$adjusted_means$adjusted_mean, tolerance = 1e-8)))

  if (has_spats) {
    tsd$z <- tsd$y + 4 * (sin(tsd$row / 2) + cos(tsd$col / 3))
    blocks <- analyze_augmented_prep_two_stage(tsd, "z")
    spatial <- analyze_augmented_prep_two_stage(tsd, "z", stage1 = "spatial")
    say("field trend added: accuracy blocks stage 1 %.3f, SpATS stage 1 %.3f",
        accuracy(blocks$adjusted_means, g), accuracy(spatial$adjusted_means, g))
    record(sec, "two-stage spatial", "SpATS stage 1 more accurate under a trend",
           accuracy(spatial$adjusted_means, g) > accuracy(blocks$adjusted_means, g))
  }
}

# =============================================================================
# 6. Resolvable (spatial) row-column designs
# =============================================================================
if (run("rrc")) {
  banner("6. Resolvable spatial row-column designs (Williams, John and Whitaker, 2006)")
  sec <- "resolvable row-column"

  table1 <- list(
    rbind(c(19, 34, 27, 3, 15, 12, 2), c(40, 4, 7, 20, 21, 33, 18), c(23, 35, 17, 29, 22, 28, 38),
          c(10, 14, 31, 39, 32, 5, 42), c(13, 6, 16, 36, 11, 37, 41), c(25, 8, 1, 26, 9, 30, 24)),
    rbind(c(38, 34, 26, 33, 7, 10, 11), c(1, 13, 15, 35, 5, 3, 18), c(32, 20, 17, 41, 8, 24, 19),
          c(12, 30, 6, 39, 23, 21, 16), c(4, 31, 42, 25, 2, 37, 29), c(36, 22, 40, 27, 9, 14, 28)),
    rbind(c(20, 26, 23, 37, 32, 27, 35), c(10, 22, 16, 18, 42, 8, 12), c(28, 4, 15, 39, 7, 13, 24),
          c(1, 41, 14, 34, 30, 29, 40), c(2, 5, 36, 17, 25, 21, 11), c(6, 19, 33, 9, 3, 38, 31)))
  table2 <- list(
    rbind(c(11, 8, 33, 23, 10, 6, 29), c(39, 28, 17, 4, 2, 21, 38), c(16, 32, 41, 37, 42, 9, 13),
          c(1, 7, 19, 5, 3, 36, 18), c(27, 14, 34, 30, 40, 31, 25), c(26, 24, 12, 35, 20, 22, 15)),
    rbind(c(6, 2, 31, 20, 13, 14, 5), c(42, 36, 26, 19, 35, 25, 38), c(24, 11, 10, 9, 40, 1, 4),
          c(30, 22, 33, 18, 16, 28, 27), c(21, 12, 32, 15, 17, 37, 8), c(3, 23, 34, 39, 7, 29, 41)),
    rbind(c(5, 28, 20, 30, 8, 7, 26), c(2, 19, 14, 12, 39, 31, 37), c(34, 11, 35, 27, 6, 15, 40),
          c(18, 4, 41, 3, 25, 10, 22), c(24, 33, 21, 13, 36, 16, 29), c(17, 42, 23, 1, 32, 38, 9)))
  as_df <- function(tab) do.call(rbind, lapply(1:3, function(q) {
    data.frame(replicate = q, row_in_rep = rep(1:6, 7), col_in_rep = rep(1:7, each = 6),
               all_entries = paste0("T", as.vector(tab[[q]])), stringsAsFactors = FALSE)
  }))
  d1 <- as_df(table1)
  d2 <- as_df(table2)

  e1 <- resolvable_row_column_efficiency(d1)
  say("Table 1 design: E_rc %.4f (paper 0.6941), E_spatial %.4f (paper 0.9204), U %.4f (paper 0.9268)",
      e1$efficiency_rc, e1$efficiency_spatial, e1$upper_bound_spatial)
  record(sec, "Table 1", "row-column efficiency 0.6941", close_to(e1$efficiency_rc, 0.6941))
  record(sec, "Table 1", "spatial efficiency 0.9204", close_to(e1$efficiency_spatial, 0.9204))
  record(sec, "Table 1", "spatial upper bound 0.9268", close_to(e1$upper_bound_spatial, 0.9268))
  nb <- e1$neighbours
  record(sec, "Table 1", "treatment 1 has 9 neighbours (4 row, 5 column)",
         nb$neighbours[nb$treatment == "T1"] == 9 && nb$row_neighbours[nb$treatment == "T1"] == 4)
  record(sec, "Table 1", "treatments 19 and 38 have 8 neighbours",
         all(nb$neighbours[nb$treatment %in% c("T19", "T38")] == 8))
  record(sec, "Table 1", "9 treatments on corner/edge positions in all 3 replicates", sum(nb$reps_on_edge == 3) == 9)

  cat("  Table 3 (relative efficiency of the LV design to the AR design):\n")
  sets <- rbind(c(35, 60, 10, 0.5, 0.5), c(5, 10, 10, 0.6, 0.5), c(5, 10, 10, 0.2, 0.8), c(70, 60, 1, 0.9, 0.9),
                c(10, 10, 1, 0.9, 0.9), c(35, 30, 0.5, 0.8, 0.8), c(2, 2, 0.5, 0.8, 0.8), c(0, 0, 0.5, 0.8, 0.8))
  paper_re <- c(0.9567, 0.9645, 0.9609, 1.0434, 1.0337, 1.0379, 1.0154, 1.0054)
  sets_to_run <- if (quick) c(1, 4, 8) else seq_len(nrow(sets))
  for (i in sets_to_run) {
    ar <- list(eta_r = sets[i, 1], eta_c = sets[i, 2], psi = sets[i, 3], rho_r = sets[i, 4], rho_c = sets[i, 5])
    re <- resolvable_row_column_efficiency(d2, ar = ar)$ar_criterion / resolvable_row_column_efficiency(d1, ar = ar)$ar_criterion
    cat(sprintf("    set %d: RE %.4f (paper %.4f)\n", i, re, paper_re[i]))
    record(sec, paste("Table 3, set", i), "relative efficiency as published", close_to(re, paper_re[i]))
  }

  laid <- resolvable_row_column_design(paste0("T", 1:42), 3, 6, 7, start_design = d1, n_iter = 1, spatial = FALSE)
  record(sec, "start design", "laying out the Table 1 design keeps its efficiency",
         close_to(laid$summary$efficiency_rc, 0.6941) && close_to(laid$summary$efficiency_spatial, 0.9204))

  gen <- safely(sec, "Example 2", "design generated",
                resolvable_row_column_design(paste0("T", 1:42), 3, 6, 7, n_starts = if (quick) 1 else 3, seed = 1))
  if (!is.null(gen)) {
    s <- gen$summary
    say("Generated Example 2: E_rc %.4f (%.1f%% of the paper's CycDesigN design), E_rows %.4f, E_cols %.4f, E_spatial %.4f",
        s$efficiency_rc, 100 * s$efficiency_rc / 0.6941, s$efficiency_rows, s$efficiency_cols, s$efficiency_spatial)
    d <- gen$design
    record(sec, "Example 2", "resolvable: every treatment once per replicate", all(table(d$all_entries, d$replicate) == 1))
    record(sec, "Example 2", "row-column efficiency within 2% of the paper's design", s$efficiency_rc >= 0.98 * 0.6941,
           sprintf("%.4f", s$efficiency_rc))
    record(sec, "Example 2", "spatial efficiency within 1% of the paper's design", s$efficiency_spatial >= 0.99 * 0.9204,
           sprintf("%.4f", s$efficiency_spatial))
    record(sec, "Example 2", "spatial efficiency below its upper bound", s$efficiency_spatial <= s$upper_bound_spatial + 1e-8)
    rnd <- d
    rnd$all_entries <- unlist(lapply(split(d$all_entries, d$replicate), sample))
    record(sec, "Example 2", "better than a random resolvable design",
           s$efficiency_rc > resolvable_row_column_efficiency(rnd)$efficiency_rc)
  }

  lat <- resolvable_row_column_design(paste0("T", 1:42), 3, 6, 7, latinize = TRUE, n_starts = 1,
                                      n_iter = if (quick) 2000 else NULL, seed = 2)
  record(sec, "latinized", "no treatment twice in a long column",
         lat$summary$latin_ok && all(table(lat$design$all_entries, lat$design$long_col) <= 1))
  tl <- resolvable_row_column_design(paste0("T", 1:30), 6, 5, 6, latinize = TRUE, latin_t = 2,
                                     replicate_layout = "horizontal", n_starts = 1, n_iter = 1500, seed = 3)
  grp <- (tl$design$replicate - 1) %/% 2
  record(sec, "t-latinized (t = 2, side by side)", "latin within each pair of replicates",
         tl$summary$latin_ok && !any(duplicated(paste(grp, tl$design$all_entries, tl$design$row_in_rep))))
  fil <- resolvable_row_column_design(paste0("P", 1:59), 4, 6, 10, latinize = TRUE, n_starts = 1,
                                      n_iter = if (quick) 1500 else 4000, seed = 1)
  record(sec, "Example 1 layout: 59 provenances in 4 x (6 x 10)", "one filler plot per replicate",
         fil$summary$n_fillers == 1 && sum(fil$design$is_filler) == 4)
  say("Latinized 6 x 10 design for 59 provenances: E_rc %.4f, E_spatial %.4f (U %.4f)",
      fil$summary$efficiency_rc, fil$summary$efficiency_spatial, fil$summary$upper_bound_spatial)
  p <- safely(sec, "plot", "plot draws", plot_resolvable_row_column_design(lat))
  if (!is.null(p)) record(sec, "plot", "plot draws", nrow(p) == 126)

  # Analyses
  dd <- lat$design
  set.seed(5)
  g <- stats::setNames(stats::rnorm(42, sd = 2), paste0("T", 1:42))
  dd$y <- 50 + g[dd$all_entries] + c(0, 2, -1)[dd$replicate] + stats::rnorm(18, sd = 1.2)[(dd$replicate - 1) * 6 + dd$row_in_rep] +
    stats::rnorm(21, sd = 0.8)[(dd$replicate - 1) * 7 + dd$col_in_rep] + stats::rnorm(nrow(dd))
  rc <- analyze_resolvable_row_column(dd, "y", model = "rc", compare = FALSE)
  if (has_lme4) {
    ref <- suppressMessages(lme4::lmer(y ~ factor(replicate) + all_entries + (1 | replicate:row_in_rep) +
                                         (1 | replicate:col_in_rep), data = dd))
    vc <- as.data.frame(lme4::VarCorr(ref))
    p <- stats::setNames(rc$variance_parameters$estimate, rc$variance_parameters$parameter)
    record(sec, "RC model", "variance components equal lme4",
           close_to(p[["row_variance"]], vc$vcov[vc$grp == "replicate:row_in_rep"], 1e-3) &&
             close_to(p[["column_variance"]], vc$vcov[vc$grp == "replicate:col_in_rep"], 1e-3) &&
             close_to(p[["sigma2"]], vc$vcov[vc$grp == "Residual"], 1e-3))
    fe <- lme4::fixef(ref)
    lv <- sort(unique(dd$all_entries))
    got <- rc$adjusted_means$adjusted_mean[match(lv, rc$adjusted_means$genotype)]
    record(sec, "RC model", "treatment contrasts equal lme4",
           isTRUE(all.equal(got - got[1], unname(c(0, fe[paste0("all_entries", lv[-1])])), tolerance = 1e-4)))
  }

  # LV + RC on data simulated under the LV model: the REML criterion equals a
  # brute-force REML, the optimum is at least as good as the true parameters,
  # and the likelihood-ratio test prefers LV over RC.
  lvd <- resolvable_row_column_design(paste0("T", 1:64), 4, 8, 8, n_starts = 1, n_iter = 300, spatial = FALSE, seed = 2)$design
  mats <- rrc_variance_matrices(factor(lvd$replicate), lvd$row_in_rep, lvd$col_in_rep, 8, 8)
  truth_lv <- c(eta_r = 6, eta_c = 5, phi_r = 4, phi_c = 3)
  R_ind <- stats::model.matrix(~ 0 + factor(lvd$replicate))
  lin_V <- function(th) diag(nrow(lvd)) + th[["eta_r"]] * mats$eta_r + th[["eta_c"]] * mats$eta_c +
    th[["phi_r"]] * mats$phi_r + th[["phi_c"]] * mats$phi_c
  ch <- chol(lin_V(truth_lv) + 50 * tcrossprod(R_ind))
  X <- stats::model.matrix(~ factor(replicate) + factor(all_entries), lvd)
  brute <- function(th, y) {
    V <- lin_V(th) + 50 * tcrossprod(R_ind)
    Vi <- solve(V)
    XtVX <- t(X) %*% Vi %*% X
    P <- Vi - Vi %*% X %*% solve(XtVX, t(X) %*% Vi)
    m <- nrow(X) - ncol(X)
    quad <- drop(t(y) %*% P %*% y)
    drop(determinant(V)$modulus + determinant(XtVX)$modulus - determinant(crossprod(X))$modulus) +
      m * log(quad / m) + m * (1 + log(2 * pi))
  }
  spec_lv <- list(terms = c("eta_r", "eta_c", "phi_r", "phi_c"), ar = FALSE)
  n_lv <- if (quick) 4 else 10
  est <- matrix(NA, n_lv, 4)
  crit_ok <- opt_ok <- lrt_sig <- logical(n_lv)
  for (i in seq_len(n_lv)) {
    y <- 10 + drop(t(ch) %*% stats::rnorm(nrow(lvd)))
    eng <- rrc_reml_engine(y, X, mats, factor(lvd$replicate), lvd$row_in_rep, lvd$col_in_rep)
    f <- eng$fit(spec_lv)
    est[i, ] <- f$parameters[names(truth_lv)]
    crit_ok[i] <- close_to(f$deviance, brute(f$parameters, y), 1e-6)
    opt_ok[i] <- f$deviance <= eng$neg2ll(truth_lv, spec_lv) + 1e-6
    rc_fit <- eng$fit(list(terms = c("eta_r", "eta_c"), ar = FALSE))
    lrt_sig[i] <- stats::pchisq(rc_fit$deviance - f$deviance, df = 2, lower.tail = FALSE) < 0.05
  }
  say("LV + RC on simulated LV data (%d simulations): mean eta_r %.2f (6), eta_c %.2f (5), phi_r %.2f (4), phi_c %.2f (3)",
      n_lv, mean(est[, 1]), mean(est[, 2]), mean(est[, 3]), mean(est[, 4]))
  say("  LRT prefers LV over RC in %d of %d simulations", sum(lrt_sig), n_lv)
  record(sec, "LV + RC model", "REML criterion equals brute-force REML", all(crit_ok))
  record(sec, "LV + RC model", "REML optimum at least as good as the true parameters", all(opt_ok))
  record(sec, "LV + RC model", "LRT detects the LV process in most simulations", mean(lrt_sig) >= 0.75,
         sprintf("%d / %d", sum(lrt_sig), n_lv))

  # AR data: the model comparison prefers AR + RC over RC.
  ard <- lat$design
  mats <- rrc_variance_matrices(factor(ard$replicate), ard$row_in_rep, ard$col_in_rep, 6, 7)
  Var <- diag(nrow(ard)) + 3 * mats$same_rep * 0.9^mats$dr * 0.9^mats$dc
  set.seed(11)
  ard$y <- 20 + g[ard$all_entries] + drop(t(chol(Var)) %*% stats::rnorm(nrow(ard)))
  cmp <- analyze_resolvable_row_column(ard, "y", model = "ar", long_col = "long_col")
  cat("  Model comparison on data with an AR1 x AR1 field process:\n")
  print(as.data.frame(cmp$model_comparison[, c("model", "n_variance_parameters", "deviance", "AIC", "lrt_p_value")]),
        row.names = FALSE, digits = 4)
  record(sec, "model comparison", "AR + RC has a lower deviance than RC on AR data",
         cmp$model_comparison$deviance[5] < cmp$model_comparison$deviance[1])
  record(sec, "model comparison", "nested LV models never fit worse than RC",
         all(cmp$model_comparison$deviance[2:4] <= cmp$model_comparison$deviance[1] + 1e-4))
  record(sec, "latinized analysis", "long columns fitted as fixed", !anyNA(cmp$adjusted_means$adjusted_mean))

  expect_error_msg(sec, "too many replicates for latinization",
    resolvable_row_column_design(paste0("T", 1:20), 6, 4, 5, latinize = TRUE), "latin")
  expect_error_msg(sec, "too many treatments",
    resolvable_row_column_design(paste0("T", 1:21), 2, 4, 5), "more treatments")
}

# =============================================================================
# 7. Multi-environment tools
# =============================================================================
if (run("env")) {
  banner("7. Multi-environment tools (Breedbase files, Hartley's Fmax, Lin grouping)")
  sec <- "multi-environment"
  pheno_csv <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
  pheno_xlsx <- system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools")
  locs_csv <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
  if (!nzchar(pheno_csv)) {
    pheno_csv <- file.path("inst", "extdata", "breedbase_pheno.csv")
    pheno_xlsx <- file.path("inst", "extdata", "breedbase_pheno.xlsx")
    locs_csv <- file.path("inst", "extdata", "breedbase_locations.csv")
  }
  trait <- "fresh root yield|CO_334:0000013"

  raw <- read_breedbase_file(pheno_csv)
  raw_x <- read_breedbase_file(pheno_xlsx)
  # The workbook holds two of the four trials, one per sheet.
  key <- function(x) paste(x$locationName, x$germplasmName, x$replicate, x[[trait]])
  record(sec, "read files", "Excel sheets combined and equal to the CSV rows of those trials",
         setequal(unique(raw_x$.sheet), c("ITH_2025_yield", "GEN_2025_yield")) &&
           setequal(key(raw_x), key(raw[raw$locationName %in% unique(raw_x$locationName), ])))
  say("Breedbase file: %d plots, %d environments, %d genotypes", nrow(raw), length(unique(raw$locationName)),
      length(unique(raw$germplasmName)))
  record(sec, "column detection", "genotype column guessed",
         identical(guess_column(raw, c("accession_name", "germplasmName")), "germplasmName"))

  pheno <- prepare_breedbase_pheno(raw, trait_col = trait)
  env <- analyze_environment_anova(pheno)
  cat("  Per-environment ANOVA:\n")
  print(as.data.frame(env[, c("environment", "n_genotypes", "ms_error", "mean", "cv_percent", "p_value")]),
        row.names = FALSE, digits = 4)
  ref <- stats::anova(stats::lm(value ~ replicate + genotype,
                                data = pheno[pheno$environment == env$environment[1], ]))
  record(sec, "environment ANOVA", "error mean square equals lm()",
         close_to(env$ms_error[1], ref$`Mean Sq`[nrow(ref)], 1e-8))

  vc <- run_breedbase_variance_check(pheno_csv, trait_col = trait, n_sim = if (quick) 2000 else 20000, export = FALSE)
  say("Hartley's Fmax over all environments: Fmax %.2f, critical %.2f -> %s",
      vc$overall_hartley$fmax_observed, vc$overall_hartley$fmax_critical,
      if (isTRUE(vc$overall_hartley$homogeneous)) "homogeneous" else "heterogeneous")
  record(sec, "Hartley Fmax", "Fmax = max / min error mean square",
         close_to(vc$overall_hartley$fmax_observed, max(vc$env_anova$ms_error) / min(vc$env_anova$ms_error), 1e-8))
  record(sec, "Hartley Fmax", "compatible groups found", nrow(vc$compatible_groups) > 0)

  strat <- run_environment_stratification(pheno_csv, trait = "fresh root yield")
  cat("  Lin (1982) environment groups:\n")
  print(as.data.frame(strat$group_summary[, c("group_id", "locations", "n_env", "ms_ge", "f_value", "p_value")]),
        row.names = FALSE, digits = 4)
  record(sec, "Lin grouping", "two groups: Geneva + Ithaca and Ibadan + Ubiaja",
         nrow(strat$group_summary) == 2 &&
           setequal(strat$group_summary$locations, c("Geneva, Ithaca", "Ibadan, Ubiaja")))
  locations <- read_locations(locs_csv, phenotype_data = strat)
  record(sec, "locations", "coordinates read for every environment", nrow(locations) == 4)
  m <- safely(sec, "map", "map built", map_environment_groups(strat, locations))
  if (!is.null(m)) record(sec, "map", "leaflet map of environment groups", inherits(m, "leaflet"))
  out_dir <- file.path(tempdir(), "breedtools-validate-json")
  invisible(run_environment_stratification(pheno_csv, trait = "fresh root yield", export_json = TRUE, output_dir = out_dir))
  record(sec, "export", "JSON written to output_dir",
         file.exists(file.path(out_dir, "environment_stratification.json")))
}

# =============================================================================
# 8. Statistical calibration
# =============================================================================
if (run("calib")) {
  banner("8. Statistical calibration")
  # Under H0 (no genotype differences) F-test p-values should be Uniform(0,1).
  # Binomial 99.9% interval for n_sims tests is used as the pass band.
  calibrate <- function(label, design_data, analyze, sim_y, p_of, n = n_sims) {
    sec <- "calibration"
    band <- stats::qbinom(c(0.0005, 0.9995), n, 0.05) / n
    p_null <- p_alt <- numeric(n)
    for (i in seq_len(n)) {
      design_data$y <- sim_y(design_data, genotype_sd = 0)
      p_null[i] <- p_of(analyze(design_data))
      design_data$y <- sim_y(design_data, genotype_sd = 3)
      p_alt[i] <- p_of(analyze(design_data))
    }
    size <- mean(p_null < 0.05)
    power <- mean(p_alt < 0.05)
    ks <- stats::ks.test(p_null, "punif")$p.value
    record(sec, label, "type I error ~ 5% under H0", size >= band[1] && size <= band[2],
           sprintf("rate = %.3f (band %.3f-%.3f)", size, band[1], band[2]))
    record(sec, label, "null p-values uniform (KS p > 0.001)", ks > 0.001, sprintf("KS p = %.3g", ks))
    record(sec, label, "power > 0.9 with genotype sd = 3 x error sd", power > 0.9, sprintf("power = %.3f", power))
    say("%-44s size %.3f  power %.3f  KS p %.3g", label, size, power, ks)
  }
  geno_p <- function(a) a$anova$p_value[a$anova$term == "genotype"]

  rc_sim <- function(d, genotype_sd) {
    g <- unique(d$all_entries)
    ge <- stats::setNames(stats::rnorm(length(g), sd = genotype_sd), g)
    10 + 0.5 * d$row + 0.3 * d$col + 2 * (d$block %% 2) + ge[d$all_entries] + stats::rnorm(nrow(d))
  }
  calibrate("row-column 12x12, 4x4 blocks, 4 checks", rc_make(12, 12, 4, 4, 4)$design,
            function(x) analyze_augmented_row_column_design(x, "y"), rc_sim, geno_p)
  ex3 <- augmented_row_column_design(paste0("E", 1:173), LETTERS[1:4], rows_per_rowgroup = c(6, 5, 6, 5),
                                     cols_per_colgroup = c(3, 4, 3),
                                     checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3)),
                                     n_candidates = 2, seed = 1)
  calibrate("row-column Piepho & Williams Example 3", ex3$design,
            function(x) analyze_augmented_row_column_design(x, "y"), rc_sim, geno_p)

  prep_sim <- function(d, genotype_sd) {
    g <- unique(d$all_entries)
    ge <- stats::setNames(stats::rnorm(length(g), sd = genotype_sd), g)
    blk <- interaction(d$location, d$block, drop = TRUE)
    be <- stats::setNames(stats::rnorm(nlevels(blk), sd = 2), levels(blk))
    20 + 5 * as.integer(factor(d$location)) + be[as.character(blk)] + ge[d$all_entries] + stats::rnorm(nrow(d))
  }
  prep_design <- augmented_prep_design(paste0("G", 1:60), c("A", "B"), block_size = 8, prep_prop = 0.3,
                                       field_cols = 10, n_candidates = 20, seed = 3)$design
  calibrate("p-rep 60 entries, 2 locations (search)", prep_design,
            function(x) analyze_augmented_prep_design(x, "y"), prep_sim, geno_p)
  alpha_design <- augmented_prep_design(paste0("G", 1:60), c("A", "B", "C"), method = "alpha", block_size = 6,
                                        seed = 3)$design
  calibrate("p-rep 60 entries, 3 locations (alpha)", alpha_design,
            function(x) analyze_augmented_prep_design(x, "y"), prep_sim, geno_p)

  rrc_design <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5, n_starts = 1, n_iter = 500, seed = 1)$design
  rrc_sim <- function(d, genotype_sd) {
    g <- unique(d$all_entries)
    ge <- stats::setNames(stats::rnorm(length(g), sd = genotype_sd), g)
    rr <- stats::rnorm(12)[(d$replicate - 1) * 4 + d$row_in_rep]
    cc <- stats::rnorm(15)[(d$replicate - 1) * 5 + d$col_in_rep]
    10 + d$replicate + rr + cc + ge[d$all_entries] + stats::rnorm(nrow(d))
  }
  calibrate("resolvable row-column 20 x 3 reps (intra-block F)", rrc_design,
            function(x) analyze_resolvable_row_column(x, "y", compare = FALSE), rrc_sim,
            function(a) a$anova$p_value[a$anova$term == "Treatments"], n = max(50, n_sims %/% 2))

  # Adjusted means: unbiased with honest SEs (95% t-intervals cover the truth).
  coverage <- function(label, design_data, analyze, sim) {
    sec <- "calibration"
    n_rep <- max(20, n_sims %/% 5)
    hits <- z <- numeric(0)
    for (i in seq_len(n_rep)) {
      s <- sim(design_data)
      design_data$y <- s$y
      a <- analyze(design_data)
      m <- a$adjusted_means[!is.na(a$adjusted_means$adjusted_mean), ]
      truth <- s$true_mean[m$genotype]
      q <- stats::qt(0.975, a$diagnostics$anova_error_df)
      hits <- c(hits, abs(m$adjusted_mean - truth) <= q * m$std_error)
      z <- c(z, (m$adjusted_mean - truth) / m$std_error)
    }
    cover <- mean(hits)
    lo <- if (quick) 0.92 else 0.93
    hi <- if (quick) 0.98 else 0.97
    record(sec, label, sprintf("adjusted-mean 95%% CI coverage in [%.2f, %.2f]", lo, hi), cover >= lo && cover <= hi,
           sprintf("coverage = %.3f", cover))
    record(sec, label, "adjusted means unbiased (|mean z| < 0.05)", abs(mean(z)) < 0.05, sprintf("mean z = %.3f", mean(z)))
    say("%-44s CI coverage %.3f  mean z %.3f", label, cover, mean(z))
  }
  rc_cov_sim <- function(d) {
    g <- unique(d$all_entries)
    ge <- stats::setNames(stats::rnorm(length(g), sd = 2), g)
    design_effect <- 10 + 0.5 * d$row + 0.3 * d$col + 2 * (d$block %% 2)
    list(y = design_effect + ge[d$all_entries] + stats::rnorm(nrow(d)), true_mean = mean(design_effect) + ge)
  }
  coverage("row-column 12x12, 4x4 blocks, 4 checks", rc_make(12, 12, 4, 4, 4)$design,
           function(x) analyze_augmented_row_column_design(x, "y"), rc_cov_sim)
  coverage("row-column 24x16, 6x4 blocks, 4 checks", rc_make(24, 16, 6, 4, 4)$design,
           function(x) analyze_augmented_row_column_design(x, "y"), rc_cov_sim)
  prep_cov_sim <- function(d) {
    g <- unique(d$all_entries)
    ge <- stats::setNames(stats::rnorm(length(g), sd = 2), g)
    blk <- interaction(d$location, d$block, drop = TRUE)
    be <- stats::setNames(stats::rnorm(nlevels(blk), sd = 2), levels(blk))
    design_effect <- 20 + 5 * as.integer(factor(d$location)) + be[as.character(blk)]
    list(y = design_effect + ge[d$all_entries] + stats::rnorm(nrow(d)), true_mean = mean(design_effect) + ge)
  }
  coverage("p-rep 60 entries, 2 locations", prep_design, function(x) analyze_augmented_prep_design(x, "y"), prep_cov_sim)
}

# =============================================================================
# 9. Reproducibility
# =============================================================================
if (run("repro")) {
  banner("9. Reproducibility")
  sec <- "reproducibility"
  same <- function(label, make, field = "design") {
    a <- make(11)
    b <- make(11)
    c <- make(12)
    record(sec, label, "same seed -> identical design", identical(a[[field]], b[[field]]))
    record(sec, label, "different seed -> different design", !identical(a[[field]]$all_entries, c[[field]]$all_entries))
  }
  same("row-column (latin)", function(s) rc_make(12, 12, 4, 4, 4, seed = s))
  same("row-column (optimal)", function(s) rc_make(12, 12, 4, 4, 4, seed = s, method = "optimal", n_candidates = 1))
  same("p-rep (search)", function(s) augmented_prep_design(paste0("G", 1:24), c("A", "B"), block_size = 6, n_candidates = 10, seed = s))
  same("p-rep (alpha)", function(s) augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha", block_size = 3, seed = s))
  same("p-rep (alpha2)", function(s) augmented_prep_design(paste0("G", 1:42), c("A", "B", "C"), method = "alpha2", block_size = 4, seed = s))
  same("resolvable row-column", function(s) resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5, n_starts = 1, n_iter = 200, seed = s))
  set.seed(99)
  r1 <- stats::runif(1)
  set.seed(99)
  invisible(rc_make(12, 12, 4, 4, 4, seed = 5))
  invisible(augmented_prep_design(paste0("G", 1:20), c("A", "B"), method = "alpha", block_size = 3, seed = 5))
  invisible(resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5, n_starts = 1, n_iter = 100, seed = 5))
  r2 <- stats::runif(1)
  record(sec, "all generators", "seed argument does not disturb the global RNG", r1 == r2)
}

# =============================================================================
# 10. Timing
# =============================================================================
if (run("timing")) {
  banner("10. Timing")
  time_it <- function(label, expr, limit) {
    t <- system.time(force(expr))[["elapsed"]]
    record("timing", label, sprintf("runs in < %gs", limit), t < limit, sprintf("%.2fs", t))
    say("%-60s %6.2fs", label, t)
  }
  big <- NULL
  time_it("row-column design 40x40, 1000 candidates (latin)", big <- rc_make(40, 40, 8, 8, 5, n_candidates = 1000), 10)
  big$design$y <- stats::rnorm(nrow(big$design))
  time_it("row-column analysis 40x40 (1600 plots)", analyze_augmented_row_column_design(big$design, "y"), 5)
  huge <- NULL
  time_it("row-column design 60x60, 200 candidates (latin)", huge <- rc_make(60, 60, 10, 10, 6, n_candidates = 200), 10)
  huge$design$y <- stats::rnorm(nrow(huge$design))
  time_it("row-column analysis 60x60 (3600 plots)", analyze_augmented_row_column_design(huge$design, "y"), 10)
  ex2 <- NULL
  time_it("row-column design 36x18, D-optimal search", ex2 <- rc_make(36, 18, 6, 6, 4, method = "optimal", n_candidates = 3), 60)
  if (has_lme4) {
    ex2$design$y <- stats::rnorm(nrow(ex2$design))
    time_it("row-column analysis 36x18, random blocks + entries",
            analyze_augmented_row_column_design(ex2$design, "y", random_blocks = TRUE, entries = "random"), 30)
  }
  if (has_spats) {
    time_it("row-column analysis 36x18, SpATS spatial", analyze_augmented_row_column_design(ex2$design, "y", spatial = TRUE), 30)
  }
  pbig <- NULL
  time_it("p-rep design 300 entries, 3 locs (search)",
          pbig <- augmented_prep_design(paste0("G", 1:300), c("A", "B", "C"), block_size = 12, field_cols = 20,
                                        n_candidates = 1000, seed = 1), 15)
  pbig$design$y <- stats::rnorm(nrow(pbig$design))
  time_it("p-rep analysis 300 entries, 3 locs", analyze_augmented_prep_design(pbig$design, "y"), 10)
  time_it("p-rep design 120 entries, 6 locs, blocks of 8 (search)",
          augmented_prep_design(sprintf("G%03d", 1:120), paste0("L", 1:6), block_size = 8, n_candidates = 200, seed = 1), 30)
  time_it("p-rep design 400 entries, 4 locs (alpha2)",
          augmented_prep_design(paste0("G", 1:400), paste0("L", 1:4), method = "alpha2", block_size = 10, seed = 1), 30)
  time_it("p-rep two-stage analysis 300 entries, 3 locs", analyze_augmented_prep_two_stage(pbig$design, "y"), 30)
  rrc_big <- NULL
  time_it("resolvable row-column design 100 treatments, 2 x (10x10)",
          rrc_big <- resolvable_row_column_design(paste0("T", 1:100), 2, 10, 10, n_starts = 1, seed = 1), 60)
  rrc_big$design$y <- stats::rnorm(nrow(rrc_big$design))
  time_it("resolvable row-column analysis with model comparison (200 plots)",
          analyze_resolvable_row_column(rrc_big$design, "y", model = "lv"), 60)
}

# =============================================================================
# Summary
# =============================================================================
grDevices::dev.off()
res <- do.call(rbind, results)
banner("Summary")
summary_tab <- stats::aggregate(pass ~ section, data = res, FUN = function(x) c(passed = sum(x), total = length(x)))
summary_tab <- data.frame(section = summary_tab$section, passed = summary_tab$pass[, 1], total = summary_tab$pass[, 2])
print(summary_tab, row.names = FALSE)

failed <- res[!res$pass, ]
cat(sprintf("\n%d checks, %d failed (%.0f s)\n", nrow(res), nrow(failed),
            as.numeric(difftime(Sys.time(), script_start, units = "secs"))))
if (nrow(failed) > 0) {
  cat("\nFailures (first 30):\n")
  print(utils::head(failed[, c("section", "scenario", "check", "detail")], 30), row.names = FALSE, right = FALSE)
}
if (!is.null(out_file)) {
  utils::write.csv(res, out_file, row.names = FALSE)
  cat("\nAll results written to", out_file, "\n")
}
if (!interactive()) quit(status = if (nrow(failed) > 0) 1L else 0L)
