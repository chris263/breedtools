#!/usr/bin/env Rscript
# =============================================================================
# breedtools design and analysis validation
#
# Checks every design generator and analysis function end to end:
#
#   1. Augmented row-column design
#      - layout sweep: every analyzable layout must produce a valid design and
#        every non-analyzable layout must be rejected with a clear error
#      - structural properties of each generated design
#      - analysis equals a full lm() + anova() fit (complete and missing data)
#      - ANOVA df equal the published design formulas
#   2. Augmented p-rep design
#      - structural properties over a set of scenarios
#      - analysis equals a full lm() + anova() fit; reported df are estimable
#   3. Statistical calibration (both designs)
#      - genotype F-test rejects ~5% of the time under H0 (no genotype effect)
#      - genotype F-test has high power when genotype effects exist
#   4. Reproducibility (same seed -> same design; different seeds -> different)
#   5. Timing of design generation and analysis on large trials
#
# Usage (from the package root):
#   Rscript scripts/validate_designs.R            # full run
#   Rscript scripts/validate_designs.R --quick    # smaller sweep / fewer sims
#   Rscript scripts/validate_designs.R --out results.csv
#
# Exits with status 1 if any check fails, so it can be used in CI.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
out_file <- if ("--out" %in% args) args[which(args == "--out") + 1] else NULL

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(breedtools))
}

set.seed(20261102)
n_sims <- if (quick) 100 else 400
grDevices::pdf(NULL)

# -----------------------------------------------------------------------------
# Result recorder
# -----------------------------------------------------------------------------
results <- list()

record <- function(section, scenario, check, pass, detail = "") {
  results[[length(results) + 1L]] <<- data.frame(
    section = section,
    scenario = scenario,
    check = check,
    pass = isTRUE(pass),
    detail = as.character(detail),
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

banner <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

# Reference sequential ANOVA from stats::lm(); terms in the given order.
reference_anova <- function(formula, data) {
  tab <- stats::anova(stats::lm(formula, data = data))
  data.frame(df = tab$Df, ss = tab$`Sum Sq`, p = tab$`Pr(>F)`)
}

# =============================================================================
# 1. Augmented row-column design
# =============================================================================
banner("1. Augmented row-column design")

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

rc_label <- function(R, C, rb, cb, nc) sprintf("%dx%d field, %dx%d blocks, %d checks", R, C, rb, cb, nc)

rc_make <- function(R, C, rb, cb, nc, seed = 1, n_candidates = 50) {
  cap <- augmented_row_column_capacity(R, C, rb, cb, nc)
  augmented_row_column_design(
    treatments = paste0("G", seq_len(cap$entry_plots)),
    controls = paste0("Check", seq_len(nc)),
    rows_in_field = R, cols_in_field = C,
    rows_per_block = rb, cols_per_block = cb,
    n_candidates = n_candidates, seed = seed
  )
}

# Structural checks on one row-column design. Returns TRUE if all pass.
rc_check_structure <- function(design, R, C, rb, cb, nc, scenario) {
  d <- design$design
  ck <- d[d$is_control == 1, ]
  en <- d[d$is_control == 0, ]
  n_blocks <- (R / rb) * (C / cb)
  sec <- "row-column structure"

  ok <- c(
    record(sec, scenario, "plot count = rows x cols", nrow(d) == R * C),
    record(sec, scenario, "every (row, col) used exactly once",
           !anyDuplicated(paste(d$row, d$col)) && all(d$row %in% seq_len(R)) && all(d$col %in% seq_len(C))),
    record(sec, scenario, "plot numbers are 1..n", setequal(d$plots, seq_len(R * C))),
    record(sec, scenario, "block id matches row/col groups",
           all(d$rowgroup == (d$row - 1) %/% rb + 1) &&
             all(d$colgroup == (d$col - 1) %/% cb + 1) &&
             all(d$block == (d$rowgroup - 1) * (C / cb) + d$colgroup)),
    record(sec, scenario, "each check exactly once per block",
           nrow(ck) == n_blocks * nc && all(table(ck$all_entries, ck$block) == 1)),
    record(sec, scenario, "checks in a block use distinct rows and cols",
           all(vapply(split(ck, ck$block), function(b) !anyDuplicated(b$row) && !anyDuplicated(b$col), logical(1)))),
    record(sec, scenario, "same check never twice in a field row",
           !any(table(ck$all_entries, ck$row) > 1)),
    record(sec, scenario, "same check never twice in a field column",
           !any(table(ck$all_entries, ck$col) > 1)),
    record(sec, scenario, "every field row has a check", setequal(unique(ck$row), seq_len(R))),
    record(sec, scenario, "every field column has a check", setequal(unique(ck$col), seq_len(C))),
    record(sec, scenario, "each entry appears exactly once",
           !anyDuplicated(en$all_entries) &&
             setequal(en$all_entries, design$input$treatments)),
    record(sec, scenario, "flags consistent (type / is_control / rep)",
           all((d$type == "check") == (d$is_control == 1)) &&
             all(d$rep[d$is_control == 1] == d$block[d$is_control == 1]) &&
             all(d$rep[d$is_control == 0] == 1)),
    record(sec, scenario, "breedbase_design matches design",
           identical(design$breedbase_design$all_entries, d$all_entries) &&
             identical(design$breedbase_design$plots, d$plots))
  )

  # Serpentine order: consecutive plots are physically adjacent.
  o <- d[order(d$plots), ]
  steps <- abs(diff(o$row)) + abs(diff(o$col))
  ok <- c(ok, record(sec, scenario, "serpentine plots are adjacent", all(steps == 1)))

  all(ok)
}

# Analysis checks on one row-column design.
rc_check_analysis <- function(design, scenario, error_df_expected) {
  sec <- "row-column analysis"
  d <- design$design
  # Simulated data with real row, column, block, and genotype effects.
  geno_eff <- stats::setNames(stats::rnorm(length(unique(d$all_entries)), sd = 2), unique(d$all_entries))
  d$y <- 50 + 0.4 * d$row - 0.2 * d$col + d$block %% 3 + geno_eff[d$all_entries] + stats::rnorm(nrow(d))

  formula <- y ~ rg + rn + cg + cn + bk + g
  to_ref <- function(x) {
    data.frame(
      y = x$y,
      rg = factor(x$rowgroup),
      rn = interaction(x$rowgroup, x$row, drop = TRUE),
      cg = factor(x$colgroup),
      cn = interaction(x$colgroup, x$col, drop = TRUE),
      bk = interaction(x$rowgroup, x$colgroup, drop = TRUE),
      g = factor(x$all_entries)
    )
  }

  for (missing in c(FALSE, TRUE)) {
    x <- d
    tag <- if (missing) "missing plots" else "complete data"
    if (missing) {
      x$y[sample(nrow(x), max(3, round(0.03 * nrow(x))))] <- NA
    }

    a <- safely(sec, scenario, paste("analysis runs,", tag), analyze_augmented_row_column_design(x, "y"))
    if (is.null(a)) next
    ref <- reference_anova(formula, to_ref(x[!is.na(x$y), ]))

    record(sec, scenario, paste("df equal lm(),", tag), isTRUE(all.equal(a$anova$df, ref$df)),
           paste(a$anova$df, collapse = ","))
    record(sec, scenario, paste("sums of squares equal lm(),", tag),
           isTRUE(all.equal(a$anova$sum_sq, ref$ss, tolerance = 1e-8)))
    record(sec, scenario, paste("p-values equal lm(),", tag),
           isTRUE(all.equal(a$anova$p_value[1:6], ref$p[1:6], tolerance = 1e-6)))
    record(sec, scenario, paste("df sum to corrected total,", tag),
           sum(a$anova_formatted$Df[1:7]) == a$anova_formatted$Df[8] &&
             a$anova_formatted$Df[8] == sum(!is.na(x$y)) - 1)
    record(sec, scenario, paste("residuals = y - fitted,", tag),
           isTRUE(all.equal(a$plot_results$residual, x$y - a$plot_results$fitted)))
    compare_adjusted(sec, scenario, tag, a, reference_adjusted(formula, to_ref(x[!is.na(x$y), ])))

    if (!missing) {
      record(sec, scenario, "error df = capacity error_df",
             a$diagnostics$anova_error_df == error_df_expected,
             paste(a$diagnostics$anova_error_df, "vs", error_df_expected))
      record(sec, scenario, "design-formula df = fitted df",
             isTRUE(all.equal(a$anova_df$display_df, a$anova$df)))
      record(sec, scenario, "no aliased design effects", !a$diagnostics$rank_deficient)
    }
  }
}

# --- 1a. Layout sweep --------------------------------------------------------
cat("Sweeping row-column layouts...\n")
dims_rows <- if (quick) c(8, 12, 16, 24) else c(8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 36)
block_sizes <- if (quick) c(3, 4, 6, 8) else c(3, 4, 5, 6, 8, 10)
checks <- if (quick) 2:5 else 2:6

sweep <- expand.grid(R = dims_rows, C = dims_rows, rb = block_sizes, cb = block_sizes, nc = checks)
sweep <- sweep[sweep$R %% sweep$rb == 0 & sweep$C %% sweep$cb == 0, ]

n_valid <- 0L
n_rejected <- 0L
sweep_time <- system.time({
  for (i in seq_len(nrow(sweep))) {
    s <- sweep[i, ]
    scen <- rc_label(s$R, s$C, s$rb, s$cb, s$nc)
    cap <- tryCatch(augmented_row_column_capacity(s$R, s$C, s$rb, s$cb, s$nc), error = function(e) NULL)
    if (is.null(cap) || cap$entry_plots < 1) next

    # Is the layout analyzable? (all conditions documented in the generator)
    nbr <- s$R / s$rb
    nbc <- s$C / s$cb
    feasible <- cap$error_df > 0 &&
      nbc <= s$rb && nbr <= s$cb &&
      s$nc * nbc >= s$rb && s$nc * nbr >= s$cb

    design <- tryCatch(rc_make(s$R, s$C, s$rb, s$cb, s$nc, seed = i, n_candidates = 20),
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

    n_valid <- n_valid + 1L
    record("row-column sweep", scen, "feasible layout is generated", TRUE)
    rc_check_structure(design, s$R, s$C, s$rb, s$cb, s$nc, scen)
    # Full analysis comparison on every 5th valid layout to keep runtime low.
    if (n_valid %% 5 == 1) rc_check_analysis(design, scen, cap$error_df)
  }
})[["elapsed"]]
cat(sprintf("  %d valid layouts, %d infeasible layouts rejected (%.1fs)\n", n_valid, n_rejected, sweep_time))

# --- 1b. Named scenarios, always analyzed ------------------------------------
named_rc <- list(
  c(12, 12, 4, 4, 4),
  c(36, 18, 6, 6, 4),
  c(24, 16, 6, 4, 4),
  c(40, 40, 8, 8, 5)
)
for (cfg in named_rc) {
  scen <- do.call(rc_label, as.list(cfg))
  design <- safely("row-column named", scen, "design generated", do.call(rc_make, as.list(cfg)))
  if (is.null(design)) next
  cap <- do.call(augmented_row_column_capacity, as.list(unname(cfg)))
  rc_check_structure(design, cfg[1], cfg[2], cfg[3], cfg[4], cfg[5], scen)
  rc_check_analysis(design, scen, cap$error_df)
}

# --- 1c. Input validation ----------------------------------------------------
sec <- "row-column validation"
expect_error_msg <- function(scenario, expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  record(sec, scenario, paste("errors with", shQuote(pattern)), !is.na(msg) && grepl(pattern, msg), msg)
}
expect_error_msg("wrong treatment count",
  augmented_row_column_design(paste0("G", 1:100), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
  "Expected 108 treatments")
expect_error_msg("duplicated treatments",
  augmented_row_column_design(c(paste0("G", 1:107), "G1"), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
  "duplicated")
expect_error_msg("treatment also a control",
  augmented_row_column_design(c(paste0("G", 1:107), "C1"), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2),
  "must not contain the same name")
expect_error_msg("field not divisible by block",
  augmented_row_column_capacity(10, 12, 4, 4, 4), "divisible")
expect_error_msg("no error df",
  augmented_row_column_design(paste0("G", 1:28), paste0("C", 1:2), 6, 6, 3, 3, n_candidates = 2),
  "no error degrees of freedom")

# --- 1d. Plot function -------------------------------------------------------
design <- rc_make(12, 12, 4, 4, 4)
for (fill in c("type", "block", "rowgroup", "colgroup", "is_control")) {
  p <- safely("row-column plot", fill, "plot draws", plot_augmented_row_column_design(design, fill = fill, show_block_labels = TRUE))
  if (!is.null(p)) record("row-column plot", paste("fill =", fill), "plot draws", nrow(p) == 144)
}

# =============================================================================
# 2. Augmented p-rep design
# =============================================================================
banner("2. Augmented p-rep design")

prep_scenarios <- list(
  list(name = "12 entries, 3 locs, block 5",
       args = list(treatments = paste0("G", 1:12), locations = c("Loc1", "Loc2", "Loc3"),
                   block_size = 5, field_cols = 5)),
  list(name = "12 entries, uneven per-location settings",
       args = list(treatments = paste0("G", 1:12), locations = c("Loc1", "Loc2", "Loc3"),
                   block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 6),
                   n_blocks_per_location = c(Loc1 = 4, Loc2 = 4, Loc3 = 4),
                   n_duplicate_per_location = c(Loc1 = 4, Loc2 = 2, Loc3 = 4),
                   replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3, Loc3 = 2),
                   field_cols = 8)),
  list(name = "60 entries, 2 locs, prep_prop 0.3",
       args = list(treatments = paste0("G", 1:60), locations = c("A", "B"),
                   block_size = 8, prep_prop = 0.3, field_cols = 10)),
  list(name = "100 entries, 4 locs, default split",
       args = list(treatments = paste0("G", 1:100), locations = paste0("L", 1:4),
                   block_size = 10, field_cols = 10)),
  list(name = "48 entries, 2 locs, block layout",
       args = list(treatments = paste0("G", 1:48), locations = c("A", "B"),
                   block_size = 6, n_duplicate_per_location = 12,
                   rows_per_block = 2, blocks_per_row = 2)),
  list(name = "200 entries, 1 loc, 25% replicated",
       args = list(treatments = paste0("G", 1:200), locations = "Station",
                   block_size = 10, prep_prop = 0.25, field_cols = 20))
)

prep_check <- function(sc) {
  sec <- "p-rep structure"
  scen <- sc$name
  design <- safely(sec, scen, "design generated",
                   do.call(augmented_prep_design, c(sc$args, list(n_candidates = if (quick) 20 else 100, seed = 7))))
  if (is.null(design)) return(invisible(NULL))
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
    expected <- ifelse(names(counts) %in% repeated, reps[[loc]], 1L)
    rep_ok <- rep_ok && all(counts == expected)
  }
  record(sec, scen, "replication per location as requested", rep_ok)
  record(sec, scen, "no entry twice in the same block",
         !anyDuplicated(paste(d$location, d$block, d$all_entries)))

  sizes <- table(paste(d$location, d$block))
  by_loc <- split(as.integer(sizes), sub(" .*", "", names(sizes)))
  record(sec, scen, "block sizes within a location differ by <= 1",
         all(vapply(by_loc, function(s) max(s) - min(s) <= 1, logical(1))))
  if (is.null(inp$n_blocks_per_location)) {
    record(sec, scen, "block sizes <= block_size",
           all(vapply(names(by_loc), function(l) max(by_loc[[l]]) <= inp$block_size[[l]], logical(1))))
  }
  record(sec, scen, "field coordinates unique per location",
         !anyDuplicated(paste(d$location, d$row, d$col)))
  record(sec, scen, "plot numbers unique", !anyDuplicated(d$plots))
  conc <- augmented_prep_concurrence(d, treatments = rownames(design$concurrence_matrix))
  record(sec, scen, "concurrence matrix consistent", identical(unname(conc), unname(design$concurrence_matrix)))

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
                           blk = interaction(xr$location, xr$block, drop = TRUE),
                           g = factor(xr$all_entries))
    ref <- if (length(unique(ref_data$loc)) > 1) {
      reference_anova(y ~ loc + blk + g, ref_data)
    } else {
      rbind(data.frame(df = 0, ss = 0, p = NA), reference_anova(y ~ blk + g, ref_data))
    }
    fitted_df <- a$anova$df[match(c("location", "block_nested", "genotype", "Residuals"), a$anova$term)]
    fitted_df[is.na(fitted_df)] <- 0
    record(sec, scen, paste("fitted df equal lm(),", tag), isTRUE(all.equal(fitted_df, ref$df)),
           paste(fitted_df, collapse = ","))
    ss <- a$anova$sum_sq[match(c("location", "block_nested", "genotype", "Residuals"), a$anova$term)]
    ss[is.na(ss)] <- 0
    record(sec, scen, paste("sums of squares equal lm(),", tag),
           isTRUE(all.equal(ss, ref$ss, tolerance = 1e-8)))
    ref_formula <- if (length(unique(ref_data$loc)) > 1) y ~ loc + blk + g else y ~ blk + g
    compare_adjusted(sec, scen, tag, a, reference_adjusted(ref_formula, ref_data))
    # The formatted table uses design-formula df; they must be estimable.
    record(sec, scen, paste("reported df = estimable df,", tag),
           isTRUE(all.equal(a$anova_formatted$Df[1:4], ref$df)),
           paste("reported", paste(a$anova_formatted$Df[1:4], collapse = ","),
                 "| estimable", paste(ref$df, collapse = ",")))
  }

  p <- safely("p-rep plot", scen, "plot draws", plot_augmented_prep_design(design))
  if (!is.null(p)) record("p-rep plot", scen, "plot draws", TRUE)
}

for (sc in prep_scenarios) prep_check(sc)

# =============================================================================
# 3. Statistical calibration
# =============================================================================
banner("3. Statistical calibration")
# Under H0 (no genotype differences) the genotype F-test p-values should be
# Uniform(0,1): rejection rate at 5% ~ 0.05. Wrong df show up here.
# Binomial 99.9% interval for n_sims tests is used as the pass band.

calibrate <- function(label, design_data, analyze, sim_y, genotype_term = "genotype") {
  sec <- "calibration"
  band <- stats::qbinom(c(0.0005, 0.9995), n_sims, 0.05) / n_sims

  p_null <- p_alt <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    design_data$y <- sim_y(design_data, genotype_sd = 0)
    p_null[i] <- with(analyze(design_data), anova$p_value[anova$term == genotype_term])
    design_data$y <- sim_y(design_data, genotype_sd = 3)
    p_alt[i] <- with(analyze(design_data), anova$p_value[anova$term == genotype_term])
  }
  size <- mean(p_null < 0.05)
  power <- mean(p_alt < 0.05)
  ks <- stats::ks.test(p_null, "punif")$p.value
  record(sec, label, "type I error ~ 5% under H0", size >= band[1] && size <= band[2],
         sprintf("rate = %.3f (band %.3f-%.3f)", size, band[1], band[2]))
  record(sec, label, "null p-values uniform (KS p > 0.001)", ks > 0.001, sprintf("KS p = %.3g", ks))
  record(sec, label, "power > 0.9 with genotype sd = 3 x error sd", power > 0.9, sprintf("power = %.3f", power))
  cat(sprintf("  %-40s size %.3f  power %.3f  KS p %.3g\n", label, size, power, ks))
}

rc_sim <- function(d, genotype_sd) {
  g <- unique(d$all_entries)
  ge <- stats::setNames(stats::rnorm(length(g), sd = genotype_sd), g)
  10 + 0.5 * d$row + 0.3 * d$col + 2 * (d$block %% 2) + ge[d$all_entries] + stats::rnorm(nrow(d))
}
calibrate("row-column 12x12, 4x4 blocks, 4 checks",
          rc_make(12, 12, 4, 4, 4)$design,
          function(x) analyze_augmented_row_column_design(x, "y"), rc_sim)
calibrate("row-column 24x16, 6x4 blocks, 4 checks",
          rc_make(24, 16, 6, 4, 4)$design,
          function(x) analyze_augmented_row_column_design(x, "y"), rc_sim)

prep_sim <- function(d, genotype_sd) {
  g <- unique(d$all_entries)
  ge <- stats::setNames(stats::rnorm(length(g), sd = genotype_sd), g)
  blk <- interaction(d$location, d$block, drop = TRUE)
  be <- stats::setNames(stats::rnorm(nlevels(blk), sd = 2), levels(blk))
  20 + 5 * as.integer(factor(d$location)) + be[as.character(blk)] + ge[d$all_entries] + stats::rnorm(nrow(d))
}
prep_design <- augmented_prep_design(paste0("G", 1:60), c("A", "B"), block_size = 8, prep_prop = 0.3,
                                     field_cols = 10, n_candidates = 20, seed = 3)$design
calibrate("p-rep 60 entries, 2 locations",
          prep_design, function(x) analyze_augmented_prep_design(x, "y"), prep_sim)

# Adjusted means: unbiased with honest SEs. With known effects the true
# adjusted mean of a genotype is its effect plus the average design effect,
# so 95% t-intervals should cover it ~95% of the time.
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
  record(sec, label, "adjusted-mean 95% CI coverage in [0.93, 0.97]", cover >= 0.93 && cover <= 0.97,
         sprintf("coverage = %.3f", cover))
  record(sec, label, "adjusted means unbiased (|mean z| < 0.05)", abs(mean(z)) < 0.05,
         sprintf("mean z = %.3f", mean(z)))
  cat(sprintf("  %-40s CI coverage %.3f  mean z %.3f\n", label, cover, mean(z)))
}

rc_cov_sim <- function(d) {
  g <- unique(d$all_entries)
  ge <- stats::setNames(stats::rnorm(length(g), sd = 2), g)
  design_effect <- 10 + 0.5 * d$row + 0.3 * d$col + 2 * (d$block %% 2)
  list(y = design_effect + ge[d$all_entries] + stats::rnorm(nrow(d)),
       true_mean = mean(design_effect) + ge)
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
  list(y = design_effect + ge[d$all_entries] + stats::rnorm(nrow(d)),
       true_mean = mean(design_effect) + ge)
}
coverage("p-rep 60 entries, 2 locations", prep_design,
         function(x) analyze_augmented_prep_design(x, "y"), prep_cov_sim)

# =============================================================================
# 4. Reproducibility
# =============================================================================
banner("4. Reproducibility")
sec <- "reproducibility"
a <- rc_make(12, 12, 4, 4, 4, seed = 11)$design
b <- rc_make(12, 12, 4, 4, 4, seed = 11)$design
c <- rc_make(12, 12, 4, 4, 4, seed = 12)$design
record(sec, "row-column", "same seed -> identical design", identical(a, b))
record(sec, "row-column", "different seed -> different design", !identical(a$all_entries, c$all_entries))
before <- stats::runif(1); set.seed(99); r1 <- stats::runif(1)
set.seed(99); invisible(rc_make(12, 12, 4, 4, 4, seed = 5)); r2 <- stats::runif(1)
record(sec, "row-column", "seed argument does not disturb global RNG", r1 == r2)

pa <- augmented_prep_design(paste0("G", 1:24), c("A", "B"), block_size = 6, n_candidates = 10, seed = 4)$design
pb <- augmented_prep_design(paste0("G", 1:24), c("A", "B"), block_size = 6, n_candidates = 10, seed = 4)$design
record(sec, "p-rep", "same seed -> identical design", identical(pa, pb))

# =============================================================================
# 5. Timing
# =============================================================================
banner("5. Timing")
time_it <- function(label, expr, limit) {
  t <- system.time(force(expr))[["elapsed"]]
  record("timing", label, sprintf("runs in < %gs", limit), t < limit, sprintf("%.2fs", t))
  cat(sprintf("  %-55s %6.2fs\n", label, t))
}
big <- NULL
time_it("row-column design 40x40, 1000 candidates",
        big <- rc_make(40, 40, 8, 8, 5, n_candidates = 1000), 10)
big$design$y <- stats::rnorm(nrow(big$design))
time_it("row-column analysis 40x40 (1600 plots)",
        analyze_augmented_row_column_design(big$design, "y"), 5)
huge <- NULL
time_it("row-column design 60x60, 200 candidates",
        huge <- rc_make(60, 60, 10, 10, 6, n_candidates = 200), 10)
huge$design$y <- stats::rnorm(nrow(huge$design))
time_it("row-column analysis 60x60 (3600 plots)",
        analyze_augmented_row_column_design(huge$design, "y"), 10)
pbig <- NULL
time_it("p-rep design 300 entries, 3 locs, 1000 candidates",
        pbig <- augmented_prep_design(paste0("G", 1:300), c("A", "B", "C"), block_size = 12,
                                      field_cols = 20, n_candidates = 1000, seed = 1), 10)
pbig$design$y <- stats::rnorm(nrow(pbig$design))
time_it("p-rep design 30 entries x 4 locs (no early stop), 1000 cand.",
        augmented_prep_design(paste0("G", 1:30), paste0("L", 1:4), block_size = 10, prep_prop = 0.2,
                              n_candidates = 1000, seed = 1), 10)
time_it("p-rep analysis 300 entries, 3 locs",
        analyze_augmented_prep_design(pbig$design, "y"), 10)

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
cat(sprintf("\n%d checks, %d failed\n", nrow(res), nrow(failed)))
if (nrow(failed) > 0) {
  cat("\nFailures (first 30):\n")
  print(utils::head(failed[, c("section", "scenario", "check", "detail")], 30), row.names = FALSE, right = FALSE)
}

if (!is.null(out_file)) {
  utils::write.csv(res, out_file, row.names = FALSE)
  cat("\nAll results written to", out_file, "\n")
}

if (!interactive()) quit(status = if (nrow(failed) > 0) 1L else 0L)
