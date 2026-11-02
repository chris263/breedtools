# =============================================================================
# breedtools: example workflow for trial design and analysis
#
# Open this file in RStudio with:
#   file.edit(system.file("examples", "breedtools_workflow.R", package = "breedtools"))
#
# Run it section by section (Ctrl/Cmd + Enter, or Ctrl/Cmd + Alt + T to run
# the current section). Use the outline pane (Ctrl/Cmd + Shift + O) to jump
# between sections.
#
#   Part A - Augmented row-column design (unreplicated entries + checks)
#   Part B - Augmented p-rep design (partially replicated entries, multi-location)
#
# The phenotypes in this script are SIMULATED so the example runs anywhere.
# Replace those steps with your own field data.
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install once (pick the option that applies to you):
# install.packages("remotes")
# install.packages("breedtools")                    # from CRAN (once released)
# remotes::install_github("chris263/breedtools")    # development version
# devtools::install("path/to/breedtools")           # from a local copy

library(breedtools)

# Folder for files this script writes (design sheets, phenotype files).
# Change it to a project folder of your own, e.g. "~/trials/2026".
out_dir <- file.path(tempdir(), "breedtools_example")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir


# =============================================================================
# PART A - AUGMENTED ROW-COLUMN DESIGN
#
# Use when you have many unreplicated entries (e.g. early generation lines)
# and a few replicated checks. The field is split into rectangular blocks;
# every block gets one plot of each check, and row/column effects are
# estimated from the checks.
# =============================================================================


# ---- A1. Plan the layout: will it be analyzable? ----------------------------

# A row-column layout needs enough check plots to estimate row, column, and
# block effects AND leave error degrees of freedom. Check before you plant.
augmented_row_column_capacity(
  rows_in_field  = 12,
  cols_in_field  = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_controls     = 4
)
# entry_plots = how many test entries fit; error_df must be > 0.

# Compare a few candidate layouts for a 12 x 12 field.
layouts <- expand.grid(
  rows_per_block = c(3, 4, 6),
  cols_per_block = c(3, 4, 6),
  n_controls     = c(2, 3, 4, 5)
)
layout_table <- do.call(rbind, lapply(seq_len(nrow(layouts)), function(i) {
  tryCatch(
    augmented_row_column_capacity(12, 12, layouts$rows_per_block[i],
                                  layouts$cols_per_block[i], layouts$n_controls[i]),
    error = function(e) NULL
  )
}))
# Analyzable layouts, most entry plots first. There is a trade-off: more
# check plots mean fewer entries but more error df, i.e. more reliable F-tests
# and spatial adjustment. Aim for at least ~10 error df when you can.
layout_table <- layout_table[layout_table$error_df > 0, ]
layout_table[order(-layout_table$entry_plots),
             c("rows_per_block", "cols_per_block", "n_controls",
               "n_blocks", "entry_plots", "check_plots", "error_df")]


# ---- A2. Generate the design ------------------------------------------------

# Your entry and check names (usually read from a seed list / spreadsheet).
entries <- sprintf("LINE_%03d", 1:108)
checks  <- c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D")

rc <- augmented_row_column_design(
  treatments     = entries,
  controls       = checks,
  rows_in_field  = 12,
  cols_in_field  = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  plot_type      = "serpentine",   # plot numbering; or "cartesian"
  n_candidates   = 1000,           # candidate check layouts to compare
  seed           = 2026            # keep the seed to reproduce the design
)

rc$summary          # size, error df, design score
head(rc$design)     # plot-level field book


# ---- A3. Inspect the design --------------------------------------------------

# One plot of each check per block:
table(rc$design$all_entries[rc$design$is_control == 1],
      rc$design$block[rc$design$is_control == 1])

# Every field row and column contains a check:
check_plots <- rc$design[rc$design$is_control == 1, ]
all(1:12 %in% check_plots$row) && all(1:12 %in% check_plots$col)


# ---- A4. Plot the field map --------------------------------------------------

plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
plot_augmented_row_column_design(rc, label = "none", fill = "block",
                                 show_block_labels = TRUE)

# Save the map as a PDF to take to the field:
pdf(file.path(out_dir, "row_column_field_map.pdf"), width = 11, height = 9)
plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
dev.off()


# ---- A5. Export the design ---------------------------------------------------

# Full field book (for data collection):
write.csv(rc$design, file.path(out_dir, "row_column_fieldbook.csv"), row.names = FALSE)

# Breedbase upload columns (plots, block, all_entries, rep, is_control):
write.csv(rc$breedbase_design, file.path(out_dir, "row_column_breedbase.csv"),
          row.names = FALSE)

# Tip: augmented_row_column_design(..., output_file = "trial.design") writes a
# tab-delimited copy directly.


# ---- A6. After harvest: bring in phenotypes ----------------------------------

# In practice you read the file you collected, e.g.
#   pheno <- read_breedbase_file("my_trial_phenotypes.xlsx")
# Here we SIMULATE a yield file with a fertility gradient down the rows and
# a few missing plots, write it, and read it back the same way.
set.seed(1)
sim <- rc$design[, c("plots", "all_entries")]
sim$yield <- round(
  60 + 0.8 * rc$design$row - 0.3 * rc$design$col +       # field trend
    ifelse(rc$design$is_control == 1, 5, 0) +            # checks yield more
    rnorm(nrow(sim), sd = 3), 1
)
sim$yield[c(7, 50, 101)] <- NA                           # lost plots
write.csv(sim, file.path(out_dir, "row_column_phenotypes.csv"), row.names = FALSE)

pheno <- read_breedbase_file(file.path(out_dir, "row_column_phenotypes.csv"))

# Join phenotypes to the design by plot number:
rc_data <- merge(rc$design, pheno[, c("plots", "yield")], by = "plots")
head(rc_data)


# ---- A7. Analyze --------------------------------------------------------------

rc_fit <- analyze_augmented_row_column_design(rc_data, response = "yield")

rc_fit                       # ANOVA + top adjusted means (notes if not estimable)
rc_fit$anova_formatted       # ANOVA as a data frame
rc_fit$diagnostics           # error df, sigma, R^2, mean SED, LSD, warnings

# ---- A8. Adjusted entry means: rank and select -------------------------------

# Adjusted means (BLUEs) remove the row, column, and block effects of the plot
# each entry landed in, estimated from the checks. Rank on these, not on the
# raw plot values: an entry in a poor corner of the field is not penalized.
rc_means <- rc_fit$adjusted_means
head(rc_means, 10)

# How much did spatial adjustment reorder the entries?
entries_only <- rc_means[rc_means$is_control == 0 & !is.na(rc_means$adjusted_mean), ]
cor(entries_only$raw_mean, entries_only$adjusted_mean, method = "spearman")

# Entries that beat the best check by more than one LSD (5%):
best_check <- max(rc_means$adjusted_mean[rc_means$is_control == 1], na.rm = TRUE)
lsd <- rc_fit$diagnostics$lsd_5pct
entries_only[entries_only$adjusted_mean > best_check + lsd,
             c("rank", "genotype", "adjusted_mean", "std_error")]

# Residual checks:
op <- par(mfrow = c(1, 2))
plot(rc_fit$plot_results$fitted, rc_fit$plot_results$residual,
     xlab = "Fitted", ylab = "Residual", main = "Residuals vs fitted")
abline(h = 0, lty = 2)
qqnorm(rc_fit$plot_results$residual[rc_fit$plot_results$is_control == 1],
       main = "Q-Q plot (check plots)")
qqline(rc_fit$plot_results$residual[rc_fit$plot_results$is_control == 1])
par(op)

# Save results:
write.csv(rc_fit$anova_formatted, file.path(out_dir, "row_column_anova.csv"), row.names = FALSE)
write.csv(rc_means, file.path(out_dir, "row_column_adjusted_means.csv"), row.names = FALSE)
write.csv(rc_fit$plot_results, file.path(out_dir, "row_column_plot_results.csv"), row.names = FALSE)


# =============================================================================
# PART B - AUGMENTED P-REP DESIGN
#
# Use for multi-location trials without replicated checks: every entry is
# planted once at every location, and a subset of entries gets extra plots
# at one location. Plots are grouped into incomplete blocks so that pairs of
# entries rarely share a block more than once.
# =============================================================================


# ---- B1. Generate the design --------------------------------------------------

prep_entries <- sprintf("G%03d", 1:120)
sites        <- c("North", "Central", "South")

prep <- augmented_prep_design(
  treatments      = prep_entries,
  locations       = sites,
  block_size      = 10,          # max plots per incomplete block
  prep_prop       = 0.2,         # 20% of entries repeated at each location
  replications_per_repeated_entry = 2,
  field_cols      = 10,          # field columns, for row/col coordinates
  serpentine      = TRUE,
  n_candidates    = 1000,        # stops early if an optimal design is found
  seed            = 2026
)

prep$summary                 # concurrence diagnostics, plots, blocks
head(prep$design)
prep$repeated_groups         # which entries are repeated where


# ---- B2. Inspect the design ---------------------------------------------------

# Plots per entry (one per location + extra plots where repeated):
table(table(prep$design$all_entries))

# Block sizes by location:
table(prep$design$location, prep$design$block)

# Pairwise concurrence: how often two entries share a block (ideally <= 1)
conc <- prep$concurrence_matrix
table(conc[upper.tri(conc)])


# ---- B3. Plot and export ------------------------------------------------------

plot_augmented_prep_design(prep, locations = "North", fill = "entry_type")
plot_augmented_prep_design(prep, locations = "South", label = "block", fill = "block")

pdf(file.path(out_dir, "prep_field_maps.pdf"), width = 11, height = 9)
plot_augmented_prep_design(prep)     # one page per location
dev.off()

write.csv(prep$design, file.path(out_dir, "prep_fieldbook.csv"), row.names = FALSE)
write.csv(prep$breedbase_design, file.path(out_dir, "prep_breedbase.csv"), row.names = FALSE)


# ---- B4. After harvest: phenotypes -------------------------------------------

# SIMULATED: location and block effects + genotype differences + noise.
set.seed(2)
genotype_effect <- setNames(rnorm(length(prep_entries), sd = 4), prep_entries)
location_effect <- c(North = 0, Central = 8, South = -5)
block_id        <- interaction(prep$design$location, prep$design$block)
block_effect    <- setNames(rnorm(nlevels(block_id), sd = 2), levels(block_id))

prep_data <- prep$design
prep_data$yield <- 50 +
  location_effect[prep_data$location] +
  block_effect[as.character(block_id)] +
  genotype_effect[prep_data$all_entries] +
  rnorm(nrow(prep_data), sd = 3)
prep_data$yield[sample(nrow(prep_data), 6)] <- NA     # a few lost plots


# ---- B5. Analyze --------------------------------------------------------------

prep_fit <- analyze_augmented_prep_design(prep_data, response = "yield")

prep_fit                     # ANOVA + top adjusted means
prep_fit$diagnostics

# Adjusted means (BLUEs) correct for the location and block of every plot,
# so entries are compared fairly even though only some are repeated.
prep_means <- prep_fit$adjusted_means
head(prep_means, 10)

# With simulated data we know the truth: adjusted means track it better.
cor(prep_means$raw_mean, genotype_effect[prep_means$genotype], use = "complete.obs")
cor(prep_means$adjusted_mean, genotype_effect[prep_means$genotype], use = "complete.obs")

# Select the top 10% of entries:
top <- prep_means[which(prep_means$rank <= ceiling(0.1 * nrow(prep_means))), ]
top[, c("rank", "genotype", "n_observations", "adjusted_mean", "std_error")]

# A single location can be analyzed on its own:
north_fit <- analyze_augmented_prep_design(
  prep_data[prep_data$location == "North", ], response = "yield"
)
north_fit$anova_formatted

write.csv(prep_fit$anova_formatted, file.path(out_dir, "prep_anova.csv"), row.names = FALSE)
write.csv(prep_means, file.path(out_dir, "prep_adjusted_means.csv"), row.names = FALSE)


# ---- Done ---------------------------------------------------------------------

list.files(out_dir)
# browseURL(out_dir)   # open the output folder
