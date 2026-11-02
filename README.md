# breedtools

Design and analysis of augmented plant breeding trials.

- **Augmented row-column designs**: unreplicated entries plus replicated
  checks. Each check appears once per block, and every field row and column
  contains a check.
- **Augmented p-rep designs**: multi-location, partially replicated trials
  in incomplete blocks with low pairwise concurrence.
- **Analysis**: ANOVA with degrees of freedom that are estimable from the
  observed plots, plus **adjusted entry means (BLUEs)** with standard
  errors, the mean SED and the 5% LSD.
- **Field maps and 'Breedbase' export**: a plot-level field book and upload
  files.
- **Multi-environment tools**: read 'Breedbase' phenotype files, run
  Hartley's Fmax variance check, and group environments by GxE.

## Installation

```r
# install.packages("remotes")
remotes::install_github("chris263/breedtools")
```

## Example

```r
library(breedtools)

rc <- augmented_row_column_design(
  treatments = sprintf("LINE_%03d", 1:108),
  controls = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"),
  rows_in_field = 12, cols_in_field = 12,
  rows_per_block = 4, cols_per_block = 4,
  seed = 2026
)
plot_augmented_row_column_design(rc)

field <- rc$design
field$yield <- rnorm(nrow(field), 60, 3)   # your phenotypes here
fit <- analyze_augmented_row_column_design(field, response = "yield")
fit                      # ANOVA and top adjusted means
fit$adjusted_means       # all entries ranked by adjusted mean
```

To see the full workflow, run `vignette("breedtools")`, or open the RStudio
script:

```r
file.edit(system.file("examples", "breedtools_workflow.R", package = "breedtools"))
```
