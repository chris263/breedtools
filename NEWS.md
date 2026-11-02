# breedtools 0.1.0

* First CRAN release.
* `augmented_row_column_design()`:
  * Check placement now uses covering Latin rectangles, so every field row
    and column gets a check.
  * Layouts whose effects cannot be estimated, or that leave no error
    degrees of freedom, are rejected.
  * The design is much faster to generate.
* `augmented_row_column_capacity()` reports the error degrees of freedom of
  a layout.
* `analyze_augmented_row_column_design()`:
  * Degrees of freedom are corrected (rows within row groups `k - gk`,
    columns within column groups `s - gs`), and the ANOVA uses the df that
    are estimable from the observed plots.
  * The model is fitted by absorbing unreplicated entries, so it is fast for
    large fields.
* `augmented_prep_design()` uses a faster candidate search with incremental
  concurrence scoring. The search stops early once an optimal design is
  found.
* `analyze_augmented_prep_design()`:
  * Handles single-location trials.
  * Reports df that are estimable when plots are missing.
* Both analyses now return adjusted entry means (BLUEs) with standard errors
  (`$adjusted_means`), and report the mean SED and the 5% LSD in
  `$diagnostics`.
* New vignette, `vignette("breedtools")`, and an RStudio workflow script in
  `inst/examples`.
