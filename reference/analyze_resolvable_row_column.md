# Analyze a resolvable row-column trial with RC, LV + RC or AR + RC models

Fits the variance models of Williams, John and Whitaker (2006) by
residual maximum likelihood (REML), with replicates and treatments
fixed: \$\$E(y) = R\pi + X\tau,\quad var(y) = \sigma^2 V.\$\$

- `model = "rc"`:

  The standard resolvable row-column mixed model: \\V = I + \eta_R P_R +
  \eta_C P_C\\, i.e. random rows and columns within replicates.

- `model = "lv"`:

  The two-dimensional linear variance model (LV + RC, Eq. 2): \\V = I +
  \eta_R P_R + \eta_C P_C - (\phi_R F_R + \phi_C F_C + \phi\_{RC}
  F\_{RC})\\, with LV matrices for plots in the same row, in the same
  column, and in neither. `lv_terms` chooses which \\\phi\\ terms are
  fitted and `equal_phi = TRUE` uses one \\\phi\\ for rows and columns
  (square plots).

- `model = "ar"`:

  Rows and columns plus a separable first-order autoregressive process
  (AR + RC, Eq. 4): \\V = I + \eta_R P_R + \eta_C P_C + \psi D(\rho_C)
  \otimes D(\rho_R)\\ within replicates.

In every model the identity term is the independent error ("nugget").
With `compare = TRUE` the nested sequence RC, LV with one \\\phi\\, LV
with \\\phi_R, \phi_C\\, LV with \\\phi\_{RC}\\ as well, and AR + RC are
all fitted and compared by REML deviance (likelihood-ratio tests for the
nested LV terms, and AIC), which is the sequential fitting process the
paper recommends for deciding whether spatial terms are needed.

## Usage

``` r
analyze_resolvable_row_column(
  data,
  response,
  model = c("rc", "lv", "ar"),
  lv_terms = c("phi_r", "phi_c"),
  equal_phi = FALSE,
  long_col = NULL,
  compare = TRUE,
  genotype_col = "all_entries",
  replicate_col = "replicate",
  row_col = "row_in_rep",
  col_col = "col_in_rep"
)
```

## Arguments

- data:

  A design data frame (for example `design$design` from
  [`resolvable_row_column_design()`](https://chris263.github.io/breedtools/reference/resolvable_row_column_design.md))
  or a design object, with a response.

- response:

  Name of the numeric response column.

- model:

  Variance model reported in `adjusted_means`: `"rc"`, `"lv"`, or
  `"ar"`.

- lv_terms:

  For `model = "lv"`, the spatial terms to fit: any of `"phi_r"`,
  `"phi_c"`, `"phi_rc"`.

- equal_phi:

  Logical. For `model = "lv"`, use a single \\\phi\\ for rows and
  columns.

- long_col:

  Optional name of a column with long-column (or long-row) labels to fit
  as fixed effects, as in latinized designs.

- compare:

  Logical. If `TRUE`, also fit and compare the sequence of models.

- genotype_col, replicate_col, row_col, col_col:

  Column names.

## Value

An object of class `resolvable_row_column_analysis`, a list with

- plot_results:

  One row per plot: the data, the fitted value (fixed effects plus the
  predicted row, column and spatial effects) and the residual
  \\V^{-1}(y - X\hat\beta)\\.

- adjusted_means:

  Treatment BLUEs (average over replicates) with standard errors, ranked
  best first.

- variance_parameters:

  The REML estimates: \\\sigma^2\\, the \\\eta\\, \\\phi\\, \\\psi\\ and
  \\\rho\\ parameters, and the implied row and column variance
  components \\\sigma^2\eta_R/s\\ and \\\sigma^2\eta_C/k\\.

- model_comparison:

  With `compare = TRUE`, the deviance table.

- anova:

  Intra-block ANOVA from the fixed model
  `replicate + rows(replicate) + columns(replicate) + treatment`.

- diagnostics:

  REML deviance, number of variance parameters, mean SED and approximate
  5% LSD.

## Details

For latinized designs the long columns (or long rows) can be fitted as
fixed effects with `long_col`.

## References

Williams, E. R., John, J. A. and Whitaker, D. (2006). Construction of
resolvable spatial row-column designs. Biometrics 62, 103-108.
[doi:10.1111/j.1541-0420.2005.00393.x](https://doi.org/10.1111/j.1541-0420.2005.00393.x)

## Examples

``` r
design <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5,
                                       n_starts = 1, n_iter = 300, seed = 1)
trial <- design$design
set.seed(2)
trial$y <- 10 + 0.3 * trial$row + stats::rnorm(nrow(trial))
fit <- analyze_resolvable_row_column(trial, "y", model = "lv")
fit$model_comparison
#> # A tibble: 5 × 7
#>   model          n_variance_parameters deviance   AIC change_vs_rc compared_with
#>   <chr>                          <int>    <dbl> <dbl>        <dbl> <chr>        
#> 1 RC                                 3     127.  133.         0    NA           
#> 2 RC + LV (sing…                     4     123.  131.         3.89 RC           
#> 3 RC + LV (phi_…                     5     121.  131.         5.46 RC           
#> 4 RC + LV (phi_…                     6     121.  133.         5.46 RC + LV (phi…
#> 5 RC + AR1 x AR1                     6     127.  139.         0    NA           
#> # ℹ 1 more variable: lrt_p_value <dbl>
head(fit$adjusted_means)
#> # A tibble: 6 × 5
#>    rank genotype n_observations adjusted_mean std_error
#>   <int> <chr>             <int>         <dbl>     <dbl>
#> 1     1 T16                   3          14.3      5.53
#> 2     2 T13                   3          12.9      5.53
#> 3     3 T8                    3          12.4      5.53
#> 4     4 T7                    3          12.3      5.54
#> 5     5 T19                   3          12.3      5.54
#> 6     6 T17                   3          12.2      5.53
```
