# Analyze an augmented row-column design

Fits the fixed-effects model:
`y = mean + rowgroup + colgroup + rowgroup:colgroup + rowgroup:row + colgroup:col + genotype + error`.
Sums of squares are sequential, in the order above. Degrees of freedom
are the ranks actually estimable from the observed (non-missing) plots.
For a complete design they equal the augmented row-column formulas,
reported in `anova_df`, which use field rows (`k`), field columns (`s`),
row groups (`gk`), column groups (`gs`), check cultivars (`vc`), and
unreplicated entries (`ve`): row groups `gk - 1`, rows within row groups
`k - gk`, column groups `gs - 1`, columns within column groups `s - gs`,
blocks `(gk - 1)(gs - 1)`, genotypes `vc + ve - 1`, and error as the
remainder of `k * s - 1`.

## Usage

``` r
analyze_augmented_row_column_design(
  data,
  response,
  genotype_col = "all_entries",
  rowgroup_col = "rowgroup",
  colgroup_col = "colgroup",
  row_col = "row",
  col_col = "col",
  type_col = "type",
  control_col = "is_control",
  random_blocks = FALSE,
  entries = c("fixed", "random"),
  spatial = FALSE,
  spatial_nseg = NULL
)
```

## Arguments

- data:

  A complete augmented row-column design data frame, usually
  `design$design` from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  with a response column added.

- response:

  Name of the numeric response column.

- genotype_col:

  Name of the genotype/treatment column.

- rowgroup_col:

  Name of the row-group column.

- colgroup_col:

  Name of the column-group column.

- row_col:

  Name of the field row column.

- col_col:

  Name of the field column column.

- type_col:

  Optional entry type column. Used only for genotype summaries.

- control_col:

  Optional control flag column. Identifies the checks for
  `entries = "random"`; otherwise used only for genotype summaries.

- random_blocks:

  Logical. If `TRUE`, rows within row groups, columns within column
  groups and blocks are random effects, which recovers inter-block
  information (Wolfinger et al., 1997; Möhring et al., 2015). Row and
  column groups stay fixed. Fitted by REML with 'lme4'.

- entries:

  `"fixed"` (default) or `"random"`. With `"random"`, checks are fixed
  and the unreplicated entries are random effects from a common breeding
  population (Piepho, Williams and Fleck, 2006); entry values are then
  BLUPs. Fitted by REML with 'lme4'.

- spatial:

  Logical. If `TRUE`, the plot errors get a spatial add-on: a
  two-dimensional P-spline surface over rows and columns fitted with the
  'SpATS' package, with row and column groups fixed and rows, columns
  and blocks random. Genotypes are fixed.

- spatial_nseg:

  Optional number of P-spline segments for columns and rows,
  `c(col, row)`. Defaults to half the number of columns and rows (at
  least 3, at most 40).

## Value

An object of class `augmented_row_column_analysis`: a list with

- input, formula, model:

  The column mapping, model formula, and an `lm` fit on the check plots.

- anova, anova_df, anova_formatted:

  Sequential ANOVA (estimable df), design df formulas, and a formatted
  ANOVA table.

- plot_results:

  The input data with fitted values and residuals.

- genotype_summary:

  One row per genotype with raw, fitted, and adjusted means and the
  standard error of the adjusted mean.

- adjusted_means:

  Genotypes ranked by adjusted mean (best first).

- diagnostics:

  Error df, residual standard deviation, R-squared, `mean_sed` (mean
  standard error of a difference between two adjusted means), `lsd_5pct`
  (the matching least significant difference), and flags for rank
  deficiency or a saturated fit.

- variance_components, mixed_model:

  For mixed or spatial fits: the REML variance components and the
  `lmerMod` or `SpATS` fit; otherwise `NULL`.

- message:

  A note when F-tests or adjusted means are not estimable.

## Details

Unreplicated entries each fit their own plot exactly, so the model is
fitted without a dense genotype design matrix: the residual comes from
the replicated (check) plots alone. `model` is the corresponding `lm`
fit on the check plots, with aliased columns removed. This keeps the
analysis fast for large trials.

## Mixed and spatial models

The ANOVA table always comes from the fixed-effects model of Piepho and
Williams (2016, Eq. 1). With `random_blocks`, `entries = "random"` or
`spatial`, the adjusted means, their standard errors, `mean_sed` and
`lsd_5pct` come from the corresponding REML fit instead. Standard errors
are exact prediction standard errors from the mixed-model equations for
the estimated variance components; `lsd_5pct` uses the fixed-model error
df as an approximation.

## Adjusted means

With the default fixed-effects model, genotype effects are fixed, so the
adjusted means are best linear unbiased estimates (BLUEs). Each is the
genotype effect plus the average row, column, and block effect over all
observed plots, i.e. the expected value of the genotype in an average
plot of the field. For an unreplicated entry this is its observed value
corrected for the row, column, and block effects of its plot, estimated
from the checks; its standard error combines the plot error and the
uncertainty of that correction. Use `adjusted_mean` (not `raw_mean`) to
rank and select entries.

## Examples

``` r
design <- augmented_row_column_design(
  treatments = paste0("G", seq_len(108)),
  controls = paste0("Check", seq_len(4)),
  rows_in_field = 12,
  cols_in_field = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_candidates = 3,
  seed = 123
)

design_data <- design$design
design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
analysis <- analyze_augmented_row_column_design(design_data, response = "yield")
analysis
#> Augmented row-column design analysis
#> 
#>                                 Source  Df      Sum Sq   Mean Sq  F value
#> 1                           Row groups   2    67.75024  33.87512 1.576604
#> 2       Rows, nested within row groups   9  1150.86304 127.87367 5.951450
#> 3                        Column groups   2    46.28071  23.14036 1.076990
#> 4 Columns, nested within column groups   9   260.04622  28.89402 1.344775
#> 5  Row groups x column groups (blocks)   4   340.33792  85.08448 3.959971
#> 6                            Genotypes 111 11383.48967 102.55396 4.773030
#> 7                                Error   6   128.91681  21.48614       NA
#> 8                      Corrected total 143 13377.68461        NA       NA
#>       Pr(>F) Significance
#> 1 0.28166563             
#> 2 0.02085954            *
#> 3 0.39842335             
#> 4 0.37093044             
#> 5 0.06585042             
#> 6 0.02684046            *
#> 7         NA             
#> 8         NA             
#> 
#> Top adjusted means (112 genotypes; see $adjusted_means):
#>  rank genotype n_observations adjusted_mean std_error
#>     1      G48              1         164.3     9.297
#>     2      G92              1         163.4     8.880
#>     3       G8              1         162.5     9.245
#>     4      G99              1         152.2     9.184
#>     5      G61              1         147.6     8.931
#>     6      G59              1         142.5     9.029
#>     7      G67              1         140.1     7.737
#>     8      G46              1         136.4     6.874
#>     9       G1              1         135.8     9.212
#>    10       G5              1         135.8     7.312
#> Mean SED = 11.27; LSD (5%) = 27.57
head(analysis$adjusted_means)
#> # A tibble: 6 × 8
#>    rank genotype type  is_control n_observations raw_mean adjusted_mean
#>   <int> <chr>    <chr>      <int>          <int>    <dbl>         <dbl>
#> 1     1 G48      entry          0              1     120.          164.
#> 2     2 G92      entry          0              1     125.          163.
#> 3     3 G8       entry          0              1     106.          162.
#> 4     4 G99      entry          0              1     110.          152.
#> 5     5 G61      entry          0              1     112.          148.
#> 6     6 G59      entry          0              1     105.          143.
#> # ℹ 1 more variable: std_error <dbl>
```
