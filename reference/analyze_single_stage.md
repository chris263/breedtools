# Single-stage analysis of several trials

Fits all trials (locations, environments) together in one mixed model
and gives the genotype adjusted means across trials, to compare with (or
use instead of) the two-stage analysis of
[`analyze_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_two_stage.md):
\$\$y = \mu + \text{trial} + \text{genotype} + \text{rep(trial)} +
\text{block(rep, trial)} + \text{genotype} \times \text{trial} + e,\$\$
with trials and genotypes fixed, and replicates, blocks and the
genotype-by-trial interaction random (fitted with 'lme4'). Genotype
means average the trial effects with equal weights. Without 'lme4' the
model is fitted with all effects fixed and no genotype-by-trial term.

## Usage

``` r
analyze_single_stage(
  data,
  response,
  genotype_col = "accession",
  location_col = "trial_name",
  replicate_col = NULL,
  block_col = NULL,
  genotype_by_location = c("random", "none"),
  spatial = FALSE,
  row_col = "row",
  col_col = "col",
  spatial_nseg = NULL
)
```

## Arguments

- data:

  A data frame with one row per plot, e.g. from
  [`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md).

- response:

  Name of the numeric response column.

- genotype_col, location_col:

  Genotype and trial (location) columns.

- replicate_col, block_col:

  Optional replicate and block columns; by default `"rep"`/`"replicate"`
  and `"block"` when the data have them.

- genotype_by_location:

  `"random"` (default) or `"none"`.

- spatial:

  Logical. Fit the single-stage model with 'SpATS' (needs the field
  `row` and `col` of every plot).

- row_col, col_col:

  Field row and column columns (for `spatial = TRUE`).

- spatial_nseg:

  Optional number of P-spline segments, `c(col, row)`.

## Value

An object of class `single_stage_analysis`: `adjusted_means` (ranked,
with standard errors and the number of trials), `adjusted_means_vcov`,
`variance_components`, `plot_results` (fitted values and residuals),
`diagnostics` and the fitted `model`.

## Details

With `spatial = TRUE` the model is fitted with 'SpATS': genotypes and
trials fixed, replicates within trials fixed, blocks and the
genotype-by-trial interaction random, and a smooth P-spline surface over
the field rows and columns. SpATS fits one surface per model, so the
trials are placed side by side in one coordinate system, separated by
empty columns: each trial keeps its own spatial trend and the trends do
not run from one trial into another.

One error variance is shared by all trials here; when trials differ much
in precision, the weighted two-stage analysis is the better choice
(Möhring and Piepho, 2009).

## References

Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in
two-stage analysis of plant breeding trials. Crop Science 49, 1977-1988.
[doi:10.2135/cropsci2009.02.0083](https://doi.org/10.2135/cropsci2009.02.0083)

## Examples

``` r
sites <- c("North", "South")
met <- do.call(rbind, lapply(1:2, function(i) {
  d <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                   n_candidates = 5, seed = i)$design
  d$location <- sites[i]
  d
}))
set.seed(1)
g <- stats::setNames(stats::rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
met$yield <- 50 + g[met$all_entries] + stats::rnorm(nrow(met))
one <- analyze_single_stage(met, "yield", genotype_col = "all_entries", location_col = "location")
#> boundary (singular) fit: see help('isSingular')
one
#> Single-stage analysis (mixed (lme4)): 288 plots, 2 trials, 112 genotypes
#> 
#> Variance components:
#>                 component variance
#>          genotype x trial  0.13311
#>  blocks within replicates  0.07254
#>  replicates within trials  0.00000
#>                  residual  0.86458
#> 
#> Top adjusted means (112 genotypes):
#>  rank genotype n_locations adjusted_mean std_error
#>     1      G61           2         58.55    0.7256
#>     2      G70           2         56.69    0.7235
#>     3      G95           2         55.69    0.7251
#>     4      G31           2         55.14    0.7222
#>     5      G56           2         54.95    0.7247
#>     6       C2           2         54.94    0.3442
#>     7       G4           2         54.87    0.7252
#>     8     G106           2         54.77    0.7256
#>     9      G11           2         54.72    0.7247
#>    10      G68           2         54.50    0.7261
#> Mean SED = 1.007
```
