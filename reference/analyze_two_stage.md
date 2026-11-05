# Two-stage analysis of a multi-location trial, for any design

Analyses each location on its own (stage 1) and combines the genotype
means across locations (stage 2), weighting them by their precision
(Möhring and Piepho, 2009). Works for every design of the package grown
at several locations: augmented row-column, augmented p-rep and
resolvable row-column designs. For augmented p-rep designs,
[`analyze_augmented_prep_two_stage()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md)
gives the same analysis with the p-rep block model of the paper.

## Usage

``` r
analyze_two_stage(
  data,
  response,
  stage1 = c("spatial", "blocks", "design"),
  weighting = c("full", "diagonal", "none"),
  genotype_by_location = c("random", "none"),
  genotype_col = "all_entries",
  location_col = "location",
  row_col = "row",
  col_col = "col",
  replicate_col = NULL,
  block_col = NULL,
  spatial_nseg = NULL
)
```

## Arguments

- data:

  A data frame with all locations (one row per plot), with the design
  columns and a response column.

- response:

  Name of the numeric response column.

- stage1:

  `"spatial"`, `"blocks"` or `"design"`.

- weighting:

  `"full"`, `"diagonal"` or `"none"`.

- genotype_by_location:

  `"random"` or `"none"`.

- genotype_col, location_col:

  Genotype and location columns.

- row_col, col_col:

  Field row and column coordinates (spatial stage 1).

- replicate_col, block_col:

  Optional replicate (fixed) and block (random) columns of the spatial
  stage 1. By default `"replicate"` and `"block"` when the data have
  them.

- spatial_nseg:

  Optional number of P-spline segments, `c(col, row)`.

## Value

An object of class `two_stage_analysis`: a list with `adjusted_means`
(combined, ranked), `variance_components`, `locations` (one row per
location), `heterogeneity` (Hartley's Fmax), `stage1_means`,
`stage1_fits` (per location, with `plot_results`, for
[`first_stage_values()`](https://chris263.github.io/breedtools/reference/first_stage_values.md)
and
[`compare_first_stage()`](https://chris263.github.io/breedtools/reference/compare_first_stage.md))
and `diagnostics`.

## Stage 1

- `stage1 = "spatial"` (default):

  A 'SpATS' spatial model at every location (Rodríguez-Álvarez et al.,
  2018): genotypes fixed, a smooth two-dimensional P-spline surface over
  the field rows and columns, replicates fixed and blocks random when
  `replicate_col` and `block_col` are given.

- `stage1 = "blocks"`:

  A linear model at every location that needs no design object:
  replicates and blocks within replicates fixed (when `replicate_col`
  and `block_col` are given), genotypes fixed. For data read from
  Breedbase without field coordinates.

- `stage1 = "design"`:

  The design's own model at every location:
  [`analyze_augmented_row_column_design()`](https://chris263.github.io/breedtools/reference/analyze_augmented_row_column_design.md)
  for augmented row-column designs,
  [`analyze_resolvable_row_column()`](https://chris263.github.io/breedtools/reference/analyze_resolvable_row_column.md)
  (row-column model) for resolvable row-column designs, and the block
  model of
  [`analyze_augmented_prep_design()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_design.md)
  for p-rep designs. The design is recognised from its columns.

The stage-1 genotype means are BLUEs with their full covariance matrix.

## Stage 2

As in
[`analyze_augmented_prep_two_stage()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md):
fixed locations and genotypes, random genotype-by-location effects
(REML), and stage-1 errors with the stage-1 covariance
(`weighting = "full"`), only the stage-1 variances (`"diagonal"`) or a
common error (`"none"`, for comparison).

## References

Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in
two-stage analysis of plant breeding trials. Crop Science 49, 1977-1988.
[doi:10.2135/cropsci2009.02.0083](https://doi.org/10.2135/cropsci2009.02.0083)

Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A. and Eilers, P.
H. C. (2018). Correcting for spatial heterogeneity in plant breeding
experiments with P-splines. Spatial Statistics 23, 52-71.
[doi:10.1016/j.spasta.2017.10.003](https://doi.org/10.1016/j.spasta.2017.10.003)

## Examples

``` r
# The same augmented row-column layout, randomized at three locations.
sites <- c("North", "Central", "South")
trials <- lapply(seq_along(sites), function(i) {
  d <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                   12, 12, 4, 4, n_candidates = 5, seed = i)$design
  d$location <- sites[i]
  d
})
met <- do.call(rbind, trials)
set.seed(1)
g <- stats::setNames(stats::rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
met$yield <- 50 + g[met$all_entries] + 0.5 * met$row + stats::rnorm(nrow(met))
fit <- analyze_two_stage(met, "yield", stage1 = "design")
fit
#> Two-stage analysis
#> Stage 1: design model per location; stage 2: full weighting, genotype x location random
#> 
#> Locations:
#>  location         stage1_model n_plots n_genotypes residual_variance
#>     North augmented row-column     144         112            1.2926
#>   Central augmented row-column     144         112            0.2687
#>     South augmented row-column     144         112            2.1303
#>  residual_df mean_sed
#>            6    2.856
#>            6    1.129
#>            6    3.119
#> 
#> Hartley Fmax = 7.93 (5% critical value 8.40): no evidence of heterogeneity
#> 
#> Stage-2 variance components:
#>            component variance
#>  genotype x location   0.1736
#> 
#> Top adjusted means (112 genotypes):
#>  rank genotype n_locations adjusted_mean std_error
#>     1      G61           3         60.53    0.6259
#>     2      G70           3         59.92    0.6283
#>     3     G106           3         58.88    0.6310
#>     4      G95           3         58.72    0.6281
#>     5       G4           3         58.56    0.6277
#>     6      G56           3         58.30    0.6283
#>     7       C2           3         58.13    0.3125
#>     8      G31           3         58.11    0.6256
#>     9      G11           3         57.70    0.6273
#>    10     G108           3         57.52    0.6251
#> Mean SED = 0.876; LSD (5%) = 1.726
# \donttest{
if (requireNamespace("SpATS", quietly = TRUE)) {
  spatial <- analyze_two_stage(met, "yield")      # SpATS at every location
  spatial$locations
}
#> # A tibble: 3 × 8
#>   location stage1_model    n_plots n_genotypes residual_variance residual_df
#>   <chr>    <chr>             <int>       <int>             <dbl>       <dbl>
#> 1 North    spatial (SpATS)     144         112             0.837        25.0
#> 2 Central  spatial (SpATS)     144         112             0.654        28.9
#> 3 South    spatial (SpATS)     144         112             1.20         29.0
#> # ℹ 2 more variables: mean_sed <dbl>, spatial_effective_dim <dbl>
# }
```
