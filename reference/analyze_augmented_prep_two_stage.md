# Two-stage analysis of a multi-location augmented p-rep trial

Analyses each location separately (stage 1) and combines the
per-location adjusted means across locations (stage 2), weighting them
by their precision. This is the approach suggested by Williams, Piepho
and Whitaker (2011, Section 3) for adding spatial error components and
for heterogeneous error variances across locations: every location gets
its own error variance (and, optionally, its own spatial trend), and the
stage-2 weights carry that information into the combined means (Möhring
and Piepho, 2009).

## Usage

``` r
analyze_augmented_prep_two_stage(
  data,
  response,
  stage1 = c("blocks", "spatial"),
  weighting = c("full", "diagonal", "none"),
  genotype_by_location = c("random", "none"),
  genotype_col = "all_entries",
  location_col = "location",
  block_col = "block",
  replicate_col = "replicate",
  row_col = "row",
  col_col = "col",
  spatial_nseg = NULL,
  spatial_random_blocks = TRUE
)
```

## Arguments

- data:

  A design data frame (usually `design$design` from
  [`augmented_prep_design()`](https://chris263.github.io/breedtools/reference/augmented_prep_design.md)),
  or a design object, with a response column.

- response:

  Name of the numeric response column.

- stage1:

  Stage-1 model: `"blocks"` or `"spatial"`.

- weighting:

  Stage-2 weighting: `"full"`, `"diagonal"`, or `"none"`.

- genotype_by_location:

  `"random"` (default) to estimate a genotype-by-location variance in
  stage 2, or `"none"` to omit it.

- genotype_col, location_col, block_col, replicate_col:

  Column names, as in
  [`analyze_augmented_prep_design()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_design.md).

- row_col, col_col:

  Field row and column coordinate columns, used by `stage1 = "spatial"`.

- spatial_nseg:

  Optional number of P-spline segments for columns and rows,
  `c(col, row)`. Defaults to half the number of distinct columns and
  rows of each location (at least 3, at most 40), which is much faster
  than one segment per row or column and usually as accurate.

- spatial_random_blocks:

  Logical. If `TRUE` (default), incomplete blocks are included as a
  random effect in the spatial model.

## Value

An object of class `augmented_prep_two_stage`, a list with

- adjusted_means:

  Combined genotype means ranked best first, with standard errors and
  the number of locations.

- variance_components:

  Stage-2 REML genotype-by-location variance.

- locations:

  One row per location: stage-1 model, plots, genotypes, residual
  variance and df, mean SED and, for spatial fits, the effective
  dimension of the spatial surface.

- heterogeneity:

  Hartley's Fmax ratio of the stage-1 residual variances with a
  simulated 5% critical value.

- stage1_means:

  The stage-1 genotype means and standard errors per location.

- stage1_fits:

  The stage-1 fits (analysis objects or `SpATS` objects), named by
  location.

- diagnostics:

  Settings, number of stage-2 means, REML log-likelihood, mean SED and
  approximate 5% LSD.

## Stage 1

- `stage1 = "blocks"`:

  The incomplete-block model of
  [`analyze_augmented_prep_design()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_design.md)
  at each location (replicates and blocks fixed, genotypes fixed).

- `stage1 = "spatial"`:

  A spatial model fitted with the 'SpATS' package (Rodríguez-Álvarez et
  al., 2018): genotypes fixed, a smooth two-dimensional P-spline surface
  over the field rows and columns (`PSANOVA`), replicates fixed and
  incomplete blocks random. It needs the `row` and `col` field
  coordinates (use `field_cols` in
  [`augmented_prep_design()`](https://chris263.github.io/breedtools/reference/augmented_prep_design.md)).

The stage-1 genotype means are BLUEs with their full covariance matrix.
Each location needs some error degrees of freedom, i.e. enough p-rep
entries for its blocks. In spatial fits part of the plot-to-plot noise
is taken up by the smooth surface, so the reported residual variance is
somewhat smaller than in the block model; the stage-1 covariance matrix
accounts for this, and the Hartley test on spatial residual variances is
approximate.

## Stage 2

The stage-1 means \\\bar y\_{gl}\\ are modelled as \\\bar y\_{gl} =
\mu + \lambda_l + \gamma_g + (\gamma\lambda)\_{gl} + e\_{gl}\\, with
fixed locations and genotypes, random genotype-by-location effects with
variance \\\sigma^2\_{GL}\\ (estimated by REML) and stage-1 errors \\e\\
whose covariance is taken from stage 1:

- `weighting = "full"`:

  the full stage-1 covariance matrix of each location (most accurate);

- `weighting = "diagonal"`:

  only the stage-1 variances (weights \\1 / se^2\\);

- `weighting = "none"`:

  an unweighted analysis of the stage-1 means with a common error, which
  ignores heterogeneity (for comparison).

Genotype adjusted means average the location effects with equal weights.

## References

Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented p-rep
designs. Biometrical Journal 53, 19-27.
[doi:10.1002/bimj.201000102](https://doi.org/10.1002/bimj.201000102)

Möhring, J. and Piepho, H.-P. (2009). Comparison of weighting in
two-stage analysis of plant breeding trials. Crop Science 49, 1977-1988.
[doi:10.2135/cropsci2009.02.0083](https://doi.org/10.2135/cropsci2009.02.0083)

Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A. and Eilers, P.
H. C. (2018). Correcting for spatial heterogeneity in plant breeding
experiments with P-splines. Spatial Statistics 23, 52-71.
[doi:10.1016/j.spasta.2017.10.003](https://doi.org/10.1016/j.spasta.2017.10.003)

## See also

[`analyze_augmented_prep_design()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_design.md)
for the one-stage analysis.

## Examples

``` r
design <- augmented_prep_design(
  treatments = sprintf("G%02d", 1:60),
  locations = c("A", "B", "C"),
  method = "alpha",
  block_size = 6,
  field_cols = 9,
  seed = 1
)
trial <- design$design
set.seed(2)
# Location C is twice as noisy as A and B.
trial$yield <- 50 + stats::rnorm(nrow(trial), sd = ifelse(trial$location == "C", 4, 2))

fit <- analyze_augmented_prep_two_stage(trial, response = "yield")
fit
#> Two-stage augmented p-rep analysis
#> Stage 1: blocks model per location; stage 2: full weighting, genotype x location random
#> 
#> Locations:
#>  location n_plots n_genotypes residual_variance residual_df mean_sed
#>         A      80          60             7.259           7    4.908
#>         B      80          60             2.267           7    2.744
#>         C      80          60            16.970           7    7.338
#> 
#> Hartley Fmax = 7.49 (5% critical value 6.93): error variances differ between locations
#> 
#> Stage-2 variance components:
#>            component variance
#>  genotype x location     1.39
#> 
#> Top adjusted means (60 genotypes):
#>  rank genotype n_locations adjusted_mean std_error
#>     1      G33           3         53.16     1.612
#>     2      G54           3         53.16     1.582
#>     3      G46           3         52.59     1.580
#>     4      G39           3         52.10     1.551
#>     5      G26           3         52.00     1.529
#>     6      G40           3         51.86     1.460
#>     7      G10           3         51.59     1.615
#>     8      G25           3         51.57     1.473
#>     9      G06           3         51.52     1.581
#>    10      G56           3         51.36     1.466
#> Mean SED = 2.162; LSD (5%) = 4.281
fit$locations
#> # A tibble: 3 × 8
#>   location stage1_model n_plots n_genotypes residual_variance residual_df
#>   <chr>    <chr>          <int>       <int>             <dbl>       <int>
#> 1 A        blocks            80          60              7.26           7
#> 2 B        blocks            80          60              2.27           7
#> 3 C        blocks            80          60             17.0            7
#> # ℹ 2 more variables: mean_sed <dbl>, spatial_effective_dim <dbl>

# \donttest{
if (requireNamespace("SpATS", quietly = TRUE)) {
  spatial <- analyze_augmented_prep_two_stage(trial, "yield", stage1 = "spatial")
  spatial$locations
}
#> # A tibble: 3 × 8
#>   location stage1_model n_plots n_genotypes residual_variance residual_df
#>   <chr>    <chr>          <int>       <int>             <dbl>       <dbl>
#> 1 A        spatial           80          60              4.01        14.0
#> 2 B        spatial           80          60              2.79        12.0
#> 3 C        spatial           80          60             17.0         15.5
#> # ℹ 2 more variables: mean_sed <dbl>, spatial_effective_dim <dbl>
# }
```
