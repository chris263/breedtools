# Run Breedbase multi-location phenotype variance check

Reads a Breedbase phenotype file, runs individual ANOVA by environment,
performs Hartley Fmax tests, and identifies compatible environment
groups for joint analysis.

## Usage

``` r
run_breedbase_variance_check(
  file,
  trait_col,
  alpha = 0.05,
  location_col = NULL,
  genotype_col = NULL,
  rep_col = NULL,
  block_col = NULL,
  max_group_size = Inf,
  n_sim = 1e+05,
  seed = 123,
  export = FALSE,
  output_dir = NULL
)
```

## Arguments

- file:

  Path to a `.xlsx`, `.xls`, or `.csv` phenotype file, or an already
  loaded data frame.

- trait_col:

  Trait column to analyze.

- alpha:

  Significance level for Hartley Fmax test.

- location_col:

  Optional location/environment column.

- genotype_col:

  Optional genotype/accession column.

- rep_col:

  Optional replicate column.

- block_col:

  Optional block column.

- max_group_size:

  Maximum environment group size to test.

- n_sim:

  Number of simulations for Hartley critical value.

- seed:

  Random seed for the Hartley critical value simulation. The user's
  global random number generator state is not modified.

- export:

  Logical. If `TRUE`, writes CSV files to `output_dir`.

- output_dir:

  Directory where CSV files should be written. Defaults to `NULL`, in
  which case nothing is written; it must be supplied when
  `export = TRUE`.

## Value

A list containing cleaned data, ANOVA results, Hartley test, all group
tests, compatible groups, maximal compatible groups, and detected
columns.

## Examples

``` r
pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")

result <- run_breedbase_variance_check(
  file = pheno_file,
  trait_col = "fresh root yield|CO_334:0000013",
  n_sim = 2000
)

result$env_anova
#> # A tibble: 4 × 12
#>   environment n_genotypes r_eff df_genotype ms_genotype df_error ms_error  mean
#>   <chr>             <int> <dbl>       <int>       <dbl>    <int>    <dbl> <dbl>
#> 1 Geneva               12     2          11        20.6       11     1.14  31.5
#> 2 Ibadan               12     2          11        39.7       11     6.45  22.1
#> 3 Ithaca               12     2          11        14.3       11     1.25  27.7
#> 4 Ubiaja               12     2          11        54.3       11     4.78  24.3
#> # ℹ 4 more variables: cv_percent <dbl>, sigma2_g <dbl>, f_value <dbl>,
#> #   p_value <dbl>
result$maximal_compatible_groups
#> # A tibble: 2 × 9
#>   group_id     k environment_list        min_ms_error max_ms_error fmax_observed
#>   <chr>    <int> <chr>                          <dbl>        <dbl>         <dbl>
#> 1 G9           3 Geneva | Ithaca | Ubia…         1.14         4.78          4.20
#> 2 G5           2 Ibadan | Ubiaja                 4.78         6.45          1.35
#> # ℹ 3 more variables: fmax_critical <dbl>, alpha <dbl>, homogeneous <lgl>

# Export CSV files to a temporary directory
out_dir <- file.path(tempdir(), "variance_check")
result <- run_breedbase_variance_check(
  file = pheno_file,
  trait_col = "fresh root yield|CO_334:0000013",
  n_sim = 2000,
  export = TRUE,
  output_dir = out_dir
)
list.files(out_dir)
#> [1] "hartley_all_group_tests_fresh root yield|CO_334:0000013.csv"          
#> [2] "hartley_maximal_compatible_groups_fresh root yield|CO_334:0000013.csv"
#> [3] "individual_anova_fresh root yield|CO_334:0000013.csv"                 
unlink(out_dir, recursive = TRUE)
```
