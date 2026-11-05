# Second-stage genotype values from analyses

Collects the adjusted means of analyses (one per trait) into one table
with one row per accession, ready for
[`save_breedbase_analysis()`](https://chris263.github.io/breedtools/reference/save_breedbase_analysis.md).

## Usage

``` r
second_stage_values(fits, trait = NULL)
```

## Arguments

- fits:

  One analysis, or a named list of analyses, one per trait (the names
  are the Breedbase trait names). Any breedtools analysis with
  `adjusted_means`, e.g.
  [`analyze_augmented_prep_two_stage()`](https://chris263.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md).

- trait:

  Trait name when `fits` is a single analysis.

## Value

A tibble with `accession` and one column per trait. Its attribute
`"estimate"` is `"blup"` when the genotypes were random, else `"blue"`.

## Examples

``` r
prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
d <- prep$design
set.seed(1)
d$yield <- rnorm(nrow(d), 25, 3)
fit <- analyze_augmented_prep_two_stage(d, "yield")
head(second_stage_values(fit, trait = "fresh root yield"))
#> # A tibble: 6 × 2
#>   accession `fresh root yield`
#>   <chr>                  <dbl>
#> 1 E50                     30.0
#> 2 E45                     29.9
#> 3 E39                     28.6
#> 4 E11                     27.8
#> 5 E2                      27.6
#> 6 E28                     27.5
```
