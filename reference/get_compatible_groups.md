# Get compatible environment groups

Get compatible environment groups

## Usage

``` r
get_compatible_groups(group_tests)
```

## Arguments

- group_tests:

  Output from
  [`make_environment_groups()`](https://solgenomics.github.io/breedtools/reference/make_environment_groups.md).

## Value

A tibble containing only homogeneous groups.

## Examples

``` r
env_anova <- data.frame(
  environment = c("Loc1", "Loc2", "Loc3"),
  ms_error = c(10, 12, 30),
  df_error = c(20, 20, 20)
)

group_tests <- make_environment_groups(env_anova, n_sim = 2000)
get_compatible_groups(group_tests)
#> # A tibble: 3 × 9
#>   group_id     k environment_list   min_ms_error max_ms_error fmax_observed
#>   <chr>    <int> <chr>                     <dbl>        <dbl>         <dbl>
#> 1 G4           3 Loc1 | Loc2 | Loc3           10           30           3  
#> 2 G1           2 Loc1 | Loc2                  10           12           1.2
#> 3 G3           2 Loc2 | Loc3                  12           30           2.5
#> # ℹ 3 more variables: fmax_critical <dbl>, alpha <dbl>, homogeneous <lgl>
```
