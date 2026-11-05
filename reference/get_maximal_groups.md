# Select maximal compatible environment groups

Removes smaller compatible groups that are fully contained inside larger
compatible groups.

## Usage

``` r
get_maximal_groups(compatible_groups)
```

## Arguments

- compatible_groups:

  Output from
  [`get_compatible_groups()`](https://chris263.github.io/breedtools/reference/get_compatible_groups.md).

## Value

A tibble with maximal compatible groups only.

## Examples

``` r
env_anova <- data.frame(
  environment = c("Loc1", "Loc2", "Loc3"),
  ms_error = c(10, 12, 30),
  df_error = c(20, 20, 20)
)

group_tests <- make_environment_groups(env_anova, n_sim = 2000)
get_maximal_groups(get_compatible_groups(group_tests))
#> # A tibble: 1 × 9
#>   group_id     k environment_list   min_ms_error max_ms_error fmax_observed
#>   <chr>    <int> <chr>                     <dbl>        <dbl>         <dbl>
#> 1 G4           3 Loc1 | Loc2 | Loc3           10           30             3
#> # ℹ 3 more variables: fmax_critical <dbl>, alpha <dbl>, homogeneous <lgl>
```
