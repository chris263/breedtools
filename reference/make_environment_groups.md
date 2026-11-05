# Test all compatible environment groups using Hartley Fmax

Tests all combinations of environments and identifies groups with
homogeneous residual variances.

## Usage

``` r
make_environment_groups(
  env_anova,
  alpha = 0.05,
  max_group_size = Inf,
  n_sim = 1e+05,
  seed = 123
)
```

## Arguments

- env_anova:

  A data frame from
  [`analyze_environment_anova()`](https://chris263.github.io/breedtools/reference/analyze_environment_anova.md).

- alpha:

  Significance level.

- max_group_size:

  Maximum environment group size to test.

- n_sim:

  Number of simulations for Hartley critical value.

- seed:

  Random seed.

## Value

A tibble with all tested groups.

## Examples

``` r
env_anova <- data.frame(
  environment = c("Loc1", "Loc2", "Loc3"),
  ms_error = c(10, 12, 30),
  df_error = c(20, 20, 20)
)

make_environment_groups(env_anova, n_sim = 2000)
#> # A tibble: 4 × 9
#>   group_id     k environment_list   min_ms_error max_ms_error fmax_observed
#>   <chr>    <int> <chr>                     <dbl>        <dbl>         <dbl>
#> 1 G4           3 Loc1 | Loc2 | Loc3           10           30           3  
#> 2 G1           2 Loc1 | Loc2                  10           12           1.2
#> 3 G3           2 Loc2 | Loc3                  12           30           2.5
#> 4 G2           2 Loc1 | Loc3                  10           30           3  
#> # ℹ 3 more variables: fmax_critical <dbl>, alpha <dbl>, homogeneous <lgl>
```
