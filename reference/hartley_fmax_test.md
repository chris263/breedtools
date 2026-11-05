# Hartley Fmax test for residual variance homogeneity

Hartley Fmax test for residual variance homogeneity

## Usage

``` r
hartley_fmax_test(env_anova, alpha = 0.05, n_sim = 1e+05, seed = 123)
```

## Arguments

- env_anova:

  A data frame containing environment-level ANOVA results. It must
  include columns `environment`, `ms_error`, and `df_error`.

- alpha:

  Significance level. Default is `0.05`.

- n_sim:

  Number of simulations for the critical value.

- seed:

  Random seed passed to
  [`hartley_critical_sim()`](https://chris263.github.io/breedtools/reference/hartley_critical_sim.md).
  The user's global random number generator state is not modified.

## Value

A tibble with observed Fmax, critical Fmax, and homogeneity decision.

## Examples

``` r
env_anova <- data.frame(
  environment = c("Loc1", "Loc2", "Loc3"),
  ms_error = c(10, 12, 15),
  df_error = c(20, 20, 20)
)

hartley_fmax_test(env_anova, n_sim = 2000)
#> # A tibble: 1 × 8
#>   environments     k min_ms_error max_ms_error fmax_observed fmax_critical alpha
#>   <chr>        <int>        <dbl>        <dbl>         <dbl>         <dbl> <dbl>
#> 1 Loc1, Loc2,…     3           10           15           1.5          3.05  0.05
#> # ℹ 1 more variable: homogeneous <lgl>
```
