# Simulate Hartley Fmax critical value

Estimates the Hartley Fmax critical value using simulation. This version
can handle different residual degrees of freedom across environments.

## Usage

``` r
hartley_critical_sim(df_error_vec, alpha = 0.05, n_sim = 1e+05, seed = 123)
```

## Arguments

- df_error_vec:

  Numeric vector of residual degrees of freedom.

- alpha:

  Significance level. Default is `0.05`.

- n_sim:

  Number of simulations. Default is `100000`.

- seed:

  Random seed used for the simulation. The user's global random number
  generator state is restored on exit. Use `NULL` to draw from the
  current random number stream without setting a seed.

## Value

A numeric critical value.

## Examples

``` r
hartley_critical_sim(c(20, 20, 20), alpha = 0.05, n_sim = 2000)
#> [1] 3.053145
```
