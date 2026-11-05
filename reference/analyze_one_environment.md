# Analyze one environment using ANOVA

Runs an individual ANOVA for one environment. If more than one replicate
is present, the model is `value ~ replicate + genotype`; otherwise, the
model is `value ~ genotype`.

## Usage

``` r
analyze_one_environment(df)
```

## Arguments

- df:

  A data frame containing one environment with columns `environment`,
  `genotype`, `replicate`, and `value`.

## Value

A tibble with ANOVA summary statistics.

## Examples

``` r
df <- data.frame(
  environment = "Loc1",
  genotype = factor(rep(c("G1", "G2", "G3"), times = 2)),
  replicate = factor(rep(1:2, each = 3)),
  value = c(10, 12, 15, 11, 13, 14)
)

analyze_one_environment(df)
#> # A tibble: 1 × 12
#>   environment n_genotypes r_eff df_genotype ms_genotype df_error ms_error  mean
#>   <chr>             <int> <dbl>       <int>       <dbl>    <int>    <dbl> <dbl>
#> 1 Loc1                  3     2           2        8.00        2    0.667  12.5
#> # ℹ 4 more variables: cv_percent <dbl>, sigma2_g <dbl>, f_value <dbl>,
#> #   p_value <dbl>
```
