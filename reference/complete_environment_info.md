# Complete environment metadata

Complete environment metadata

## Usage

``` r
complete_environment_info(env_info)
```

## Arguments

- env_info:

  Data frame with environment, location, trial, and year.

## Value

Data frame with an added `environment_label` column.

## Examples

``` r
env_info <- data.frame(
  environment = c("Loc1_T1", "Loc2_T2"),
  location = c("Loc1", "Loc2"),
  trial = c("T1", "T2"),
  year = c("2025", NA)
)

complete_environment_info(env_info)
#>   environment location trial year environment_label
#> 1     Loc1_T1     Loc1    T1 2025  Loc1 / T1 / 2025
#> 2     Loc2_T2     Loc2    T2              Loc2 / T2
```
