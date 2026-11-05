# Remove outliers from trait data (IQR rule)

Finds outliers with
[`find_outliers()`](https://chris263.github.io/breedtools/reference/find_outliers.md)
and sets those values to `NA`, so the plots stay in the design but the
values are left out of the analysis. The plots and values removed are
printed and returned with the data.

## Usage

``` r
remove_outliers(
  data,
  traits,
  multiplier = 1.5,
  by = if ("location" %in% names(data)) "location" else NULL,
  id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries",
    "accession", "row", "col"),
  quiet = FALSE
)
```

## Arguments

- data:

  A data frame with one row per plot.

- traits:

  Names of the numeric trait columns.

- multiplier:

  The IQR multiplier \\k\\.

- by:

  Optional grouping column(s); by default `"location"` when the data
  have it.

- id_cols:

  Columns that identify the plot in the report (those present in
  `data`).

- quiet:

  Logical. If `FALSE` (default), print what was removed.

## Value

`data` with the outlying values set to `NA`. The removed values are in
`attr(result, "outliers")` (as from
[`find_outliers()`](https://chris263.github.io/breedtools/reference/find_outliers.md)).

## Examples

``` r
d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
                yield = c(rnorm(18, 50, 2), 80, 10))
clean <- remove_outliers(d, "yield")
#> Removed 2 of 20 values as outliers (outside Q1 - 1.5 IQR and Q3 + 1.5 IQR):
#>  plots all_entries trait value lower upper side
#>     19         G19 yield    80 43.13 54.86 high
#>     20         G20 yield    10 43.13 54.86  low
attr(clean, "outliers")
#> # A tibble: 2 × 7
#>   plots all_entries trait value lower upper side 
#>   <int> <chr>       <chr> <dbl> <dbl> <dbl> <chr>
#> 1    19 G19         yield    80  43.1  54.9 high 
#> 2    20 G20         yield    10  43.1  54.9 low  
```
