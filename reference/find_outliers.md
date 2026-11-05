# Find outliers in trait data (IQR rule)

Flags values below \\Q_1 - k \cdot IQR\\ or above \\Q_3 + k \cdot IQR\\,
where \\Q_1\\, \\Q_3\\ are the quartiles and \\k\\ is `multiplier`
(Tukey's fences; 1.5 by default, 3 for only the extreme values). The
rule is applied within each group of `by`, usually each location,
because locations have different means.

## Usage

``` r
find_outliers(
  data,
  traits,
  multiplier = 1.5,
  by = if ("location" %in% names(data)) "location" else NULL,
  id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries",
    "accession", "row", "col")
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

## Value

A tibble with one row per outlier: the identifying columns, the group,
`trait`, `value`, the fences `lower` and `upper`, and `side` (`"low"` or
`"high"`).

## Examples

``` r
d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
                yield = c(rnorm(18, 50, 2), 80, 10))
find_outliers(d, "yield")
#> # A tibble: 2 × 8
#>   plots all_entries trait value lower upper side   .row
#>   <int> <chr>       <chr> <dbl> <dbl> <dbl> <chr> <int>
#> 1    19 G19         yield    80  40.5  58.4 high     19
#> 2    20 G20         yield    10  40.5  58.4 low      20
```
