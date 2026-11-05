# List the analyzable layouts of an augmented row-column field

Tries every combination of block size and number of checks for a field
and returns the layouts that can be analysed, so a layout can be chosen
without writing a loop over
[`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md).
Each block gets one plot of every check. Suggested blocks have at least
as many rows and at least as many columns as there are checks, leaving
room for each check to occupy a distinct row and column.

## Usage

``` r
augmented_row_column_layouts(
  rows_in_field,
  cols_in_field,
  rows_per_block = NULL,
  cols_per_block = NULL,
  n_controls = 2:6,
  min_error_df = 1,
  n_entries = NULL
)
```

## Arguments

- rows_in_field, cols_in_field:

  Field size.

- rows_per_block, cols_per_block:

  Block sizes to try. By default every divisor of the field size from 2
  up to the field size.

- n_controls:

  Numbers of checks to try.

- min_error_df:

  Smallest error df to keep.

- n_entries:

  Optional number of entries to plant. Only layouts with at least that
  many entry plots are kept, and `spare_plots` shows how many entry
  plots would be left over.

## Value

A tibble with one row per analyzable layout, sorted by entry plots (most
first) and then error df: block sizes, number of checks, blocks, entry
and check plots, `check_fraction`, `error_df`, `method` (the
construction that
[`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md)
uses for it: `"latin"` or `"optimal"`) and, when `n_entries` is given,
`spare_plots`.

## Details

There is a trade-off: more check plots mean fewer entries but more error
degrees of freedom, i.e. more reliable F-tests and spatial adjustment.
Aim for at least about 10 error df when the field allows it (Clarke and
Stefanova, 2011, recommend about 50 for large trials).

## See also

[`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md),
[`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md)

## Examples

``` r
# Suggested analyzable layouts of a 12 x 12 field with 2 to 5 checks:
augmented_row_column_layouts(12, 12, n_controls = 2:5)
#> # A tibble: 9 × 9
#>   rows_per_block cols_per_block n_controls entry_plots error_df check_plots
#>            <int>          <int>      <int>       <int>    <dbl>       <int>
#> 1              4              4          4         108        6          36
#> 2              4              3          3         108        5          36
#> 3              3              4          3         108        5          36
#> 4              4              2          2         108        2          36
#> 5              2              4          2         108        2          36
#> 6              3              3          3          96       14          48
#> 7              3              2          2          96        9          48
#> 8              2              3          2          96        9          48
#> 9              2              2          2          72       23          72
#> # ℹ 3 more variables: n_blocks <int>, check_fraction <dbl>, method <chr>

# Layouts that fit 100 entries with at least 10 error df:
augmented_row_column_layouts(12, 12, n_controls = 2:6, min_error_df = 10, n_entries = 100)
#> # A tibble: 0 × 10
#> # ℹ 10 variables: rows_per_block <int>, cols_per_block <int>, n_controls <int>,
#> #   entry_plots <int>, error_df <dbl>, check_plots <int>, n_blocks <int>,
#> #   check_fraction <dbl>, method <chr>, spare_plots <int>
```
