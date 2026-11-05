# Calculate augmented row-column design capacity

Reports the plot counts, error degrees of freedom and feasibility of an
augmented row-column layout with superimposed blocks (Piepho and
Williams, 2016). Rows are grouped into row groups and columns into
column groups; their intersections are the blocks. Groups can have equal
sizes (`rows_per_block`, `cols_per_block`) or unequal sizes
(`rows_per_rowgroup`, `cols_per_colgroup`), and blocks can hold any
number of check plots (`checks_per_block`).

## Usage

``` r
augmented_row_column_capacity(
  rows_in_field = NULL,
  cols_in_field = NULL,
  rows_per_block = NULL,
  cols_per_block = NULL,
  n_controls,
  rows_per_rowgroup = NULL,
  cols_per_colgroup = NULL,
  checks_per_block = NULL
)
```

## Arguments

- rows_in_field:

  Number of rows in the whole field.

- cols_in_field:

  Number of columns in the whole field.

- rows_per_block:

  Number of rows in each row group (equal groups).

- cols_per_block:

  Number of columns in each column group (equal groups).

- n_controls:

  Number of replicated checks/controls (`vc`).

- rows_per_rowgroup:

  Optional vector with the number of rows of each row group (`kb(i)`),
  for unequal groups. Replaces `rows_per_block`.

- cols_per_colgroup:

  Optional vector with the number of columns of each column group
  (`sb(j)`). Replaces `cols_per_block`.

- checks_per_block:

  Optional number of check plots per block (`pb(i, j)`): a single
  number, a matrix with one row per row group and one column per column
  group, or a vector by row group. Defaults to `n_controls` (every block
  a complete block of checks). Smaller values give incomplete blocks;
  larger values give extra check replication.

## Value

A one-row tibble with the field size, number of groups and blocks, check
and entry plots, `error_df`, the error degrees of freedom of the
row-column ANOVA,
`check_plots - n_controls - (rows - 1) - (cols - 1) - (row groups - 1)(column groups - 1)`
(a design is only analyzable when `error_df > 0`; Clarke and Stefanova
(2011) recommend about 50), `check_fraction`, and the smallest slack of
the row-group and column-group conditions of Piepho and Williams (2016),
which must be non-negative for the row-group and column-group
sub-designs to be estimable on their own.

## References

Piepho, H.-P. and Williams, E. R. (2016). Augmented row-column designs
for a small number of checks. Agronomy Journal 108, 2256-2262.
[doi:10.2134/agronj2016.06.0325](https://doi.org/10.2134/agronj2016.06.0325)

## Examples

``` r
augmented_row_column_capacity(
  rows_in_field = 12,
  cols_in_field = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_controls = 4
)
#> # A tibble: 1 × 17
#>   rows_in_field cols_in_field rows_per_block cols_per_block n_block_rows
#>           <int>         <int>          <int>          <int>        <int>
#> 1            12            12              4              4            3
#> # ℹ 12 more variables: n_block_cols <int>, n_blocks <int>, n_controls <int>,
#> #   total_plots <int>, check_plots <int>, entry_plots <int>, error_df <dbl>,
#> #   check_fraction <dbl>, min_rowgroup_slack <dbl>, min_colgroup_slack <dbl>,
#> #   every_row_has_check <lgl>, every_col_has_check <lgl>

# Example 3 of Piepho and Williams (2016): unequal groups and blocks.
augmented_row_column_capacity(
  n_controls = 4,
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
)
#> # A tibble: 1 × 17
#>   rows_in_field cols_in_field rows_per_block cols_per_block n_block_rows
#>           <int>         <int>          <int>          <int>        <int>
#> 1            22            10             NA             NA            4
#> # ℹ 12 more variables: n_block_cols <int>, n_blocks <int>, n_controls <int>,
#> #   total_plots <int>, check_plots <int>, entry_plots <int>, error_df <dbl>,
#> #   check_fraction <dbl>, min_rowgroup_slack <dbl>, min_colgroup_slack <dbl>,
#> #   every_row_has_check <lgl>, every_col_has_check <lgl>
```
