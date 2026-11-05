# Plot a resolvable row-column design

Draws the field with replicate boundaries.

## Usage

``` r
plot_resolvable_row_column_design(
  design,
  label = c("all_entries", "replicate", "none"),
  fill = c("replicate", "is_filler"),
  main = NULL,
  ...
)
```

## Arguments

- design:

  An object from
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)
  or its `design` data frame.

- label:

  Cell label: `"all_entries"`, `"replicate"`, or `"none"`.

- fill:

  Cell colouring: `"replicate"` or `"is_filler"`.

- main:

  Optional title.

- ...:

  Passed to
  [`plot_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_row_column_design.md).

## Value

Invisibly, the plotted data.

## Examples

``` r
design <- resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5,
                                       n_starts = 1, n_iter = 300, seed = 1)
plot_resolvable_row_column_design(design)

```
