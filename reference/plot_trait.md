# Plot a trait: boxplot, distribution and field map

Draws the trait per location (or other group) so gradients, odd plots
and outliers are easy to see:

- `"boxplot"`:

  All locations side by side. Outliers (IQR rule, as in
  [`find_outliers()`](https://solgenomics.github.io/breedtools/reference/find_outliers.md),
  or the table given in `outliers`) are drawn in red and labelled with
  their plot number; the other plots in grey.

- `"histogram"`:

  The distribution at each location.

- `"field"`:

  A map of the field with the value of every plot, with a colour bar.
  The colour scale of each location spans the values inside the IQR
  fences, so a few extreme plots do not wash out the rest of the map;
  plots beyond the scale get the end colour and are circled in red.
  Missing plots are grey.

## Usage

``` r
plot_trait(
  data,
  trait,
  by = if ("location" %in% names(data)) "location" else NULL,
  type = c("both", "histogram", "field", "boxplot", "all"),
  row_col = "row",
  col_col = "col",
  outliers = NULL,
  breaks = "Sturges",
  multiplier = 1.5,
  scale = c("location", "common"),
  plot_col = intersect(c("plots", "plot_number"), names(data))[1]
)
```

## Arguments

- data:

  A data frame with one row per plot.

- trait:

  Name of the numeric trait column.

- by:

  Optional grouping column; by default `"location"` when the data have
  it.

- type:

  `"both"` (histogram and field map), `"histogram"`, `"field"`,
  `"boxplot"` or `"all"` (boxplot, then histogram and field map).

- row_col, col_col:

  Field row and column columns (for the map).

- outliers:

  Optional table from
  [`find_outliers()`](https://solgenomics.github.io/breedtools/reference/find_outliers.md)
  or `attr(remove_outliers(...), "outliers")`: these plots are marked.
  By default the outliers are found in `data` with the IQR rule.

- breaks:

  Histogram breaks, as in
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html).

- multiplier:

  IQR multiplier for the default outliers and the colour scale.

- scale:

  `"location"` (default: one colour scale per location, so each field
  shows its own gradient) or `"common"` (one scale for all locations, to
  compare their levels).

- plot_col:

  Column with the plot numbers used as labels.

## Value

Invisibly, a tibble with a summary per group: plots, values, missing,
outliers, mean, sd, min, max and CV (%).

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
d <- rc$design
d$yield <- 50 + 0.8 * d$row + rnorm(nrow(d), sd = 2)
d$yield[c(5, 40)] <- c(120, 10)
plot_trait(d, "yield", type = "boxplot")

plot_trait(d, "yield", type = "field")

```
