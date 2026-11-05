# Calculate pairwise concurrence for a p-rep design

Counts how many times each pair of entries appears together in the same
incomplete block across all locations.

## Usage

``` r
augmented_prep_concurrence(design, treatments = NULL)
```

## Arguments

- design:

  A p-rep design data frame from
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md).

- treatments:

  Optional treatment names defining matrix order.

## Value

A symmetric integer matrix of pairwise block concurrences.

## Examples

``` r
design <- augmented_prep_design(
  treatments = paste0("G", seq_len(12)),
  locations = c("Loc1", "Loc2"),
  block_size = 4,
  n_candidates = 5,
  seed = 1
)
conc <- augmented_prep_concurrence(design$design)
table(conc[upper.tri(conc)])
#> 
#>  0  1  2 
#> 19 46  1 
```
