# Efficiency of a resolvable row-column design

Computes the design criteria of Williams, John and Whitaker (2006) for a
resolvable row-column design with complete replicates: the average
efficiency factors of the row-column design, the row design and the
column design (rows and columns nested within replicates), the spatial
average efficiency factor of the two-dimensional linear variance model,
\$\$A = rI - X'(I_r \otimes W_R \otimes W_C)X,\$\$ with its upper bound
\\U = (v - t)/(v - 1)\\, \\t = tr(W_R) tr(W_C)\\, and treatment
neighbour and edge statistics. Optionally, it evaluates the average
pairwise variance criterion \\z(X)\\ under an assumed autoregressive
(AR) model, which is used to compare designs (Section 4 of the paper).

## Usage

``` r
resolvable_row_column_efficiency(
  design,
  ar = NULL,
  genotype_col = "all_entries",
  replicate_col = "replicate",
  row_col = "row_in_rep",
  col_col = "col_in_rep"
)
```

## Arguments

- design:

  A design object from
  [`resolvable_row_column_design()`](https://chris263.github.io/breedtools/reference/resolvable_row_column_design.md)
  or a data frame with replicate, row, column and treatment columns.
  Each replicate must be complete (`k * s` plots, each treatment once;
  filler plots count as separate treatments).

- ar:

  Optional named list of AR model parameters `eta_r`, `eta_c`, `psi`,
  `rho_r`, `rho_c` (Eq. 4 of the paper).

- genotype_col, replicate_col, row_col, col_col:

  Column names.

## Value

A list with `efficiency_rc`, `efficiency_rows`, `efficiency_cols`,
`efficiency_spatial`, `upper_bound_spatial`, `ar_criterion` (`NA` unless
`ar` is given), and `neighbours`, a tibble with each treatment's number
of distinct row and column neighbours and the number of replicates in
which it is on an edge or in a corner.

## Examples

``` r
design <- resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5,
                                       n_starts = 1, n_iter = 500, seed = 1)
eff <- resolvable_row_column_efficiency(design,
  ar = list(eta_r = 35, eta_c = 60, psi = 10, rho_r = 0.5, rho_c = 0.5))
eff[c("efficiency_rc", "efficiency_spatial", "ar_criterion")]
#> $efficiency_rc
#> [1] 0.5496625
#> 
#> $efficiency_spatial
#> [1] 0.8154889
#> 
#> $ar_criterion
#> [1] 4.338679
#> 
```
