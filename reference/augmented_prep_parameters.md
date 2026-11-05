# Feasible parameters for augmented p-rep constructions

Lists the parameter sets that
[`augmented_prep_design()`](https://chris263.github.io/breedtools/reference/augmented_prep_design.md)
can build for a given number of entries and locations, in the spirit of
Table 3 of Williams, Piepho and Whitaker (2011). For each number `m` of
p-rep entries per block, it gives the number of blocks per replicate
`s`, the final block size(s), the number of entries duplicated at each
location, and how many "dummy" entries are dropped when the number of
entries is not a multiple of the array size (dropping them makes some
blocks one plot smaller).

## Usage

``` r
augmented_prep_parameters(n_treatments, n_locations, max_block_size = 20)
```

## Arguments

- n_treatments:

  Number of entries (`v`).

- n_locations:

  Number of locations (`l`).

- max_block_size:

  Largest final block size to list.

## Value

A tibble with one row per method and `m`: `method`, `m`, `s`,
`block_size` (largest final block), `min_block_size`,
`blocks_per_location`, `plots_per_location` (largest),
`duplicated_per_location` (largest), `dummies`, `feasible`, and `note`.

## References

Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented p-rep
designs. Biometrical Journal 53, 19-27.
[doi:10.1002/bimj.201000102](https://doi.org/10.1002/bimj.201000102)

## Examples

``` r
augmented_prep_parameters(n_treatments = 120, n_locations = 6, max_block_size = 10)
#> # A tibble: 2 × 11
#>   method     m     s block_size min_block_size blocks_per_location
#>   <chr>  <int> <int>      <int>          <int>               <int>
#> 1 alpha      2    10          7              7                  20
#> 2 alpha2     2    10          7              7                  20
#> # ℹ 5 more variables: plots_per_location <int>, duplicated_per_location <int>,
#> #   dummies <int>, feasible <lgl>, note <chr>
```
