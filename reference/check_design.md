# Standard checks of a design

Runs the structural checks that every design of its kind should pass and
prints one `[OK]` / `[FAIL]` line per check (`[NOTE]` for informational
ones). The results are recorded for
[`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md).

## Usage

``` r
check_design(design, section = NULL)
```

## Arguments

- design:

  A design object from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
  or
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md).

- section:

  Optional label for the recorded checks. Defaults to the label set with
  [`check_section()`](https://solgenomics.github.io/breedtools/reference/check_section.md),
  or else the kind of design.

## Value

Invisibly, `TRUE` if all (non-note) checks passed.

## Details

- Augmented row-column designs:

  every plot used once, every entry planted once, every field row and
  column has a check, the requested number of check plots in each block,
  positive error df, all row, column and block effects estimable from
  the checks, and no check plots next to each other inside a block (a
  note for neighbours across block boundaries).

- Augmented p-rep designs:

  every entry at every location, no entry twice in a block, block sizes
  within a location differ by at most one, every location connected, the
  pairwise concurrence (a note when pairs share more than one block),
  and the replicate structure of resolvable designs.

- Resolvable row-column designs:

  every treatment once in every replicate, every plot used once, the
  spatial efficiency below its upper bound, and latinization when
  requested.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
check_reset()
#> Ready to start: no earlier checks to clear. The checks you run from now on will be counted together.
check_design(rc)
#> Checking the augmented row-column design
#>   [OK]   every plot used exactly once
#>   [OK]   every entry planted exactly once
#>   [OK]   every field row has a check
#>   [OK]   every field column has a check
#>   [OK]   check plots per block as requested
#>   [OK]   error df > 0 (6)
#>   [OK]   row, column and block effects estimable from the checks
#>   [OK]   checks spread evenly in each block (every check in full blocks)
#>   [OK]   even numbers of check plots in the rows and columns of each group
#>   [OK]   every check in each row and column with room for all checks
#>   [OK]   no checks next to each other inside a block (0 pairs)
#>   [OK]   no checks next to each other across block boundaries (0 pairs)
#>   [OK]   no check twice in a field row or column
check_summary()
#> 
#> 13 of 13 checks passed
#> 
#> By section:
#>  section                     checks passed failed notes
#>  augmented row-column design 13     13     0      0    
#> 
#> All checks are correct.
check_reset()
#> Earlier checks cleared (13 checks, all passed). Ready to start new checks.
```
