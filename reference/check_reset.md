# Clear the recorded checks

Removes the checks recorded so far with
[`check_result()`](https://chris263.github.io/breedtools/reference/check_result.md)
and
[`check_design()`](https://chris263.github.io/breedtools/reference/check_design.md),
so that the next
[`check_summary()`](https://chris263.github.io/breedtools/reference/check_summary.md)
only counts new checks. It does not check or change any design or data;
it only empties the list. A message says how many earlier checks were
cleared and how many of them had passed.

## Usage

``` r
check_reset()
```

## Value

Invisibly, a tibble of the checks that were cleared.

## Examples

``` r
check_reset()
#> Ready to start: no earlier checks to clear. The checks you run from now on will be counted together.
check_result("one equals one", 1 == 1)
#>   [OK]   one equals one
check_reset()
#> Earlier checks cleared (1 check, which passed). Ready to start new checks.
```
