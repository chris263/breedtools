# Record and print a check

Prints `[OK]` or `[FAIL]` followed by a description, and records the
result so that
[`check_summary()`](https://chris263.github.io/breedtools/reference/check_summary.md)
can report how many checks passed. Use it to write your own checks while
exploring a design or an analysis.

## Usage

``` r
check_result(what, ok, section = NULL, note = FALSE)
```

## Arguments

- what:

  Description of the check, e.g. `"every entry planted once"`.

- ok:

  Logical. `TRUE` when the check passes; anything else (including `NA`)
  counts as a failure.

- section:

  Optional label for grouping checks, e.g. `"p-rep"`. Defaults to the
  label set with
  [`check_section()`](https://chris263.github.io/breedtools/reference/check_section.md).

- note:

  Logical. If `TRUE`, the check is informational: it prints `[NOTE]`
  instead of `[FAIL]` when `ok` is not `TRUE` and is not counted as a
  failure.

## Value

Invisibly, `TRUE` if the check passed.

## See also

[`check_section()`](https://chris263.github.io/breedtools/reference/check_section.md),
[`check_summary()`](https://chris263.github.io/breedtools/reference/check_summary.md),
[`check_reset()`](https://chris263.github.io/breedtools/reference/check_reset.md),
[`check_design()`](https://chris263.github.io/breedtools/reference/check_design.md)

## Examples

``` r
check_reset()
#> Ready to start: no earlier checks to clear. The checks you run from now on will be counted together.
check_result("two plus two is four", 2 + 2 == 4)
#>   [OK]   two plus two is four
check_result("a deliberate failure", FALSE, section = "demo")
#>   [FAIL] a deliberate failure
check_summary()
#> 
#> 1 of 2 checks passed
#> 
#> By section:
#>  section      checks passed failed notes
#>  (no section) 1      1      0      0    
#>  demo         1      0      1      0    
#> 
#> Failed checks (not correct):
#>  section check               
#>  demo    a deliberate failure
check_reset()
#> Earlier checks cleared (2 checks: 1 passed, 1 failed). Ready to start new checks.
```
