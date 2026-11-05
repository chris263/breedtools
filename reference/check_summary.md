# Summary of the recorded checks

Reports what was checked and what is correct: how many of the checks
recorded with
[`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md)
and
[`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md)
passed, broken down by section, followed by the checks that failed and
the notes.

## Usage

``` r
check_summary(details = FALSE, print = TRUE)
```

## Arguments

- details:

  Logical. If `TRUE`, also list every check with its status.

- print:

  Logical. If `TRUE` (default), print the report.

## Value

Invisibly, a tibble of all recorded checks with columns `section`,
`check` and `status` (`"OK"`, `"FAIL"` or `"NOTE"`).

## Examples

``` r
check_reset()
#> Earlier checks cleared (4 checks, all passed). Ready to start new checks.
check_result("one equals one", 1 == 1, section = "arithmetic")
#>   [OK]   one equals one
check_result("two is even", 2 %% 2 == 0, section = "arithmetic")
#>   [OK]   two is even
check_summary()
#> 
#> 2 of 2 checks passed
#> 
#> By section:
#>  section    checks passed failed notes
#>  arithmetic 2      2      0      0    
#> 
#> All checks are correct.
check_summary(details = TRUE)
#> 
#> 2 of 2 checks passed
#> 
#> By section:
#>  section    checks passed failed notes
#>  arithmetic 2      2      0      0    
#> 
#> All checks are correct.
#> 
#> Every check:
#>  status section    check         
#>  OK     arithmetic one equals one
#>  OK     arithmetic two is even   
check_reset()
#> Earlier checks cleared (2 checks, all passed). Ready to start new checks.
```
