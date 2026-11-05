# Set the section label for the next checks

Checks recorded after this call with
[`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md)
and
[`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md)
are grouped under `name` in
[`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md).

## Usage

``` r
check_section(name)
```

## Arguments

- name:

  Section label, e.g. `"1. Row-column design"`. `NULL` removes the
  label.

## Value

Invisibly, the previous label.

## Examples

``` r
check_reset()
#> Ready to start: no earlier checks to clear. The checks you run from now on will be counted together.
check_section("arithmetic")
#> 
#> == arithmetic ==
check_result("one equals one", 1 == 1)
#>   [OK]   one equals one
check_summary()
#> 
#> 1 of 1 checks passed
#> 
#> By section:
#>  section    checks passed failed notes
#>  arithmetic 1      1      0      0    
#> 
#> All checks are correct.
check_reset()
#> Earlier checks cleared (1 check, which passed). Ready to start new checks.
```
