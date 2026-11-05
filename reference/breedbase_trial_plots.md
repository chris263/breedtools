# Read the plots of a trial saved in Breedbase

Returns the plots of a Breedbase trial (through BrAPI): plot name,
number, accession, rep, block, row and column. Use it to check that a
design in R is the one saved in Breedbase (see `breedbase_plots` in
[`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)),
or to work with a trial made elsewhere.

## Usage

``` r
breedbase_trial_plots(con, trial_name)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- trial_name:

  The trial name in Breedbase.

## Value

A tibble with `plot_name`, `plot_number`, `accession`, `rep`, `block`,
`row`, `col` and `is_check`, ordered by plot number.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
plots <- breedbase_trial_plots(con, "AYT_2026_Ibadan")
} # }
```
