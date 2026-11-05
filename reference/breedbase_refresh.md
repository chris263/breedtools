# Refresh the Breedbase phenotype cache

Values saved through BrAPI
([`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md),
[`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md))
are in the database at once (the BrAPI search finds them), but the trial
pages and the search wizard of Breedbase read a cache ("materialized
views") that the BrAPI upload does not refresh. This asks Breedbase to
refresh it, like the "Update" button of the search wizard. The save
functions call it after saving. If it cannot be reached, the values are
still saved: they appear after the next refresh of the site (usually
overnight), or when you click "Update" on the Search Wizard page.

## Usage

``` r
breedbase_refresh(con, views = c("phenotypes", "fullview"), quiet = FALSE)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- views:

  Which cache to refresh: `"phenotypes"` (default) or `"fullview"`.

- quiet:

  Logical. If `FALSE`, print what happened.

## Value

Invisibly, `TRUE` if Breedbase accepted the refresh, else `FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
breedbase_refresh(con)
} # }
```
