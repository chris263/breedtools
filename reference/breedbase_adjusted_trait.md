# Find or create the Breedbase trait for adjusted values

First-stage values need their own trait, so raw and adjusted values stay
apart. Breedbase makes such traits by composing the measured trait with
a statistic term, e.g. "fresh root yield" + "Adjusted Means from Spatial
Correction using SpATS R" gives the composed trait
`"fresh root yield|Adjusted Means from Spatial Correction using SpATS R"`
(ontology `COMP`). This function finds that trait, or creates it with
`create = TRUE` through the same calls as the "Compose a new trait" page
of Breedbase.

## Usage

``` r
breedbase_adjusted_trait(con, traits, statistic = NULL, create = FALSE)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- traits:

  Measured trait names, as in Breedbase.

- statistic:

  The statistic term, with or without its id, as listed by
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md).
  By default it is chosen from this Breedbase's own list of statistics:
  the spatial correction term (SpATS). It is always checked against the
  server.

- create:

  Logical. Create the traits that do not exist yet.

## Value

A tibble with `trait`, `adjusted_trait`, `adjusted_trait_db_id` and
`status` (`"exists"`, `"created"` or `"missing"`).

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
breedbase_adjusted_trait(con, "fresh root yield")                  # look only
breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)   # create if needed
} # }
```
