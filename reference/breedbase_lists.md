# List breeding programs, folders and locations in Breedbase

Use these to find the exact names to give to
[`save_breedbase()`](https://chris263.github.io/breedtools/reference/save_breedbase.md).

## Usage

``` r
breedbase_programs(con)

breedbase_folders(con, program)

breedbase_locations(con)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://chris263.github.io/breedtools/reference/breedbase_login.md).

- program:

  For `breedbase_folders()`: the breeding program name.

## Value

A tibble with names and database ids.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
breedbase_programs(con)
breedbase_folders(con, program = "IITA")
breedbase_locations(con)
} # }
```
