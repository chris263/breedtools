# Load location coordinates

Loads location metadata from a `.xlsx`, `.xls`, or `.csv` file, or from
an already loaded data frame. The returned table has standardized
columns: `location`, `latitude`, and `longitude`.

## Usage

``` r
read_locations(
  file,
  location_col = NULL,
  latitude_col = NULL,
  longitude_col = NULL,
  phenotype_data = NULL,
  phenotype_location_col = NULL,
  match_by = "location",
  reference_locations = NULL,
  ignore_case = FALSE,
  sheet = NULL,
  ...
)
```

## Arguments

- file:

  Path to a location file, or an already loaded data frame.

- location_col:

  Optional location name column.

- latitude_col:

  Optional latitude column.

- longitude_col:

  Optional longitude column.

- phenotype_data:

  Optional raw phenotype data frame, prepared phenotype data frame, or
  result from
  [`run_environment_stratification()`](https://chris263.github.io/breedtools/reference/run_environment_stratification.md).
  If supplied, location names are extracted and checked against the
  location file.

- phenotype_location_col:

  Optional location column in `phenotype_data`.

- match_by:

  Column to extract when `phenotype_data` is a stratification result.
  Use `"location"`, `"environment"`, or `"environment_label"`.

- reference_locations:

  Optional character vector of expected location names. If supplied,
  every reference location must be present in `file`.

- ignore_case:

  Logical. If `TRUE`, location matching ignores case.

- sheet:

  Optional Excel sheet name or number passed to
  [`read_breedbase_file()`](https://chris263.github.io/breedtools/reference/read_breedbase_file.md).

- ...:

  Additional arguments passed to
  [`read_breedbase_file()`](https://chris263.github.io/breedtools/reference/read_breedbase_file.md).

## Value

A tibble with `location`, `latitude`, and `longitude`.

## Examples

``` r
locations <- data.frame(
  location = c("Loc1", "Loc2"),
  latitude = c(35.1, 36.2),
  longitude = c(-78.9, -79.4)
)

read_locations(locations)
#> # A tibble: 2 × 3
#>   location latitude longitude
#>   <chr>       <dbl>     <dbl>
#> 1 Loc1         35.1     -78.9
#> 2 Loc2         36.2     -79.4

# Breedbase-style location file, checked against the phenotype locations
loc_file <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
read_locations(loc_file, phenotype_data = pheno)
#> # A tibble: 4 × 3
#>   location latitude longitude
#>   <chr>       <dbl>     <dbl>
#> 1 Ithaca      42.4     -76.5 
#> 2 Geneva      42.9     -77.0 
#> 3 Ibadan       7.38      3.95
#> 4 Ubiaja       6.66      6.39
```
