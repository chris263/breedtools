# Normalize Breedbase-style column names

Removes ontology suffixes such as `.CO...` or `|CO_...` and replaces
dots with spaces.

## Usage

``` r
normalize_column_name(x)
```

## Arguments

- x:

  Character vector of column names.

## Value

Normalized character vector.

## Examples

``` r
normalize_column_name(c(
  "fresh root yield|CO_334:0000013",
  "plant.height.CO_334.0000018",
  "germplasmName"
))
#> [1] "fresh root yield" "plant height"     "germplasmName"   
```
