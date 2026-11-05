# Guess a column name from possible Breedbase column names

Guess a column name from possible Breedbase column names

## Usage

``` r
guess_column(df, possible_names)
```

## Arguments

- df:

  A data frame.

- possible_names:

  A character vector of possible column names.

## Value

A column name if found, otherwise `NULL`.

## Examples

``` r
df <- data.frame(germplasmName = "G1", phenotype = 10)
guess_column(df, c("germplasmName", "accession_name"))
#> [1] "germplasmName"
```
