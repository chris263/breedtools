# Find the first matching column

Find the first matching column

## Usage

``` r
find_column(data, candidates, label = "column")
```

## Arguments

- data:

  A data frame.

- candidates:

  Candidate column names.

- label:

  Label used in the error message.

## Value

The first matching column name.

## Examples

``` r
df <- data.frame(locationName = "Loc1", germplasmName = "G1")
find_column(df, c("germplasmName", "accession_name"), "accession")
#> [1] "germplasmName"
```
