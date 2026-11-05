# Prepare Breedbase phenotype data for multi-environment analysis

This function standardizes Breedbase phenotype data into five columns:
`environment`, `genotype`, `replicate`, `block`, and `value`.

## Usage

``` r
prepare_breedbase_pheno(
  pheno_raw,
  trait_col,
  location_col = NULL,
  genotype_col = NULL,
  rep_col = NULL,
  block_col = NULL
)
```

## Arguments

- pheno_raw:

  A raw Breedbase phenotype data frame, usually read from a Breedbase
  Excel file.

- trait_col:

  Name of the trait column to analyze.

- location_col:

  Optional location/environment column name. If `NULL`, the function
  tries to guess it.

- genotype_col:

  Optional genotype/accession column name. If `NULL`, the function tries
  to guess it.

- rep_col:

  Optional replicate column name.

- block_col:

  Optional block column name.

## Value

A tibble with standardized phenotype data.

## Examples

``` r
df <- data.frame(
  location = c("Loc1", "Loc1", "Loc2", "Loc2"),
  germplasmName = c("G1", "G2", "G1", "G2"),
  rep = c(1, 1, 1, 1),
  phenotype = c(10, 12, 15, 14)
)

prepare_breedbase_pheno(df, trait_col = "phenotype")
#> # A tibble: 4 × 5
#>   environment genotype replicate block value
#>   <fct>       <fct>    <fct>     <fct> <dbl>
#> 1 Loc1        G1       1         1        10
#> 2 Loc1        G2       1         1        12
#> 3 Loc2        G1       1         1        15
#> 4 Loc2        G2       1         1        14
```
