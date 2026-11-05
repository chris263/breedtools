# Prepare phenotype data for Lin environment stratification

Reads a Breedbase-style phenotype data frame and standardizes it for
environment stratification.

## Usage

``` r
prepare_environment_stratification_data(pheno, trait, normalize_names = TRUE)
```

## Arguments

- pheno:

  A phenotype data frame.

- trait:

  Trait column name.

- normalize_names:

  Logical. If `TRUE`, normalizes column names.

## Value

A list with standardized phenotype data, environment metadata, and
selected trait column.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

head(prepared$data)
#>             environment location          trial year study_design
#> 1 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#> 2 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#> 3 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#> 4 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#> 5 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#> 6 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield              RCBD
#>   accession_name rep_number block_number row_number col_number phenotype
#> 1            G12          1            1          1          1     26.44
#> 2            G05          1            1          1          2     25.17
#> 3            G09          1            1          1          3     29.85
#> 4            G02          1            1          1          4     25.68
#> 5            G01          1            1          1          5     31.21
#> 6            G07          1            1          1          6     27.84
prepared$env_info
#>             environment location          trial year       environment_label
#> 1 Geneva_GEN_2025_yield   Geneva GEN_2025_yield      Geneva / GEN_2025_yield
#> 2 Ibadan_IBA_2025_yield   Ibadan IBA_2025_yield      Ibadan / IBA_2025_yield
#> 3 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield      Ithaca / ITH_2025_yield
#> 4 Ubiaja_UBI_2025_yield   Ubiaja UBI_2025_yield      Ubiaja / UBI_2025_yield
prepared$trait_col
#> [1] "fresh root yield"
```
