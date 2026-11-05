# Prepare design factors

Prepare design factors

## Usage

``` r
prepare_design_factors(data)
```

## Arguments

- data:

  Standardized phenotype data.

## Value

Data frame with design columns converted to factors.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

str(prepare_design_factors(prepared$data))
#> 'data.frame':    96 obs. of  11 variables:
#>  $ environment   : Factor w/ 4 levels "Geneva_GEN_2025_yield",..: 3 3 3 3 3 3 3 3 3 3 ...
#>  $ location      : chr  "Ithaca" "Ithaca" "Ithaca" "Ithaca" ...
#>  $ trial         : chr  "ITH_2025_yield" "ITH_2025_yield" "ITH_2025_yield" "ITH_2025_yield" ...
#>  $ year          : chr  "" "" "" "" ...
#>  $ study_design  : chr  "RCBD" "RCBD" "RCBD" "RCBD" ...
#>  $ accession_name: Factor w/ 12 levels "G01","G02","G03",..: 12 5 9 2 1 7 3 10 11 8 ...
#>  $ rep_number    : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ block_number  : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ row_number    : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ col_number    : Factor w/ 12 levels "1","10","11",..: 1 5 6 7 8 9 10 11 12 2 ...
#>  $ phenotype     : num  26.4 25.2 29.9 25.7 31.2 ...
```
