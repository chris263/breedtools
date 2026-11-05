# Read a Breedbase Excel phenotype file

`read_breedbase_excel()` is deprecated. Use
[`read_breedbase_file()`](https://chris263.github.io/breedtools/reference/read_breedbase_file.md)
for `.xlsx`, `.xls`, and `.csv` files.

## Usage

``` r
read_breedbase_excel(file, ...)
```

## Arguments

- file:

  Path to the Excel file.

- ...:

  Additional arguments passed to
  [`read_breedbase_file()`](https://chris263.github.io/breedtools/reference/read_breedbase_file.md).

## Value

A tibble containing combined data from all non-empty sheets.

## Examples

``` r
xlsx_file <- system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools")

# Deprecated (emits a deprecation warning): use read_breedbase_file()
pheno <- suppressWarnings(read_breedbase_excel(xlsx_file))

# Preferred
pheno <- read_breedbase_file(xlsx_file)
head(pheno)
#> # A tibble: 6 × 15
#>   studyYear programName     studyName     studyDesign locationName germplasmName
#>       <dbl> <chr>           <chr>         <chr>       <chr>        <chr>        
#> 1      2025 Example program ITH_2025_yie… RCBD        Ithaca       G12          
#> 2      2025 Example program ITH_2025_yie… RCBD        Ithaca       G05          
#> 3      2025 Example program ITH_2025_yie… RCBD        Ithaca       G09          
#> 4      2025 Example program ITH_2025_yie… RCBD        Ithaca       G02          
#> 5      2025 Example program ITH_2025_yie… RCBD        Ithaca       G01          
#> 6      2025 Example program ITH_2025_yie… RCBD        Ithaca       G07          
#> # ℹ 9 more variables: observationUnitName <chr>, replicate <dbl>,
#> #   blockNumber <dbl>, plotNumber <dbl>, rowNumber <dbl>, colNumber <dbl>,
#> #   `fresh root yield|CO_334:0000013` <dbl>,
#> #   `plant height|CO_334:0000018` <dbl>, .sheet <chr>
```
