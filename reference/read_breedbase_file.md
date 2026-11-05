# Read phenotype or design data from Excel, CSV, or tab-delimited files

Loads phenotype or field-design data from `.xlsx`, `.xls`, `.csv`,
`.tsv`, `.txt`, or `.design` files. Excel files are read from all
non-empty sheets and combined into one tibble with a `.sheet` column
identifying the source sheet. CSV files are comma-delimited; `.tsv`,
`.txt`, and `.design` files are tab-delimited.

## Usage

``` r
read_breedbase_file(file, sheet = NULL, ...)
```

## Arguments

- file:

  Path to the phenotype file.

- sheet:

  Optional Excel sheet name or number. If `NULL`, all non-empty Excel
  sheets are combined.

- ...:

  Additional arguments passed to
  [`readxl::read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html)
  for Excel files,
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) for CSV
  files, or
  [`utils::read.table()`](https://rdrr.io/r/utils/read.table.html) for
  tab-delimited files.

## Value

A tibble containing the loaded phenotype data.

## Examples

``` r
# CSV phenotype download
pheno_csv <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
head(pheno_csv)
#> # A tibble: 6 × 14
#>   studyYear programName     studyName     studyDesign locationName germplasmName
#>       <int> <chr>           <chr>         <chr>       <chr>        <chr>        
#> 1      2025 Example program ITH_2025_yie… RCBD        Ithaca       G12          
#> 2      2025 Example program ITH_2025_yie… RCBD        Ithaca       G05          
#> 3      2025 Example program ITH_2025_yie… RCBD        Ithaca       G09          
#> 4      2025 Example program ITH_2025_yie… RCBD        Ithaca       G02          
#> 5      2025 Example program ITH_2025_yie… RCBD        Ithaca       G01          
#> 6      2025 Example program ITH_2025_yie… RCBD        Ithaca       G07          
#> # ℹ 8 more variables: observationUnitName <chr>, replicate <int>,
#> #   blockNumber <int>, plotNumber <int>, rowNumber <int>, colNumber <int>,
#> #   `fresh root yield|CO_334:0000013` <dbl>,
#> #   `plant height|CO_334:0000018` <dbl>

# Excel workbook: all non-empty sheets are combined
pheno_xlsx <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools")
)
table(pheno_xlsx$.sheet)
#> 
#> GEN_2025_yield ITH_2025_yield 
#>             24             24 

# A single sheet
read_breedbase_file(
  system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools"),
  sheet = 1
)
#> # A tibble: 24 × 14
#>    studyYear programName     studyName    studyDesign locationName germplasmName
#>        <dbl> <chr>           <chr>        <chr>       <chr>        <chr>        
#>  1      2025 Example program ITH_2025_yi… RCBD        Ithaca       G12          
#>  2      2025 Example program ITH_2025_yi… RCBD        Ithaca       G05          
#>  3      2025 Example program ITH_2025_yi… RCBD        Ithaca       G09          
#>  4      2025 Example program ITH_2025_yi… RCBD        Ithaca       G02          
#>  5      2025 Example program ITH_2025_yi… RCBD        Ithaca       G01          
#>  6      2025 Example program ITH_2025_yi… RCBD        Ithaca       G07          
#>  7      2025 Example program ITH_2025_yi… RCBD        Ithaca       G03          
#>  8      2025 Example program ITH_2025_yi… RCBD        Ithaca       G10          
#>  9      2025 Example program ITH_2025_yi… RCBD        Ithaca       G11          
#> 10      2025 Example program ITH_2025_yi… RCBD        Ithaca       G08          
#> # ℹ 14 more rows
#> # ℹ 8 more variables: observationUnitName <chr>, replicate <dbl>,
#> #   blockNumber <dbl>, plotNumber <dbl>, rowNumber <dbl>, colNumber <dbl>,
#> #   `fresh root yield|CO_334:0000013` <dbl>,
#> #   `plant height|CO_334:0000018` <dbl>
```
