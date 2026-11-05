# Calculate design-aware ANOVA

Calculate design-aware ANOVA

## Usage

``` r
calculate_environment_stratification_anova(data)
```

## Arguments

- data:

  Standardized phenotype data.

## Value

ANOVA table as a data frame.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

calculate_environment_stratification_anova(prepared$data)
#>   design                       term df     sum_sq    mean_sq    f_value
#> 1   RCBD                environment  3 1220.69248 406.897492 119.505481
#> 2   RCBD             accession_name 11  509.56706  46.324278  13.605405
#> 3   RCBD   environment:block_number  4   22.10538   5.526344   1.623083
#> 4   RCBD environment:accession_name 33  907.76415  27.508005   8.079080
#> 5   RCBD                  Residuals 44  149.81312   3.404844         NA
#>        p_value
#> 1 3.611560e-21
#> 2 8.173352e-11
#> 3 1.853783e-01
#> 4 3.241667e-10
#> 5           NA
#>                                                                                                               message
#> 1 Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels.
#> 2 Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels.
#> 3 Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels.
#> 4 Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels.
#> 5 Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels.
```
