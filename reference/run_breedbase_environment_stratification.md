# Run Lin-style environment stratification from a Breedbase phenotype file

Reads a Breedbase phenotype file, prepares the phenotype data, detects
the experimental design, calculates ANOVA, and groups environments using
a Lin-style genotype-by-environment interaction test.

## Usage

``` r
run_breedbase_environment_stratification(
  file,
  trait_col,
  alpha = 0.05,
  normalize_names = TRUE
)
```

## Arguments

- file:

  Path to a `.xlsx`, `.xls`, or `.csv` phenotype file, or an already
  loaded data frame.

- trait_col:

  Trait column name.

- alpha:

  Significance level. Default is `0.05`.

- normalize_names:

  Logical. If `TRUE`, normalize Breedbase-style column names.

## Value

A list containing standardized phenotype data, environment metadata,
ANOVA results, pairwise environment tests, group summary, group
membership, ungrouped environments, summary, and message.

## Examples

``` r
pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")

result <- run_breedbase_environment_stratification(
  file = pheno_file,
  trait_col = "fresh root yield",
  alpha = 0.05
)

result$group_summary
#>   group_id                                     environments      locations
#> 1  Group_1 Geneva / GEN_2025_yield, Ithaca / ITH_2025_yield Geneva, Ithaca
#> 2  Group_2 Ibadan / IBA_2025_yield, Ubiaja / UBI_2025_yield Ibadan, Ubiaja
#>                           trials years n_env n_genotypes r_eff    ss_ge df_ge
#> 1 GEN_2025_yield, ITH_2025_yield           2          12     2 24.84937    11
#> 2 IBA_2025_yield, UBI_2025_yield           2          12     2 48.65081    11
#>      ms_ge message mse_error df_error   f_value    p_value compatible
#> 1 2.259034      OK  1.193709       22 1.8924496 0.09774358       TRUE
#> 2 4.422801      OK  5.615979       22 0.7875387 0.64998303       TRUE
result$group_membership
#>   group_id           environment location          trial year
#> 1  Group_1 Geneva_GEN_2025_yield   Geneva GEN_2025_yield     
#> 2  Group_1 Ithaca_ITH_2025_yield   Ithaca ITH_2025_yield     
#> 3  Group_2 Ibadan_IBA_2025_yield   Ibadan IBA_2025_yield     
#> 4  Group_2 Ubiaja_UBI_2025_yield   Ubiaja UBI_2025_yield     
#>         environment_label
#> 1 Geneva / GEN_2025_yield
#> 2 Ithaca / ITH_2025_yield
#> 3 Ibadan / IBA_2025_yield
#> 4 Ubiaja / UBI_2025_yield
```
