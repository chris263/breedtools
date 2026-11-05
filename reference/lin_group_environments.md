# Group environments using Lin-style GE compatibility testing

Group environments using Lin-style GE compatibility testing

## Usage

``` r
lin_group_environments(data, alpha = 0.05)
```

## Arguments

- data:

  Standardized phenotype data.

- alpha:

  Significance level.

## Value

A list with pairwise tests, group summary, group membership, and
ungrouped environments.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

lin <- lin_group_environments(prepared$data, alpha = 0.05)
lin$group_summary
#>   group_id                                 environments n_env n_genotypes r_eff
#> 1  Group_1 Geneva_GEN_2025_yield, Ithaca_ITH_2025_yield     2          12     2
#> 2  Group_2 Ibadan_IBA_2025_yield, Ubiaja_UBI_2025_yield     2          12     2
#>      ss_ge df_ge    ms_ge message mse_error df_error   f_value    p_value
#> 1 24.84937    11 2.259034      OK  1.193709       22 1.8924496 0.09774358
#> 2 48.65081    11 4.422801      OK  5.615979       22 0.7875387 0.64998303
#>   compatible
#> 1       TRUE
#> 2       TRUE
lin$group_membership
#>   group_id           environment
#> 1  Group_1 Geneva_GEN_2025_yield
#> 2  Group_1 Ithaca_ITH_2025_yield
#> 3  Group_2 Ibadan_IBA_2025_yield
#> 4  Group_2 Ubiaja_UBI_2025_yield
```
