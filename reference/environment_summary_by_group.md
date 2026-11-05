# Summarize environment names by group

Summarize environment names by group

## Usage

``` r
environment_summary_by_group(group_summary, group_membership)
```

## Arguments

- group_summary:

  Group summary table.

- group_membership:

  Group membership table.

## Value

Group summary with display metadata.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

lin <- lin_group_environments(prepared$data)
membership <- add_environment_metadata(lin$group_membership, prepared$env_info)
environment_summary_by_group(lin$group_summary, membership)
#>   group_id                                     environments      locations
#> 1  Group_1 Geneva / GEN_2025_yield, Ithaca / ITH_2025_yield Geneva, Ithaca
#> 2  Group_2 Ibadan / IBA_2025_yield, Ubiaja / UBI_2025_yield Ibadan, Ubiaja
#>                           trials years n_env n_genotypes r_eff    ss_ge df_ge
#> 1 GEN_2025_yield, ITH_2025_yield           2          12     2 24.84937    11
#> 2 IBA_2025_yield, UBI_2025_yield           2          12     2 48.65081    11
#>      ms_ge message mse_error df_error   f_value    p_value compatible
#> 1 2.259034      OK  1.193709       22 1.8924496 0.09774358       TRUE
#> 2 4.422801      OK  5.615979       22 0.7875387 0.64998303       TRUE
```
