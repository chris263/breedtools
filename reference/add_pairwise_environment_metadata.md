# Add environment metadata to pairwise Lin results

Add environment metadata to pairwise Lin results

## Usage

``` r
add_pairwise_environment_metadata(pairwise, env_info)
```

## Arguments

- pairwise:

  Pairwise result table.

- env_info:

  Environment metadata table.

## Value

Pairwise table with metadata.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

lin <- lin_group_environments(prepared$data)
add_pairwise_environment_metadata(lin$pairwise, prepared$env_info)
#>   env1_location     env1_trial env1_year               env1_name env2_location
#> 1        Geneva GEN_2025_yield           Geneva / GEN_2025_yield        Ithaca
#> 2        Ibadan IBA_2025_yield           Ibadan / IBA_2025_yield        Ubiaja
#> 3        Ibadan IBA_2025_yield           Ibadan / IBA_2025_yield        Ithaca
#> 4        Geneva GEN_2025_yield           Geneva / GEN_2025_yield        Ibadan
#> 5        Ithaca ITH_2025_yield           Ithaca / ITH_2025_yield        Ubiaja
#> 6        Geneva GEN_2025_yield           Geneva / GEN_2025_yield        Ubiaja
#>       env2_trial env2_year               env2_name n_env n_genotypes r_eff
#> 1 ITH_2025_yield           Ithaca / ITH_2025_yield     2          12     2
#> 2 UBI_2025_yield           Ubiaja / UBI_2025_yield     2          12     2
#> 3 ITH_2025_yield           Ithaca / ITH_2025_yield     2          12     2
#> 4 IBA_2025_yield           Ibadan / IBA_2025_yield     2          12     2
#> 5 UBI_2025_yield           Ubiaja / UBI_2025_yield     2          12     2
#> 6 UBI_2025_yield           Ubiaja / UBI_2025_yield     2          12     2
#>       ss_ge df_ge     ms_ge message mse_error df_error    f_value      p_value
#> 1  24.84937    11  2.259034      OK  1.193709       22  1.8924496 9.774358e-02
#> 2  48.65081    11  4.422801      OK  5.615979       22  0.7875387 6.499830e-01
#> 3 400.64496    11 36.422269      OK  3.847724       22  9.4659254 5.168582e-06
#> 4 419.74277    11 38.158433      OK  3.793680       22 10.0584219 3.095762e-06
#> 5 455.49937    11 41.409034      OK  3.016008       22 13.7297514 1.998853e-07
#> 6 466.14102    11 42.376457      OK  2.961963       22 14.3068803 1.373132e-07
#>   compatible
#> 1       TRUE
#> 2       TRUE
#> 3      FALSE
#> 4      FALSE
#> 5      FALSE
#> 6      FALSE
```
