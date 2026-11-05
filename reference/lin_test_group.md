# Test one Lin environment group

Test one Lin environment group

## Usage

``` r
lin_test_group(data, envs, alpha = 0.05)
```

## Arguments

- data:

  Standardized phenotype data.

- envs:

  Character vector of environments.

- alpha:

  Significance level.

## Value

Data frame with Lin GE test result.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

envs <- prepared$env_info$environment[1:2]
lin_test_group(prepared$data, envs, alpha = 0.05)
#>                                   environments n_env n_genotypes r_eff    ss_ge
#> 1 Geneva_GEN_2025_yield, Ibadan_IBA_2025_yield     2          12     2 419.7428
#>   df_ge    ms_ge message mse_error df_error  f_value      p_value compatible
#> 1    11 38.15843      OK   3.79368       22 10.05842 3.095762e-06      FALSE
```
