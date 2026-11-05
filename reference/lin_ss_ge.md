# Calculate Lin genotype-by-environment sum of squares

Calculate Lin genotype-by-environment sum of squares

## Usage

``` r
lin_ss_ge(data, envs)
```

## Arguments

- data:

  Standardized phenotype data.

- envs:

  Character vector of environments to test.

## Value

A list with a summary table and common genotype names.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

envs <- prepared$env_info$environment[1:2]
lin_ss_ge(prepared$data, envs)$summary
#>                                   environments n_env n_genotypes r_eff    ss_ge
#> 1 Geneva_GEN_2025_yield, Ibadan_IBA_2025_yield     2          12     2 419.7428
#>   df_ge    ms_ge message
#> 1    11 38.15843      OK
```
