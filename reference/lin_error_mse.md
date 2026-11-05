# Calculate residual error MSE for a Lin environment group

Calculate residual error MSE for a Lin environment group

## Usage

``` r
lin_error_mse(data, envs, common_genotypes)
```

## Arguments

- data:

  Standardized phenotype data.

- envs:

  Character vector of environments.

- common_genotypes:

  Character vector of common genotypes.

## Value

A list with MSE, residual df, and message.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

envs <- prepared$env_info$environment[1:2]
genotypes <- unique(prepared$data$accession_name)
lin_error_mse(prepared$data, envs, genotypes)
#> $mse
#> [1] 3.79368
#> 
#> $df_error
#> [1] 22
#> 
#> $message
#> [1] "Detected randomized complete block layout; blocks are fitted within environment when block has at least two levels."
#> 
```
