# Add environment metadata to result table

Add environment metadata to result table

## Usage

``` r
add_environment_metadata(results, env_info, env_col = "environment")
```

## Arguments

- results:

  Result data frame.

- env_info:

  Environment metadata table.

- env_col:

  Name of the environment column in `results`.

## Value

Joined data frame.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

membership <- data.frame(
  group_id = "Group_1",
  environment = prepared$env_info$environment[1:2]
)
add_environment_metadata(membership, prepared$env_info)
#>   group_id           environment location          trial year
#> 1  Group_1 Geneva_GEN_2025_yield   Geneva GEN_2025_yield     
#> 2  Group_1 Ibadan_IBA_2025_yield   Ibadan IBA_2025_yield     
#>         environment_label
#> 1 Geneva / GEN_2025_yield
#> 2 Ibadan / IBA_2025_yield
```
