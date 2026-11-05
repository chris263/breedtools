# Run Lin-style environment stratification

This function reads phenotype data, detects design structure, calculates
ANOVA, tests genotype-by-environment interaction for environment groups,
and returns compatible environment groups.

## Usage

``` r
run_environment_stratification(
  phenotype_file,
  trait,
  alpha = 0.05,
  sep = NULL,
  normalize_names = TRUE,
  export_json = FALSE,
  output_dir = NULL,
  prefix = "environment_stratification"
)
```

## Arguments

- phenotype_file:

  Path to a `.xlsx`, `.xls`, or `.csv` phenotype file, or an already
  loaded data frame.

- trait:

  Trait column name.

- alpha:

  Significance level. Default is `0.05`.

- sep:

  Deprecated. File delimiters are detected from the extension by
  [`read_breedbase_file()`](https://chris263.github.io/breedtools/reference/read_breedbase_file.md).

- normalize_names:

  Logical. If `TRUE`, normalize Breedbase column names.

- export_json:

  Logical. If `TRUE`, writes a JSON result file to `output_dir`.

- output_dir:

  Output directory for the JSON file. Defaults to `NULL`, in which case
  nothing is written; it must be supplied when `export_json = TRUE`.

- prefix:

  Prefix for JSON output files.

## Value

A list with summary, ANOVA results, pairwise tests, group summary, group
membership, ungrouped environments, and message.

## Examples

``` r
pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")

result <- run_environment_stratification(
  phenotype_file = pheno_file,
  trait = "fresh root yield",
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

# Export the result as JSON to a temporary directory
out_dir <- file.path(tempdir(), "stratification")
result <- run_environment_stratification(
  phenotype_file = pheno_file,
  trait = "fresh root yield",
  export_json = TRUE,
  output_dir = out_dir
)
list.files(out_dir)
#> [1] "environment_stratification.json"
unlink(out_dir, recursive = TRUE)
```
