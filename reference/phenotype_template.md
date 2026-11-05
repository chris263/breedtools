# Phenotyping template with the Breedbase plot names of a design

Makes a table with one row per plot, the plot names that
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
gives (or gave) to the plots, the design columns (accession, block, rep,
row, column) and one empty column per trait. Fill it in (in R or in a
spreadsheet) and upload it with
[`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md).

## Usage

``` r
phenotype_template(
  design,
  trial_name,
  traits = character(0),
  plot_name_pattern = "{trial}_{accession}_{block}_{rep}",
  plot_names = NULL
)
```

## Arguments

- design:

  A design from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md),
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
  or
  [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md).

- trial_name:

  The trial name used when saving the design. For p-rep designs saved as
  one trial per location, the location is added as in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md).
  Can be omitted for a prepared Breedbase trial.

- traits:

  Trait names as in Breedbase, e.g. `"fresh root yield"` or
  `"fresh root yield|CO_334:0000013"`. Can be empty (no trait columns).

- plot_name_pattern, plot_names:

  The same plot naming used in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md).

## Value

A tibble: `plot_name`, `trial_name`, `accession`, `plot_number`,
`block`, `rep`, `row`, `col`, `is_check` and one column per trait.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
sheet <- phenotype_template(rc, "AYT_2026_Ibadan",
                            traits = c("fresh root yield", "dry matter content percentage"))
head(sheet)
#> # A tibble: 6 × 11
#>   plot_name    trial_name accession plot_number block   rep   row   col is_check
#>   <chr>        <chr>      <chr>           <int> <int> <int> <int> <int> <lgl>   
#> 1 AYT_2026_Ib… AYT_2026_… G51                 1     1     1     1     1 FALSE   
#> 2 AYT_2026_Ib… AYT_2026_… G64                 2     1     1     1     2 FALSE   
#> 3 AYT_2026_Ib… AYT_2026_… G97                 3     1     1     1     3 FALSE   
#> 4 AYT_2026_Ib… AYT_2026_… C3                  4     1     1     1     4 TRUE    
#> 5 AYT_2026_Ib… AYT_2026_… G8                  5     2     1     1     5 FALSE   
#> 6 AYT_2026_Ib… AYT_2026_… C2                  6     2     2     1     6 TRUE    
#> # ℹ 2 more variables: `fresh root yield` <dbl>,
#> #   `dry matter content percentage` <dbl>
```
