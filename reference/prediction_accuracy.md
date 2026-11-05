# Accuracy of adjusted means against known values

In a simulation the true genotype values are known. This function
measures how well an analysis recovered them: the correlation between
the adjusted means and the true values.

## Usage

``` r
prediction_accuracy(
  fit,
  truth,
  genotypes = NULL,
  method = c("pearson", "spearman")
)
```

## Arguments

- fit:

  An analysis object (from any `analyze_*()` function) or a data frame
  with `genotype` and `adjusted_mean` columns.

- truth:

  Named numeric vector of true genotype values; names are the genotype
  names.

- genotypes:

  Optional subset of genotypes to use (for example only the test
  entries, without the checks).

- method:

  Correlation method: `"pearson"` (default) or `"spearman"` (rank
  correlation, i.e. how well the ranking is recovered).

## Value

A single number between -1 and 1.

## Examples

``` r
design <- augmented_prep_design(paste0("G", 1:30), c("A", "B"), block_size = 6, seed = 1)
trial <- design$design
set.seed(1)
truth <- setNames(rnorm(30, sd = 3), paste0("G", 1:30))
trial$yield <- 50 + truth[trial$all_entries] + rnorm(nrow(trial))
fit <- analyze_augmented_prep_design(trial, "yield")
prediction_accuracy(fit, truth)
#> [1] 0.9871109
```
