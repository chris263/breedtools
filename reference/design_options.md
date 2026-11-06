# Suggest trial layouts from the number of entries and available field space

Lists parameter options for the package's three design families,
starting with the information the breeder knows. Rows, columns and check
counts can be left unspecified. No randomization or optimization is
performed.

## Usage

``` r
design_options(
  n_entries,
  design = c("augmented_row_column", "prep", "resolvable"),
  n_checks = NULL,
  rows_in_field = NULL,
  cols_in_field = NULL,
  n_locations = 1L,
  n_reps = 2:3,
  block_size = c(5L, 8L, 10L, 12L, 16L, 20L),
  prep_prop = NULL,
  max_plots = NULL,
  max_dimension = NULL,
  min_error_df = NULL,
  n_options = 6L,
  n_candidates = 200L,
  rows = NULL,
  cols = NULL
)
```

## Arguments

- n_entries:

  Number of test entries, excluding checks.

- design:

  Design families to consider: `"augmented_row_column"` (augmented
  row-column), `"prep"` (partially replicated), and `"resolvable"`
  (fully replicated row-column). Defaults to all three.

- n_checks:

  Distinct check genotypes to consider, a non-negative integer or
  vector. `NULL` tries 2 through 6 for augmented row-column designs and
  zero for the other families. p-rep uses no separate checks; it is
  excluded when only positive check counts are requested. Resolvable
  designs replicate checks just like the other genotypes.

- rows_in_field, cols_in_field:

  Known dimensions of the whole field at one location (plot rows and
  plot columns, including all checks and replicates), not of a block.
  Block dimensions are suggested in the result as `rows_per_block` and
  `cols_per_block`. Each can be a positive integer, a vector of
  dimensions to consider, or `NULL` to search automatically. Supplied
  dimensions are exact candidate sizes, not upper bounds. When only one
  is supplied, the square field (the same number of rows and columns) is
  also considered. For p-rep, the last row may be partially occupied,
  but wholly empty final rows are not proposed.

- n_locations:

  Number of locations. Augmented row-column and resolvable layouts
  describe one field to randomize independently at each location. A
  p-rep option describes the combined multi-location design.

- n_reps:

  Full replication counts to try for resolvable designs (at least 2).
  Other design families use their own replication scheme.

- block_size:

  Maximum incomplete block sizes to try for p-rep.

- prep_prop:

  Proportions of entries to duplicate at each p-rep location. Counts use
  `round(n_entries * prep_prop)`, as in the generator. An entry may be
  duplicated at only one location, so the sum of these counts over
  locations must not exceed `n_entries`. `NULL` tries 15%, 20% and 25%,
  and, for several locations, a balanced share of entries per location
  rounded down to respect this budget.

- max_plots:

  Largest field area (rows times columns) per location. Defaults to four
  times the number of genotypes (increased when more full replicates are
  requested), or the largest supplied field area when both dimensions
  are given.

- max_dimension:

  Largest automatically searched row or column dimension. Defaults to
  twice the square root of `max_plots`, with a minimum of 20. Explicit
  dimensions are not restricted by this bound. When neither dimension is
  supplied, fields with an aspect ratio above 3 are excluded. Supply a
  dimension explicitly to consider longer fields.

- min_error_df:

  Minimum nominal error degrees of freedom. These are count-based values
  for the design model, assuming its effects are estimable, not degrees
  of freedom of a later spatial model. `NULL` (default) requires at
  least 1, except for the layout of Fig. 2 of Piepho and Williams
  (2016): blocks of `n_checks` x `n_checks` plots with each check
  exactly once, such as 36 entries and three checks in a 9 x 6 field.
  That layout is saturated with fixed rows, columns and blocks (error df
  = 0) and is analyzed with
  `analyze_augmented_row_column_design(random_blocks = TRUE)` or
  `entries = "random"`. A supplied value applies to every layout; 0 also
  admits other saturated layouts.

- n_options:

  Maximum options returned per design family. The shortlist includes
  different check, full-replication or p-rep counts before adding
  further layouts with the same counts.

- n_candidates:

  Search candidates to pass to the augmented row-column and p-rep
  generators. Defaults to 200. This is a user-selected search effort,
  not a dimension inferred from entry counts.

- rows, cols:

  Deprecated names of `rows_in_field` and `cols_in_field`.

## Value

A tibble with option id, design family, entry/check counts, field
dimensions, a `square_check_blocks` flag, location and replication
counts, active plots per location and in total, field area, unused
positions, check plots, repeated entries, blocks, nominal `error_df`,
method and notes. Field dimensions are named `rows_in_field` and
`cols_in_field`, matching
[`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md).
Block and replicate dimensions are included where applicable.
`n_candidates` reports the chosen search effort. `checks_per_block` is a
list-column with each augmented row-column option's check-allocation
matrix, or `NULL` for the default complete check blocks. These columns
support an explicit generator call. `rows_per_rowgroup` and
`cols_per_colgroup` are list-columns with the group sizes of
unequal-group options, or `NULL`. `generator` names the function to
call; `arguments` is a list-column of its layout arguments. Add your own
`treatments`, `controls` (augmented row-column) and `seed` before
calling it with [`do.call()`](https://rdrr.io/r/base/do.call.html).
p-rep arguments include placeholder location names, which can be
replaced with your own names. Attribute `diagnostics` reports the number
of shortlisted options and the reason for any excluded family. No
matches return a typed empty tibble and a message suggesting which
constraints to relax.

## Details

This is a parameter planner. Candidates pass plot-count, block-size and
nominal error-df checks. Augmented row-column candidates also pass the
check coverage and row/column group conditions used by the generator.
p-rep candidates have enough repeated plots to connect their blocks and
a pair-count budget compatible with (0,1) concurrence. These necessary
conditions do not certify a randomized design's connectivity or
efficiency: generate the chosen option, inspect its summary and run
[`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md).
The search is bounded and not exhaustive.

Augmented row-column options with equal row and column groups have
blocks of at least `n_checks` rows and at least `n_checks` columns,
leaving room for one copy of each check in a distinct row and column.
For example, four checks require blocks of at least 4 rows by 4 columns.
Every such block contains all checks at least once; extra positions are
assigned to extra check replication, so the number of test-entry
positions is exactly `n_entries`. Flexible allocations use
`method = "optimal"`.

Options with unequal row or column groups (Piepho and Williams, 2016,
Example 3 and Fig. 4) split the field into near-equal groups, such as
rows 6, 5, 6, 5 and columns 3, 4, 3 for a 22 x 10 field. Their blocks
hold `n_checks - 1` to `n_checks + 1` check plots, in proportion to
block size, so some blocks may be incomplete for the checks. These
options report `rows_per_rowgroup`, `cols_per_colgroup` and
`checks_per_block`, with `rows_per_block` and `cols_per_block` set to
`NA`.

Resolvable options may contain filler positions, which are reported in
`unused_plots`. Include both entries and checks in the generator's
`treatments`. Its native generator treats them all as treatments; if
checks need to be marked in Breedbase, set `is_control` in both the
generated `design` and `breedbase_design` tables before saving.

Square blocks are preferred but not required: for augmented row-column
designs, when a layout with blocks of `n_checks` rows by `n_checks`
columns containing each check exactly once exists (for example 3 x 3
blocks with checks A, B and C, as in Fig. 2 of Piepho and Williams,
2016), the smallest one is the first option and `square_check_blocks`
flags it. Blocks of that size with extra check plots are not flagged.
The other options follow in the usual order: within each resource-count
group, smaller field areas and fewer unused positions come first. Within
a field, augmented row-column layouts whose blocks hold about one check
plot per check (complete blocks for the checks, as in Piepho and
Williams, 2016) come first, then Latin constructions, squarer blocks and
less elongated fields. This ordering does not estimate statistical
precision or choose the biologically best design.

## See also

[`augmented_row_column_layouts()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_layouts.md),
[`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md),
[`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
[`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
[`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)

## Examples

``` r
# Only the number of entries is known:
options <- design_options(100)
options[, c("option", "design", "n_checks", "rows_in_field", "cols_in_field", "plots_per_location")]
#> # A tibble: 18 × 6
#>    option design         n_checks rows_in_field cols_in_field plots_per_location
#>     <int> <chr>             <int>         <int>         <int>              <dbl>
#>  1      1 augmented_row…        2             9            14                126
#>  2      2 augmented_row…        3             9            14                126
#>  3      3 augmented_row…        4            14             9                126
#>  4      4 augmented_row…        5             8            16                128
#>  5      5 augmented_row…        6            10            13                130
#>  6      6 augmented_row…        2            14             9                126
#>  7      7 prep                  0            13             9                115
#>  8      8 prep                  0            10            12                120
#>  9      9 prep                  0             9            14                125
#> 10     10 prep                  0            13             9                115
#> 11     11 prep                  0            12            10                120
#> 12     12 prep                  0             9            14                125
#> 13     13 resolvable            0            20            10                200
#> 14     14 resolvable            0            20            15                300
#> 15     15 resolvable            0            10            20                200
#> 16     16 resolvable            0            15            20                300
#> 17     17 resolvable            0            10            20                200
#> 18     18 resolvable            0            12            25                300

# A known field, four checks and exactly 108 test entries:
options <- design_options(108, design = "augmented_row_column", n_checks = 4,
                          rows_in_field = 12, cols_in_field = 12, n_candidates = 5)
selected <- options[1, ]
trial <- augmented_row_column_design(
  treatments = paste0("LINE_", seq_len(108)),
  controls = paste0("CHECK_", seq_len(selected$n_checks)),
  rows_in_field = selected$rows_in_field,
  cols_in_field = selected$cols_in_field,
  rows_per_block = selected$rows_per_block,
  cols_per_block = selected$cols_per_block,
  checks_per_block = selected$checks_per_block[[1]],
  method = selected$method,
  n_candidates = selected$n_candidates,
  seed = 2026
)

# Fully replicated options with one known dimension:
design_options(42, design = "resolvable", cols_in_field = 7, n_reps = 3)
#> # A tibble: 4 × 30
#>   option design     n_entries n_checks rows_in_field cols_in_field
#>    <int> <chr>          <int>    <int>         <int>         <int>
#> 1      1 resolvable        42        0            18             7
#> 2      2 resolvable        42        0            21             7
#> 3      3 resolvable        42        0            24             7
#> 4      4 resolvable        42        0            27             7
#> # ℹ 24 more variables: square_check_blocks <lgl>, n_locations <int>,
#> #   n_reps <int>, rows_per_block <int>, rows_per_rowgroup <list>,
#> #   cols_per_colgroup <list>, cols_per_block <int>, rows_per_rep <dbl>,
#> #   cols_per_rep <int>, block_size <int>, n_candidates <int>,
#> #   checks_per_block <list>, plots_per_location <dbl>, total_plots <dbl>,
#> #   field_plots <dbl>, unused_plots <dbl>, check_plots <int>, n_repeated <dbl>,
#> #   n_blocks <dbl>, error_df <dbl>, method <chr>, generator <chr>, …

# Only the number of field rows is known. With three checks, the first
# option is the layout of Piepho and Williams (2016, Fig. 2): a 9 x 6 field
# of six 3 x 3 blocks, each check once per block (error df = 0).
design_options(36, design = "augmented_row_column", n_checks = 3, rows_in_field = 9)
#> # A tibble: 6 × 30
#>   option design               n_entries n_checks rows_in_field cols_in_field
#>    <int> <chr>                    <int>    <int>         <int>         <int>
#> 1      1 augmented_row_column        36        3             9             6
#> 2      2 augmented_row_column        36        3             9             6
#> 3      3 augmented_row_column        36        3             9             6
#> 4      4 augmented_row_column        36        3             9             6
#> 5      5 augmented_row_column        36        3             9             7
#> 6      6 augmented_row_column        36        3             9             7
#> # ℹ 24 more variables: square_check_blocks <lgl>, n_locations <int>,
#> #   n_reps <int>, rows_per_block <int>, rows_per_rowgroup <list>,
#> #   cols_per_colgroup <list>, cols_per_block <int>, rows_per_rep <int>,
#> #   cols_per_rep <int>, block_size <int>, n_candidates <int>,
#> #   checks_per_block <list>, plots_per_location <int>, total_plots <dbl>,
#> #   field_plots <dbl>, unused_plots <dbl>, check_plots <dbl>, n_repeated <int>,
#> #   n_blocks <int>, error_df <dbl>, method <chr>, generator <chr>, …
# Require positive error df for every layout:
design_options(36, design = "augmented_row_column", n_checks = 3, rows_in_field = 9,
               min_error_df = 1)
#> # A tibble: 6 × 30
#>   option design               n_entries n_checks rows_in_field cols_in_field
#>    <int> <chr>                    <int>    <int>         <int>         <int>
#> 1      1 augmented_row_column        36        3             9             6
#> 2      2 augmented_row_column        36        3             9             6
#> 3      3 augmented_row_column        36        3             9             6
#> 4      4 augmented_row_column        36        3             9             7
#> 5      5 augmented_row_column        36        3             9             7
#> 6      6 augmented_row_column        36        3             9             7
#> # ℹ 24 more variables: square_check_blocks <lgl>, n_locations <int>,
#> #   n_reps <int>, rows_per_block <int>, rows_per_rowgroup <list>,
#> #   cols_per_colgroup <list>, cols_per_block <int>, rows_per_rep <int>,
#> #   cols_per_rep <int>, block_size <int>, n_candidates <int>,
#> #   checks_per_block <list>, plots_per_location <int>, total_plots <dbl>,
#> #   field_plots <dbl>, unused_plots <dbl>, check_plots <dbl>, n_repeated <int>,
#> #   n_blocks <int>, error_df <dbl>, method <chr>, generator <chr>, …
# Partial replication at three locations:
design_options(100, design = "prep", n_locations = 3)
#> # A tibble: 6 × 30
#>   option design n_entries n_checks rows_in_field cols_in_field
#>    <int> <chr>      <int>    <int>         <int>         <int>
#> 1      1 prep         100        0            13             9
#> 2      2 prep         100        0            10            12
#> 3      3 prep         100        0             9            14
#> 4      4 prep         100        0             7            19
#> 5      5 prep         100        0            13             9
#> 6      6 prep         100        0            12            10
#> # ℹ 24 more variables: square_check_blocks <lgl>, n_locations <int>,
#> #   n_reps <int>, rows_per_block <int>, rows_per_rowgroup <list>,
#> #   cols_per_colgroup <list>, cols_per_block <int>, rows_per_rep <int>,
#> #   cols_per_rep <int>, block_size <int>, n_candidates <int>,
#> #   checks_per_block <list>, plots_per_location <int>, total_plots <dbl>,
#> #   field_plots <dbl>, unused_plots <dbl>, check_plots <int>, n_repeated <int>,
#> #   n_blocks <dbl>, error_df <dbl>, method <chr>, generator <chr>, …
```
