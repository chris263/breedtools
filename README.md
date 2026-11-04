# breedtools

Design and analysis of augmented plant breeding trials.

- **Augmented row-column designs** (Piepho and Williams, 2016):
  unreplicated entries plus replicated checks in rows, columns and
  superimposed blocks. Row and column groups can be unequal, blocks can
  hold any number of check plots, and checks are placed by a D-optimal
  search (or covering Latin rectangles for regular layouts). Analyses can
  use random blocks, random entries (BLUPs), or a 'SpATS' spatial model.
- **Augmented p-rep designs** (Williams, Piepho and Whitaker, 2011):
  multi-location, partially replicated trials in incomplete blocks where
  pairs of entries share a block at most once. Choose the construction with
  `method`: a flexible randomized `"search"`, the paper's partial
  alpha-array construction (`"alpha"`), or its augmented alpha(2)-array
  construction (`"alpha2"`). `augmented_prep_parameters()` lists the
  feasible settings, and every design reports its efficiency factor.
- **Resolvable row-column designs** (Williams, John and Whitaker, 2006):
  fully replicated trials with rows and columns within replicates,
  optimized for row-column and spatial (linear variance) efficiency, with
  latinized options; analyses with RC, LV + RC and AR + RC models by REML.
- **Analysis**: ANOVA with degrees of freedom that are estimable from the
  observed plots, plus **adjusted entry means (BLUEs)** with standard
  errors, the mean SED and the 5% LSD. p-rep analyses fit replicates within
  locations and can add a fixed or random (REML) genotype-by-location
  interaction. A two-stage analysis (`analyze_augmented_prep_two_stage()`)
  handles heterogeneous error variances across locations and optional
  spatial adjustment with 'SpATS'.
- **Field maps and 'Breedbase'**: a plot-level field book and upload files,
  or save the design directly into a Breedbase breeding program and folder
  through BrAPI with `save_breedbase()`, after checking the location and
  accessions with `breedbase_check()`. After harvest, save the phenotypic
  data by plot and trait name with `save_breedbase_phenotypes()` (sheet from
  `phenotype_template()`), and the analysis results: first-stage plot values
  as a trait of the field plots (`save_breedbase_first_stage()`) and
  second-stage genotype means as a Breedbase analysis
  (`save_breedbase_analysis()`).
- **Multi-environment tools**: read 'Breedbase' phenotype files, run
  Hartley's Fmax variance check, and group environments by GxE.

## Installation

```r
# install.packages("remotes")
remotes::install_github("chris263/breedtools")
```

## Example

```r
library(breedtools)

rc <- augmented_row_column_design(
  treatments = sprintf("LINE_%03d", 1:108),
  controls = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"),
  rows_in_field = 12, cols_in_field = 12,
  rows_per_block = 4, cols_per_block = 4,
  seed = 2026
)
plot_augmented_row_column_design(rc)

field <- rc$design
field$yield <- rnorm(nrow(field), 60, 3)   # your phenotypes here
fit <- analyze_augmented_row_column_design(field, response = "yield")
fit                      # ANOVA and top adjusted means
fit$adjusted_means       # all entries ranked by adjusted mean
```

To see the full workflow, run `vignette("breedtools")`, or open the RStudio
script:

```r
file.edit(system.file("examples", "breedtools_workflow.R", package = "breedtools"))
```
