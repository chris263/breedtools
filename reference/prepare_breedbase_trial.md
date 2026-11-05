# Prepare an agricolae design and trial information for Breedbase

Converts an existing agricolae design and records the information needed
to save it as a new Breedbase trial. For example, generate an RCBD
(randomized complete block design) with
[`agricolae::design.rcbd()`](https://rdrr.io/pkg/agricolae/man/design.rcbd.html),
prepare it here, then call `save_breedbase(trial, con)`. Preparation is
offline; the randomized allocation, plot numbers and block assignments
are kept.

## Usage

``` r
prepare_breedbase_trial(
  design,
  trial_name = NULL,
  year = NULL,
  location = NULL,
  program = NULL,
  description = NULL,
  folder = NULL,
  checks = NULL,
  accession = NULL,
  design_type = NULL,
  trial_type = "phenotyping_trial",
  prompt = interactive(),
  ...
)
```

## Arguments

- design:

  An agricolae design, its field-book data frame, or a design already
  converted with
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md).

- trial_name:

  Name of the new trial.

- year:

  Year of the trial, as an integer or integer-valued string.

- location:

  Exact Breedbase location name, for one trial.

- program:

  Exact Breedbase breeding program name.

- description:

  Description of the trial's purpose, material and season.

- folder:

  Optional Breedbase folder name. `NULL` saves directly in the breeding
  program. Use `create_folder = TRUE` in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  if a new folder should be created.

- checks, accession, design_type:

  As in
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md).
  Supply these when converting a raw agricolae design or field book.

- trial_type:

  Breedbase trial type; defaults to `"phenotyping_trial"`.

- prompt:

  Logical. Ask for missing required trial information in the R console.
  Defaults to `TRUE` in an interactive session and `FALSE` in scripts.
  With `FALSE`, all missing fields are reported in one error.

- ...:

  Additional column mappings for
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md),
  such as `plot`, `block`, `rep`, `row` and `col`.

## Value

A `breedbase_trial` object extending `breedbase_design`, with converted
`design`, original `book` and `parameters`, `design_type` and the
recorded `trial` information. Pass the whole object to
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
or
[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md).

## Details

Required trial information is `trial_name`, `year`, `location`,
`program` and `description`. Provide it as named arguments for a
reproducible script. In RStudio, `prepare_breedbase_trial(design)` asks
for these fields when they are missing.

Inspect `trial$trial` for the trial information and `trial$design` for
the converted plot table. `phenotype_template(trial)` uses the recorded
trial name to make a field sheet with the plot names that will be saved.
`save_breedbase(trial, dry_run = TRUE)` previews the trial and plot
payload without a connection. With a connection, it checks the breeding
program, location, accessions and uniqueness of the new trial and plot
names. `save_breedbase(trial, con)` then saves through the existing
BrAPI workflow. Explicit arguments to the save function override
recorded information, for example `trial_name = "Another_trial"`.

No database names or accessions are created during preparation. Physical
row and column coordinates are optional for saving; include them in the
input field book when known, especially for later spatial analyses. They
are not inferred from agricolae plot codes.

## See also

[`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md),
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md),
[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)

## Examples

``` r
if (requireNamespace("agricolae", quietly = TRUE)) {
  entries <- c("G1", "G2", "G3", "CHECK_A")
  rcbd <- agricolae::design.rcbd(entries, r = 3, seed = 2026)
  trial <- prepare_breedbase_trial(rcbd,
    trial_name = "RCBD_2026_Ibadan", year = 2026, location = "Ibadan",
    program = "IITA", description = "RCBD with three replicates",
    checks = "CHECK_A", prompt = FALSE)
  trial$trial
  head(trial$design)
  payload <- save_breedbase(trial, dry_run = TRUE)
  sheet <- phenotype_template(trial, traits = "fresh root yield")
}
#> Dry run without a connection: returning the data that would be sent.
if (FALSE) { # \dontrun{
con <- breedbase_login("https://your-breedbase.org", "your_username")
save_breedbase(trial, con, dry_run = TRUE)  # check names in the database
save_breedbase(trial, con)                  # save the new trial
} # }
```
