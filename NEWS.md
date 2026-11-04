# breedtools 0.1.0

* First CRAN release.
* New `breedbase_check()` checks the location(s) and every entry and check
  (as accessions) in Breedbase before a design is built, plus the program,
  folder and trial name when given. All problems are listed together with
  close matches (e.g. "ibadan" -> "Ibadan"). `save_breedbase()` now reports
  location and accession problems together too.
* Breedbase searches work on slow servers: accessions are searched 25 at a
  time and a batch that times out (HTTP 504/502/503/429) is split and tried
  again, with a line after each batch (names done, percent, time). Accessions found
  are remembered for the R session, so `save_breedbase()` after
  `breedbase_check()` does not search them again.
* `save_breedbase()` survives server timeouts when writing: after a timeout
  it checks whether the folder, trial or plots were created before trying
  again (no duplicates), and a re-run continues a save that was cut off
  (an existing folder is used; a trial without plots gets its plots).
* `breedbase_check()` says plainly when the trial already exists ("choose
  another trial_name"), and that an existing folder is fine: the trial is
  saved inside it. `save_breedbase()` also uses a folder that Breedbase
  reports as existing when it tries to create it.
* The `folder` of `save_breedbase()` is optional: without it the trial is
  saved directly in the breeding program.
* `save_breedbase()` sends the `germplasmDbId` of each accession with the
  plots, which Breedbase requires ("Verify Germplasm Id." error before).
* New `save_breedbase_phenotypes()` saves phenotypic data to Breedbase by
  plot name and trait name (wide: one column per trait, or long: trait and
  value columns; optional time stamps for repeated measurements). Plots and
  traits are checked first, values are sent in batches (surviving server
  timeouts) and read back to verify. New `phenotype_template()` makes the
  phenotyping sheet with the Breedbase plot names of a design.
* Analysis results to Breedbase. First stage: `first_stage_values()` gives
  one value per field plot (adjusted mean of the genotype + plot residual,
  design and spatial effects removed) and `save_breedbase_first_stage()`
  saves it on the field plots as its own trait: the measured trait composed
  with a statistic ("fresh root yield|Adjusted Means from Spatial Correction
  using SpATS R"), found or created by `breedbase_adjusted_trait()`. Second
  stage: `second_stage_values()` collects the genotype adjusted means and
  `save_breedbase_analysis()` saves them as a Breedbase analysis, like the
  Breedbase mixed-models tool.
* `analyze_resolvable_row_column()` now returns `plot_results` with the
  fitted values and plot residuals. See `options(breedtools.germplasm_batch)` and
  `options(breedtools.breedbase_timeout)`.
* New `save_breedbase()` saves any design (augmented row-column, p-rep,
  resolvable row-column) as a trial with its plots in a 'Breedbase' database
  through BrAPI v2, in a chosen breeding program and folder. It checks the
  program, folder (optionally created), locations, accessions and trial name
  before writing, asks for confirmation, and reads the plots back to verify.
  `dry_run = TRUE` makes all checks without writing. Plot names are made
  automatically as `trialName_accession_block_rep` (`plot_name_pattern`) or
  given by you (`plot_names`), and are checked to be new in the database. `breedbase_login()`
  logs in without storing the password; `breedbase_programs()`,
  `breedbase_folders()` and `breedbase_locations()` list what exists.
* New resolvable row-column designs (Williams, John and Whitaker, 2006):
  * `resolvable_row_column_design()` builds designs for `v <= k * s`
    treatments in `r` replicates of `k x s`, optimizing the row-column, row
    and column efficiency factors and then the spatial (LV + RC) efficiency
    factor. Latinized and t-latinized designs for contiguous replicates,
    filler plots, and improving or laying out an existing design
    (`start_design`) are supported.
  * `resolvable_row_column_efficiency()` reports the efficiency factors, the
    spatial upper bound, the AR-model criterion z(X), and treatment
    neighbour and edge counts.
  * `analyze_resolvable_row_column()` fits the RC, LV + RC and AR + RC
    models by REML and compares them by deviance; long columns of latinized
    designs can be fitted as fixed.
  * `plot_resolvable_row_column_design()` draws the field.
* `augmented_row_column_design()`:
  * Check placement now uses covering Latin rectangles, so every field row
    and column gets a check.
  * Layouts whose effects cannot be estimated, or that leave no error
    degrees of freedom, are rejected.
  * The design is much faster to generate.
* New helpers for checking designs and simulations: `check_setup()` checks
  at the start of a script that breedtools, the optional 'lme4' and 'SpATS'
  packages and the output folder are ready, `check_design()` runs
  the standard checks of any design and prints `[OK]` / `[FAIL]` lines,
  `check_result()` records a check of your own, `check_section()` groups
  the following checks, `check_summary()` reports what was checked by section
  and what is correct, `check_reset()` clears the tally (saying how many
  checks it cleared), and `prediction_accuracy()`
  gives the correlation between adjusted means and known true values.
* New `augmented_row_column_layouts()` lists every analyzable layout of a
  field (block sizes and numbers of checks), most entry plots first, with
  optional `n_entries` and `min_error_df` filters.
* `augmented_row_column_capacity()` reports the error degrees of freedom of
  a layout, the share of check plots, and the row-group and column-group
  conditions of Piepho and Williams (2016).
* Augmented row-column designs now support the full strategy of Piepho and
  Williams (2016): row and column groups of unequal size
  (`rows_per_rowgroup`, `cols_per_colgroup`) and any number of check plots
  per block (`checks_per_block`), giving incomplete blocks or extra check
  replication. The new `method = "optimal"` allocates checks by a
  D-optimal local search that keeps check plots from being neighbours
  (same column in consecutive rows, or same row in consecutive columns) and
  keeps the check allocation balanced: a block with at least as many check
  plots as checks holds every check (extra plots repeat checks), an
  incomplete block holds each check at most once, rows and columns of a group
  get even numbers of check plots, and a column (or row) with room for all
  checks holds every check. Figures 3 and 4 of Piepho and Williams (2016) are
  reproduced with all their stated properties; `"latin"` keeps the covering Latin rectangle
  construction (now with evenly spaced check rows and columns, and polished
  by swapping rows or columns within groups to remove neighbouring checks), and `"auto"` picks it when it applies. `summary` reports
  the D-criterion.
* `analyze_augmented_row_column_design()` can take rows, columns and blocks
  as random to recover inter-block information (`random_blocks = TRUE`),
  model entries as random with fixed checks (`entries = "random"`, BLUPs),
  or add a spatial P-spline surface with 'SpATS' (`spatial = TRUE`).
* `analyze_augmented_row_column_design()`:
  * Degrees of freedom are corrected (rows within row groups `k - gk`,
    columns within column groups `s - gs`), and the ANOVA uses the df that
    are estimable from the observed plots.
  * The model is fitted by absorbing unreplicated entries, so it is fast for
    large fields.
* `augmented_prep_design()` uses a faster candidate search with incremental
  concurrence scoring. The search stops early once an optimal design is
  found.
* `augmented_prep_design()` follows Williams, Piepho and Whitaker (2011)
  more closely:
  * p-rep plots are spread evenly over the blocks of each location;
  * the blocks of each location are kept connected by the p-rep entries, so
    each location can be analysed on its own (a warning is given when there
    are too few p-rep plots to do so);
  * a swap-based improvement step removes pairs of entries that share more
    than one block, giving (0,1)-designs for cases such as the paper's
    six-location trial with blocks of eight; repeated concurrences are
    penalized quadratically, so when a (0,1)-design is out of reach the
    excess is spread over pairs instead of one pair meeting at every
    location, and greedy candidates that put an entry twice in a block are
    repaired instead of discarded;
  * new `resolvable = TRUE` option builds a resolvable design at each
    location (every p-rep entry once per replicate), returned in the new
    `replicate` column;
  * `summary` reports the average efficiency factor (`efficiency_factor`)
    and the number of disconnected locations; among optimal candidates the
    most efficient one is returned (`n_efficiency_candidates`).
  * new `method` argument chooses the construction: `"search"` (default),
    `"alpha"` (partial alpha-arrays, Sections 2.1-2.2 of the paper) or
    `"alpha2"` (augmented alpha(2)-arrays, Section 2.4), with
    `n_prep_per_block` for the number of p-rep entries per block (`m`).
* New `analyze_augmented_prep_two_stage()`: two-stage analysis of
  multi-location p-rep trials. Stage 1 analyses each location with its own
  error variance, either with the incomplete-block model or with a spatial
  P-spline model from the suggested 'SpATS' package (`stage1 = "spatial"`).
  Stage 2 combines the location means with full or diagonal stage-1
  weighting (Möhring and Piepho, 2009) and a REML genotype-by-location
  variance. It reports Hartley's Fmax test of heterogeneous error variances.
* `analyze_augmented_prep_design()` returns the covariance matrix of the
  adjusted means (`$adjusted_means_vcov`).
* Methods `"alpha"` and `"alpha2"` stop searching once further random starts
  no longer help, and warn when no (0,1)-design was found.
* New `augmented_prep_parameters()` lists the feasible parameter sets of the
  algebraic constructions (as in Table 3 of the paper).
* `map_environment_groups()` now uses the same maps as the 'Breedbase'
  Locations page by default (`tiles = "breedbase"`): Esri Street,
  Topographical and Satellite, with a switch between them. They need no
  account or API key and work in the RStudio viewer. `tiles` can also be
  `"street"`, `"topo"`, `"satellite"`, `"OpenStreetMap"` or a leaflet
  provider name.
* `plot_augmented_prep_design()` draws clearer maps: entries repeated at the
  location in orange with bold names (they replace the checks), single plots
  in light grey, thick black lines around blocks, a blue line and labels for
  replicates, a plain-language legend, a caption explaining the map
  (`caption = TRUE`), and locations in the order of the design. It can also
  colour plots by replicate.
* `analyze_augmented_prep_design()`:
  * Handles single-location trials.
  * Reports df that are estimable when plots are missing.
* `analyze_augmented_prep_design()` fits replicates within locations for
  resolvable designs (`replicate_col`) and can add a genotype-by-location
  interaction, fixed (`genotype_by_location = "fixed"`) or random
  (`"random"`, REML via the suggested 'lme4' package, with variance
  components in `$variance_components`).
* Both analyses now return adjusted entry means (BLUEs) with standard errors
  (`$adjusted_means`), and report the mean SED and the 5% LSD in
  `$diagnostics`.
* New tutorial, `vignette("breedtools")`, covering every design and analysis
  and the multi-environment tools, and an RStudio workflow script in
  `inst/examples`. The source package includes `scripts/validate_designs.R`,
  a step-by-step tour of every design and analysis with simple checks, and
  `scripts/run_all_checks.R`, the full automated validation against the
  published examples, `lm()`, 'lme4' and simulation.
