#' Prepare an agricolae design and trial information for Breedbase
#'
#' Converts an existing agricolae design and records the information needed
#' to save it as a new Breedbase trial. For example, generate an RCBD
#' (randomized complete block design) with `agricolae::design.rcbd()`, prepare
#' it here, then call `save_breedbase(trial, con)`. Preparation is offline;
#' the randomized allocation, plot numbers and block assignments are kept.
#'
#' @param design An agricolae design, its field-book data frame, or a design
#'   already converted with [as_breedbase_design()].
#' @param trial_name Name of the new trial.
#' @param year Year of the trial, as an integer or integer-valued string.
#' @param location Exact Breedbase location name, for one trial.
#' @param program Exact Breedbase breeding program name.
#' @param description Description of the trial's purpose, material and season.
#' @param folder Optional Breedbase folder name. `NULL` saves directly in
#'   the breeding program. Use `create_folder = TRUE` in [save_breedbase()]
#'   if a new folder should be created.
#' @param checks,accession,design_type As in [as_breedbase_design()]. Supply
#'   these when converting a raw agricolae design or field book.
#' @param trial_type Breedbase trial type; defaults to `"phenotyping_trial"`.
#' @param prompt Logical. Ask for missing required trial information in the
#'   R console. Defaults to `TRUE` in an interactive session and `FALSE` in
#'   scripts. With `FALSE`, all missing fields are reported in one error.
#' @param ... Additional column mappings for [as_breedbase_design()], such
#'   as `plot`, `block`, `rep`, `row` and `col`.
#'
#' @details
#' Required trial information is `trial_name`, `year`, `location`,
#' `program` and `description`. Provide it as named arguments for a
#' reproducible script. In RStudio, `prepare_breedbase_trial(design)` asks
#' for these fields when they are missing.
#'
#' Inspect `trial$trial` for the trial information and `trial$design` for
#' the converted plot table. `phenotype_template(trial)` uses the recorded
#' trial name to make a field sheet with the plot names that will be saved.
#' `save_breedbase(trial, dry_run = TRUE)` previews the trial and plot
#' payload without a connection. With a connection, it checks the breeding
#' program, location, accessions and uniqueness of the new trial and plot
#' names. `save_breedbase(trial, con)` then saves through the existing BrAPI
#' workflow. Explicit arguments to the save function override recorded
#' information, for example `trial_name = "Another_trial"`.
#'
#' No database names or accessions are created during preparation.
#' Physical row and column coordinates are optional for saving; include
#' them in the input field book when known, especially for later spatial
#' analyses. They are not inferred from agricolae plot codes.
#'
#' @return A `breedbase_trial` object extending `breedbase_design`, with
#'   converted `design`, original `book` and `parameters`, `design_type`
#'   and the recorded `trial` information. Pass the whole object to
#'   [save_breedbase()] or [phenotype_template()].
#' @seealso [as_breedbase_design()], [save_breedbase()], [phenotype_template()]
#' @examples
#' if (requireNamespace("agricolae", quietly = TRUE)) {
#'   entries <- c("G1", "G2", "G3", "CHECK_A")
#'   rcbd <- agricolae::design.rcbd(entries, r = 3, seed = 2026)
#'   trial <- prepare_breedbase_trial(rcbd,
#'     trial_name = "RCBD_2026_Ibadan", year = 2026, location = "Ibadan",
#'     program = "IITA", description = "RCBD with three replicates",
#'     checks = "CHECK_A", prompt = FALSE)
#'   trial$trial
#'   head(trial$design)
#'   payload <- save_breedbase(trial, dry_run = TRUE)
#'   sheet <- phenotype_template(trial, traits = "fresh root yield")
#' }
#' \dontrun{
#' con <- breedbase_login("https://your-breedbase.org", "your_username")
#' save_breedbase(trial, con, dry_run = TRUE)  # check names in the database
#' save_breedbase(trial, con)                  # save the new trial
#' }
#' @export
prepare_breedbase_trial <- function(design, trial_name = NULL, year = NULL, location = NULL,
                                    program = NULL, description = NULL, folder = NULL,
                                    checks = NULL, accession = NULL, design_type = NULL,
                                    trial_type = "phenotyping_trial", prompt = interactive(), ...) {
  if (!is.logical(prompt) || length(prompt) != 1L || is.na(prompt)) {
    stop("`prompt` must be TRUE or FALSE.", call. = FALSE)
  }
  fields <- list(trial_name = trial_name, year = year, location = location,
                 program = program, description = description)
  is_missing <- function(x) {
    is.null(x) || length(x) == 0L ||
      (length(x) == 1L && (is.na(x) || !nzchar(trimws(as.character(x)))))
  }
  missing_fields <- names(fields)[vapply(fields, is_missing, logical(1))]
  if (isTRUE(prompt)) {
    labels <- c(trial_name = "New trial name", year = "Trial year", location = "Breedbase location name",
                program = "Breedbase breeding program", description = "Trial description")
    for (field in missing_fields) fields[[field]] <- bb_prompt_trial_field(labels[[field]])
    missing_fields <- names(fields)[vapply(fields, is_missing, logical(1))]
  }
  if (length(missing_fields)) {
    stop("Give the required trial information: ", paste(missing_fields, collapse = ", "),
         ". Supply named arguments or use prompt = TRUE in an interactive session.", call. = FALSE)
  }
  text_field <- function(value, name) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
      stop("`", name, "` must be one non-empty character string.", call. = FALSE)
    }
    trimws(value)
  }
  for (field in setdiff(names(fields), "year")) fields[[field]] <- text_field(fields[[field]], field)
  y <- fields$year
  if ((!is.numeric(y) && !is.character(y)) || length(y) != 1L) {
    stop("`year` must be one positive integer year.", call. = FALSE)
  }
  y <- suppressWarnings(as.numeric(y))
  if (!is.finite(y) || y < 1 || y > .Machine$integer.max || y != floor(y)) {
    stop("`year` must be one positive integer year.", call. = FALSE)
  }
  fields$year <- as.integer(y)
  fields$trial_name <- bb_clean_name(fields$trial_name, "trial name")
  fields["folder"] <- list(if (is.null(folder)) NULL else text_field(folder, "folder"))
  fields$trial_type <- text_field(trial_type, "trial_type")
  conversion_args <- list(...)
  if (inherits(design, "breedbase_design")) {
    if (!is.null(checks) || !is.null(accession) || !is.null(design_type) || length(conversion_args)) {
      stop("This design is already converted. Supply check and column mappings when converting the original field book.",
           call. = FALSE)
    }
    converted <- design
  } else {
    converted <- as_breedbase_design(design, accession = accession, checks = checks,
                                     design_type = design_type, ...)
  }
  converted$trial <- fields
  class(converted) <- unique(c("breedbase_trial", class(converted)))
  converted
}

bb_prompt_trial_field <- function(label) {
  trimws(readline(paste0(label, ": ")))
}
