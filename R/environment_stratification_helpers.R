#' Normalize Breedbase-style column names
#'
#' Removes ontology suffixes such as `.CO...` or `|CO_...` and replaces dots
#' with spaces.
#'
#' @param x Character vector of column names.
#'
#' @return Normalized character vector.
#'
#' @examples
#' normalize_column_name(c(
#'   "fresh root yield|CO_334:0000013",
#'   "plant.height.CO_334.0000018",
#'   "germplasmName"
#' ))
#'
#' @keywords internal
#' @export
normalize_column_name <- function(x) {
  x <- gsub("\\.CO.*", "", x)
  x <- gsub("\\|CO_.*", "", x)
  x <- gsub("\\.", " ", x)
  x
}

#' Find the first matching column
#'
#' @param data A data frame.
#' @param candidates Candidate column names.
#' @param label Label used in the error message.
#'
#' @return The first matching column name.
#'
#' @examples
#' df <- data.frame(locationName = "Loc1", germplasmName = "G1")
#' find_column(df, c("germplasmName", "accession_name"), "accession")
#'
#' @keywords internal
#' @export
find_column <- function(data, candidates, label = "column") {
  matches <- candidates[candidates %in% colnames(data)]

  if (length(matches) == 0) {
    stop(
      "Could not find required column for ",
      label,
      ". Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  matches[1]
}

#' Calculate harmonic mean
#'
#' @param x Numeric vector.
#'
#' @return Harmonic mean.
#'
#' @examples
#' harmonic_mean(c(2, 3, 4, NA))
#'
#' @keywords internal
#' @export
harmonic_mean <- function(x) {
  x <- x[!is.na(x) & x > 0]

  if (length(x) == 0) {
    return(NA_real_)
  }

  length(x) / sum(1 / x)
}

clean_display_values <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- x[x != ""]
  unique(x)
}

environment_display_label <- function(environment, location, trial, year) {
  parts <- clean_display_values(c(location, trial, year))

  if (length(parts) == 0) {
    return(as.character(environment))
  }

  paste(parts, collapse = " / ")
}

#' Complete environment metadata
#'
#' @param env_info Data frame with environment, location, trial, and year.
#'
#' @return Data frame with an added `environment_label` column.
#'
#' @examples
#' env_info <- data.frame(
#'   environment = c("Loc1_T1", "Loc2_T2"),
#'   location = c("Loc1", "Loc2"),
#'   trial = c("T1", "T2"),
#'   year = c("2025", NA)
#' )
#'
#' complete_environment_info(env_info)
#'
#' @export
complete_environment_info <- function(env_info) {
  if (nrow(env_info) == 0) {
    env_info$environment_label <- character()
    return(env_info)
  }

  env_info$location <- as.character(env_info$location)
  env_info$trial <- as.character(env_info$trial)
  env_info$year <- as.character(env_info$year)

  env_info$location[is.na(env_info$location)] <- ""
  env_info$trial[is.na(env_info$trial)] <- ""
  env_info$year[is.na(env_info$year)] <- ""

  env_info$environment_label <- mapply(
    environment_display_label,
    env_info$environment,
    env_info$location,
    env_info$trial,
    env_info$year,
    USE.NAMES = FALSE
  )

  env_info
}

#' Add environment metadata to result table
#'
#' @param results Result data frame.
#' @param env_info Environment metadata table.
#' @param env_col Name of the environment column in `results`.
#'
#' @return Joined data frame.
#'
#' @examples
#' pheno <- read_breedbase_file(
#'   system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
#' )
#' prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")
#'
#' membership <- data.frame(
#'   group_id = "Group_1",
#'   environment = prepared$env_info$environment[1:2]
#' )
#' add_environment_metadata(membership, prepared$env_info)
#'
#' @export
add_environment_metadata <- function(results, env_info, env_col = "environment") {
  if (nrow(results) == 0 || !(env_col %in% colnames(results))) {
    return(results)
  }

  dplyr::left_join(
    results,
    env_info,
    by = stats::setNames("environment", env_col)
  )
}

#' Add environment metadata to pairwise Lin results
#'
#' @param pairwise Pairwise result table.
#' @param env_info Environment metadata table.
#'
#' @return Pairwise table with metadata.
#'
#' @examples
#' pheno <- read_breedbase_file(
#'   system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
#' )
#' prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")
#'
#' lin <- lin_group_environments(prepared$data)
#' add_pairwise_environment_metadata(lin$pairwise, prepared$env_info)
#'
#' @export
add_pairwise_environment_metadata <- function(pairwise, env_info) {
  if (nrow(pairwise) == 0) {
    return(pairwise)
  }

  env1_info <- env_info |>
    dplyr::rename(
      env1 = "environment",
      env1_name = "environment_label",
      env1_location = "location",
      env1_trial = "trial",
      env1_year = "year"
    )

  env2_info <- env_info |>
    dplyr::rename(
      env2 = "environment",
      env2_name = "environment_label",
      env2_location = "location",
      env2_trial = "trial",
      env2_year = "year"
    )

  pairwise |>
    dplyr::left_join(env1_info, by = "env1") |>
    dplyr::left_join(env2_info, by = "env2") |>
    dplyr::select(
      "env1_location",
      "env1_trial",
      "env1_year",
      "env1_name",
      "env2_location",
      "env2_trial",
      "env2_year",
      "env2_name",
      dplyr::everything(),
      -dplyr::any_of(c("env1", "env2", "environments"))
    )
}

#' Summarize environment names by group
#'
#' @param group_summary Group summary table.
#' @param group_membership Group membership table.
#'
#' @return Group summary with display metadata.
#'
#' @examples
#' pheno <- read_breedbase_file(
#'   system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
#' )
#' prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")
#'
#' lin <- lin_group_environments(prepared$data)
#' membership <- add_environment_metadata(lin$group_membership, prepared$env_info)
#' environment_summary_by_group(lin$group_summary, membership)
#'
#' @export
environment_summary_by_group <- function(group_summary, group_membership) {
  if (nrow(group_summary) == 0 || nrow(group_membership) == 0) {
    return(group_summary |> dplyr::select(-dplyr::any_of("environments")))
  }

  display_summary <- group_membership |>
    dplyr::group_by(.data$group_id) |>
    dplyr::summarise(
      environments = paste(clean_display_values(.data$environment_label), collapse = ", "),
      locations = paste(clean_display_values(.data$location), collapse = ", "),
      trials = paste(clean_display_values(.data$trial), collapse = ", "),
      years = paste(clean_display_values(.data$year), collapse = ", "),
      .groups = "drop"
    )

  group_summary |>
    dplyr::select(-dplyr::any_of("environments")) |>
    dplyr::left_join(display_summary, by = "group_id") |>
    dplyr::select(
      "group_id",
      "environments",
      "locations",
      "trials",
      "years",
      dplyr::everything()
    )
}
