#' Load location coordinates
#'
#' Loads location metadata from a `.xlsx`, `.xls`, or `.csv` file, or from an
#' already loaded data frame. The returned table has standardized columns:
#' `location`, `latitude`, and `longitude`.
#'
#' @param file Path to a location file, or an already loaded data frame.
#' @param location_col Optional location name column.
#' @param latitude_col Optional latitude column.
#' @param longitude_col Optional longitude column.
#' @param phenotype_data Optional raw phenotype data frame, prepared phenotype
#'   data frame, or result from `run_environment_stratification()`. If supplied,
#'   location names are extracted and checked against the location file.
#' @param phenotype_location_col Optional location column in `phenotype_data`.
#' @param match_by Column to extract when `phenotype_data` is a stratification
#'   result. Use `"location"`, `"environment"`, or `"environment_label"`.
#' @param reference_locations Optional character vector of expected location
#'   names. If supplied, every reference location must be present in `file`.
#' @param ignore_case Logical. If `TRUE`, location matching ignores case.
#' @param sheet Optional Excel sheet name or number passed to
#'   `read_breedbase_file()`.
#' @param ... Additional arguments passed to `read_breedbase_file()`.
#'
#' @return A tibble with `location`, `latitude`, and `longitude`.
#'
#' @examples
#' locations <- data.frame(
#'   location = c("Loc1", "Loc2"),
#'   latitude = c(35.1, 36.2),
#'   longitude = c(-78.9, -79.4)
#' )
#'
#' read_locations(locations)
#'
#' @export
read_locations <- function(
    file,
    location_col = NULL,
    latitude_col = NULL,
    longitude_col = NULL,
    phenotype_data = NULL,
    phenotype_location_col = NULL,
    match_by = "location",
    reference_locations = NULL,
    ignore_case = FALSE,
    sheet = NULL,
    ...
) {
  match_by <- match.arg(match_by, c("location", "environment", "environment_label"))

  loc_raw <- if (is.data.frame(file)) {
    tibble::as_tibble(file)
  } else if (is.character(file) && length(file) == 1 && !is.na(file)) {
    read_breedbase_file(file, sheet = sheet, ...)
  } else {
    stop("`file` must be a file path or a data frame.", call. = FALSE)
  }

  location_col <- location_col %||% guess_column(
    loc_raw,
    c("location", "location_name", "locationName", "name", "site", "environment")
  )
  latitude_col <- latitude_col %||% guess_column(
    loc_raw,
    c("latitude", "lat", "y", "decimalLatitude", "decimal_latitude")
  )
  longitude_col <- longitude_col %||% guess_column(
    loc_raw,
    c("longitude", "lon", "long", "lng", "x", "decimalLongitude", "decimal_longitude")
  )

  if (is.null(location_col)) {
    stop("Could not detect location name column.", call. = FALSE)
  }

  if (is.null(latitude_col)) {
    stop("Could not detect latitude column.", call. = FALSE)
  }

  if (is.null(longitude_col)) {
    stop("Could not detect longitude column.", call. = FALSE)
  }

  locations <- loc_raw |>
    dplyr::transmute(
      location = trimws(as.character(.data[[location_col]])),
      latitude = suppressWarnings(as.numeric(.data[[latitude_col]])),
      longitude = suppressWarnings(as.numeric(.data[[longitude_col]]))
    ) |>
    dplyr::filter(
      !is.na(.data$location),
      .data$location != "",
      !is.na(.data$latitude),
      !is.na(.data$longitude)
    ) |>
    dplyr::distinct(.data$location, .keep_all = TRUE)

  if (nrow(locations) == 0) {
    stop("No valid locations with latitude and longitude were found.", call. = FALSE)
  }

  if (!is.null(phenotype_data)) {
    reference_locations <- unique(c(
      reference_locations,
      extract_phenotype_locations(
        phenotype_data = phenotype_data,
        phenotype_location_col = phenotype_location_col,
        match_by = match_by
      )
    ))
  }

  if (!is.null(reference_locations)) {
    assert_locations_present(
      expected = reference_locations,
      available = locations$location,
      ignore_case = ignore_case,
      expected_label = "phenotype location"
    )
  }

  tibble::as_tibble(locations)
}

#' Create a map of environment groups
#'
#' Creates an interactive `leaflet` map where locations belonging to the same
#' environment group use the same marker color. Ungrouped environments are shown
#' with `ungrouped_color`.
#'
#' @param stratification_result Result from `run_environment_stratification()`
#'   or `run_breedbase_environment_stratification()`.
#' @param locations Location coordinates from `read_locations()`, a location
#'   file path, or a data frame with location, latitude, and longitude columns.
#' @param ignore_case Logical. If `TRUE`, location matching ignores case.
#' @param group_palette Optional vector of colors used for environment groups.
#' @param ungrouped_color Color for ungrouped environments.
#' @param match_by Name from the stratification result used to match the
#'   `location` column in `locations`. Use `"location"` when coordinates are
#'   named by location, `"environment"` when coordinates are named by the
#'   constructed environment identifier, or `"environment_label"` when
#'   coordinates are named by display labels.
#' @param strict Logical. If `TRUE`, stop when any environment location has no
#'   matching coordinates.
#'
#' @return A `leaflet` htmlwidget.
#'
#' @examples
#' \dontrun{
#' result <- run_environment_stratification("phenotype.csv", trait = "yield")
#' locations <- read_locations("locations.csv", ignore_case = TRUE)
#' map_environment_groups(result, locations, ignore_case = TRUE)
#' }
#'
#' @export
map_environment_groups <- function(
    stratification_result,
    locations,
    ignore_case = FALSE,
    group_palette = NULL,
    ungrouped_color = "#7A7A7A",
    match_by = "location",
    strict = TRUE
) {
  match_by <- match.arg(match_by, c("location", "environment", "environment_label"))

  map_data <- environment_group_location_data(stratification_result)

  if (nrow(map_data) == 0) {
    stop("No environment group data found to map.", call. = FALSE)
  }

  locs <- if (inherits(locations, "data.frame") &&
    all(c("location", "latitude", "longitude") %in% names(locations))) {
    tibble::as_tibble(locations)
  } else {
    read_locations(locations, ignore_case = ignore_case)
  }

  joined <- join_location_coordinates(
    map_data = map_data,
    locations = locs,
    ignore_case = ignore_case,
    match_by = match_by
  )

  missing_locations <- unique(joined[[match_by]][is.na(joined$latitude) | is.na(joined$longitude)])

  if (length(missing_locations) > 0 && isTRUE(strict)) {
    stop(
      "Missing coordinates for location(s): ",
      paste(missing_locations, collapse = ", "),
      ". Matching used `",
      match_by,
      "`. Available coordinate names include: ",
      paste(utils::head(unique(locs$location), 10), collapse = ", "),
      call. = FALSE
    )
  }

  joined <- joined |>
    dplyr::filter(!is.na(.data$latitude), !is.na(.data$longitude))

  if (nrow(joined) == 0) {
    stop("No mappable locations remain after matching coordinates.", call. = FALSE)
  }

  group_ids <- sort(unique(joined$group_id[joined$group_id != "Ungrouped"]))

  if (is.null(group_palette)) {
    group_palette <- stats::setNames(grDevices::hcl.colors(
      n = max(length(group_ids), 1),
      palette = "Dark 3"
    ), group_ids)
  } else {
    group_palette <- rep(group_palette, length.out = length(group_ids))
    group_palette <- stats::setNames(group_palette, group_ids)
  }

  joined$marker_color <- unname(group_palette[joined$group_id])
  joined$marker_color[joined$group_id == "Ungrouped" | is.na(joined$marker_color)] <- ungrouped_color

  leaflet::leaflet(joined) |>
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
    leaflet::addCircleMarkers(
      lng = ~longitude,
      lat = ~latitude,
      color = ~marker_color,
      fillColor = ~marker_color,
      fillOpacity = 0.85,
      radius = 7,
      stroke = TRUE,
      weight = 1,
      popup = ~paste0(
        "<strong>", htmltools::htmlEscape(environment_label), "</strong><br/>",
        "Location: ", htmltools::htmlEscape(location), "<br/>",
        "Group: ", htmltools::htmlEscape(group_id)
      ),
      label = ~paste(environment_label, group_id, sep = " - ")
    ) |>
    leaflet::addLegend(
      position = "bottomright",
      colors = c(unname(group_palette), ungrouped_color),
      labels = c(names(group_palette), "Ungrouped"),
      opacity = 0.85,
      title = "Environment group"
    )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

location_match_key <- function(x, ignore_case = FALSE) {
  x <- trimws(as.character(x))

  if (isTRUE(ignore_case)) {
    x <- tolower(x)
  }

  x
}

assert_locations_present <- function(
    expected,
    available,
    ignore_case = FALSE,
    expected_label = "location"
) {
  expected <- unique(location_match_key(expected, ignore_case = ignore_case))
  available <- unique(location_match_key(available, ignore_case = ignore_case))
  missing <- setdiff(expected, available)

  if (length(missing) > 0) {
    stop(
      "Missing coordinates for ",
      expected_label,
      "(s): ",
      paste(missing, collapse = ", "),
      ". Available coordinate names include: ",
      paste(utils::head(available, 10), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

extract_phenotype_locations <- function(
    phenotype_data,
    phenotype_location_col = NULL,
    match_by = "location"
) {
  if (is.list(phenotype_data) && !is.null(phenotype_data$env_info)) {
    env_info <- tibble::as_tibble(phenotype_data$env_info)

    if (match_by %in% names(env_info)) {
      return(clean_reference_locations(env_info[[match_by]]))
    }
  }

  if (is.list(phenotype_data) && !is.null(phenotype_data$pheno)) {
    phenotype_data <- phenotype_data$pheno
  }

  if (is.data.frame(phenotype_data)) {
    phenotype_data <- tibble::as_tibble(phenotype_data)

    if (!is.null(phenotype_location_col)) {
      if (!phenotype_location_col %in% names(phenotype_data)) {
        stop(
          "Phenotype location column not found: ",
          phenotype_location_col,
          call. = FALSE
        )
      }

      return(clean_reference_locations(phenotype_data[[phenotype_location_col]]))
    }

    if (match_by %in% names(phenotype_data)) {
      return(clean_reference_locations(phenotype_data[[match_by]]))
    }

    detected_col <- guess_column(
      phenotype_data,
      c(
        "location",
        "location_name",
        "locationName",
        "locationDbId",
        "location_db_id",
        "environment",
        "env",
        "trial_location",
        "studyLocation"
      )
    )

    if (!is.null(detected_col)) {
      return(clean_reference_locations(phenotype_data[[detected_col]]))
    }
  }

  stop(
    "Could not extract location names from `phenotype_data`. ",
    "Provide `phenotype_location_col` or `reference_locations`.",
    call. = FALSE
  )
}

clean_reference_locations <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

environment_group_location_data <- function(stratification_result) {
  if (!is.list(stratification_result) || is.null(stratification_result$group_membership)) {
    stop(
      "`stratification_result` must be a result from run_environment_stratification().",
      call. = FALSE
    )
  }

  grouped <- tibble::as_tibble(stratification_result$group_membership)

  if (nrow(grouped) > 0 && !"location" %in% names(grouped) &&
    !is.null(stratification_result$env_info)) {
    grouped <- add_environment_metadata(grouped, stratification_result$env_info)
  }

  ungrouped <- if (!is.null(stratification_result$ungrouped)) {
    tibble::as_tibble(stratification_result$ungrouped)
  } else {
    tibble::tibble()
  }

  if (nrow(ungrouped) > 0) {
    ungrouped$group_id <- "Ungrouped"

    if (!"location" %in% names(ungrouped) && !is.null(stratification_result$env_info)) {
      ungrouped <- add_environment_metadata(ungrouped, stratification_result$env_info)
    }
  }

  out <- dplyr::bind_rows(grouped, ungrouped)

  if (nrow(out) == 0) {
    return(tibble::tibble(
      group_id = character(),
      environment = character(),
      location = character(),
      environment_label = character()
    ))
  }

  if (!"location" %in% names(out)) {
    out$location <- out$environment
  }

  if (!"environment_label" %in% names(out)) {
    out$environment_label <- out$environment
  }

  out$group_id <- as.character(out$group_id)
  out$location <- as.character(out$location)
  out$environment_label <- as.character(out$environment_label)

  out |>
    dplyr::distinct(.data$group_id, .data$environment, .data$location, .data$environment_label)
}

join_location_coordinates <- function(
    map_data,
    locations,
    ignore_case = FALSE,
    match_by = "location"
) {
  map_data$.location_key <- location_match_key(map_data[[match_by]], ignore_case = ignore_case)

  locations <- locations |>
    dplyr::mutate(.location_key = location_match_key(.data$location, ignore_case = ignore_case)) |>
    dplyr::distinct(.data$.location_key, .keep_all = TRUE) |>
    dplyr::select(".location_key", "latitude", "longitude")

  map_data |>
    dplyr::left_join(locations, by = ".location_key") |>
    dplyr::select(-dplyr::all_of(".location_key"))
}
