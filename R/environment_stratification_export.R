write_environment_stratification_json <- function(result, output_dir, prefix) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, paste0(prefix, ".json"))

  jsonlite::write_json(
    result,
    path = output_file,
    dataframe = "rows",
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )

  invisible(output_file)
}
