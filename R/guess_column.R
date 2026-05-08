#' Guess a column name from possible Breedbase column names
#'
#' @param df A data frame.
#' @param possible_names A character vector of possible column names.
#'
#' @return A column name if found, otherwise `NULL`.
#'
#' @examples
#' df <- data.frame(germplasmName = "G1", phenotype = 10)
#' guess_column(df, c("germplasmName", "accession_name"))
#'
#' @export
guess_column <- function(df, possible_names) {
  nms <- names(df)
  nms_lower <- tolower(nms)

  possible_lower <- tolower(possible_names)

  hit <- which(nms_lower %in% possible_lower)

  if (length(hit) > 0) {
    return(nms[hit[1]])
  }

  NULL
}
