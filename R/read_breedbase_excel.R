#' Read a Breedbase Excel phenotype file
#'
#' Reads all sheets from a Breedbase Excel file and combines them into one
#' data frame. A `.sheet` column is added to identify the source sheet.
#'
#' @param file Path to the Excel file.
#'
#' @return A tibble containing combined data from all non-empty sheets.
#'
#' @examples
#' \dontrun{
#' pheno_raw <- read_breedbase_excel("breedbase_trials_pheno.xlsx")
#' }
#'
#' @export
read_breedbase_excel <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }

  sheets <- readxl::excel_sheets(file)

  out <- purrr::map_dfr(
    sheets,
    function(sh) {
      df <- readxl::read_excel(file, sheet = sh)

      if (nrow(df) == 0) {
        return(NULL)
      }

      dplyr::mutate(df, .sheet = sh)
    }
  )

  tibble::as_tibble(out)
}
