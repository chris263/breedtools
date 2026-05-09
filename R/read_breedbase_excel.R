#' Read a Breedbase Excel phenotype file
#'
#' `read_breedbase_excel()` is deprecated. Use `read_breedbase_file()` for
#' `.xlsx`, `.xls`, and `.csv` files.
#'
#' @param file Path to the Excel file.
#' @param ... Additional arguments passed to `read_breedbase_file()`.
#'
#' @return A tibble containing combined data from all non-empty sheets.
#'
#' @examples
#' \dontrun{
#' pheno_raw <- read_breedbase_file("breedbase_trials_pheno.xlsx")
#' }
#'
read_breedbase_excel <- function(file, ...) {
  .Deprecated("read_breedbase_file")
  read_breedbase_file(file, ...)
}
