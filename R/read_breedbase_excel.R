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
#' xlsx_file <- system.file("extdata", "breedbase_pheno.xlsx", package = "breedtools")
#'
#' # Deprecated (emits a deprecation warning): use read_breedbase_file()
#' pheno <- suppressWarnings(read_breedbase_excel(xlsx_file))
#'
#' # Preferred
#' pheno <- read_breedbase_file(xlsx_file)
#' head(pheno)
#'
#' @keywords internal
#' @export
read_breedbase_excel <- function(file, ...) {
  .Deprecated("read_breedbase_file")
  read_breedbase_file(file, ...)
}
