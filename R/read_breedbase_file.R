#' Read phenotype or design data from Excel, CSV, or tab-delimited files
#'
#' Loads phenotype or field-design data from `.xlsx`, `.xls`, `.csv`, `.tsv`,
#' `.txt`, or `.design` files. Excel files are read from all non-empty sheets
#' and combined into one tibble with a `.sheet` column identifying the source
#' sheet. CSV files are comma-delimited; `.tsv`, `.txt`, and `.design` files are
#' tab-delimited.
#'
#' @param file Path to the phenotype file.
#' @param sheet Optional Excel sheet name or number. If `NULL`, all non-empty
#'   Excel sheets are combined.
#' @param ... Additional arguments passed to `readxl::read_excel()` for Excel
#'   files, `utils::read.csv()` for CSV files, or `utils::read.table()` for
#'   tab-delimited files.
#'
#' @return A tibble containing the loaded phenotype data.
#'
#' @examples
#' \dontrun{
#' pheno <- read_breedbase_file("breedbase_trials_pheno.xlsx")
#' pheno_csv <- read_breedbase_file("breedbase_trials_pheno.csv")
#' design <- read_breedbase_file("augmented_row_column.design")
#' }
#'
#' @export
read_breedbase_file <- function(file, sheet = NULL, ...) {
  if (!is.character(file) || length(file) != 1 || is.na(file)) {
    stop("`file` must be a single file path.", call. = FALSE)
  }

  if (!file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(file))

  if (ext %in% c("xlsx", "xls")) {
    return(read_breedbase_excel_file(file, sheet = sheet, ...))
  }

  if (ext == "csv") {
    return(read_breedbase_csv_file(file, ...))
  }

  if (ext %in% c("tsv", "txt", "design")) {
    return(read_breedbase_tab_file(file, ...))
  }

  stop(
    "Unsupported file format: .",
    ext,
    ". Supported formats are .xlsx, .xls, .csv, .tsv, .txt, and .design.",
    call. = FALSE
  )
}

read_breedbase_excel_file <- function(file, sheet = NULL, ...) {
  if (!is.null(sheet)) {
    return(tibble::as_tibble(readxl::read_excel(file, sheet = sheet, ...)))
  }

  sheets <- readxl::excel_sheets(file)

  out <- purrr::map_dfr(
    sheets,
    function(sh) {
      df <- readxl::read_excel(file, sheet = sh, ...)

      if (nrow(df) == 0) {
        return(NULL)
      }

      dplyr::mutate(df, .sheet = sh)
    }
  )

  tibble::as_tibble(out)
}

read_breedbase_csv_file <- function(file, ...) {
  df <- utils::read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    ...
  )

  tibble::as_tibble(df)
}

read_breedbase_tab_file <- function(file, ...) {
  df <- utils::read.table(
    file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    ...
  )

  tibble::as_tibble(df)
}

load_or_use_data_frame <- function(x, label = "file") {
  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }

  if (is.character(x) && length(x) == 1 && !is.na(x)) {
    return(read_breedbase_file(x))
  }

  stop(
    "`",
    label,
    "` must be a file path or a data frame.",
    call. = FALSE
  )
}
