#' Function takes csv file and converts is to saved object
#'
#' @param filename A file location (csv file) with births and deaths of each animal
#'
#' @export
#' @importFrom usethis use_data
#' @examples
#' \dontrun{
#' writedata("inst/extdata/whaleData_08-06-2018.csv")
#' }
writedata <- function(filename) {
  orca = read.csv(filename, stringsAsFactors = FALSE)
  use_data(orca, overwrite=TRUE)
}
