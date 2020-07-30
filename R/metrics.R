## Integer approaches -----

#' Metric for determining the most popular value
#'
#' Computes the mode given a vector of data.
#'
#' @param x      A `vector` of data.
#' @param na.rm  A `logical` indicating if missing values (including NaN) should
#'               be removed. Default: `FALSE`
#'
#' @return
#' A single `numeric` or `integer` value.
#'
#' @section Recovery Use:
#' The mode should be used when computing latent class recovery among subjects.
#'
#' @details
#' Be forewarned that this method of obtaining the mode does not take into
#' consideration ties. That is, only the first group with the maximum value
#' is returned. If a second group also has a similar count, this group will
#' not be returned.
#'
#' @export
#' @examples
#' # Sample class-data
#' x = c(0, 0, 0, 1, 1, 2)
#' metric_mode(x)
metric_mode = function(x, na.rm = FALSE) {

  # Remove missing values
  if (na.rm) {
    x = x[!is.na(x)]
  }

  # Obtain unique values
  uniq_x = unique(x)

  # Obtain the first occurrence of a maximum value
  # Note: There are no tie-breakers here.
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}
