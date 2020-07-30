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
#' @export
#' @examples
#' # Sample class-data
#' x = c(0, 0, 0, 1, 1, 2)
#' metric_mode(x)
metric_mode = function(x, na.rm = FALSE) {

  if (na.rm) {
    x = x[!is.na(x)]
  }
  uniq_x = unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}
