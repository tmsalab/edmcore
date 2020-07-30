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


## Numerical approaches -----

center_values = function(estimate, oracle, na.rm = FALSE) {
  # Take difference across array
  base::sweep(x = estimate,
              MARGIN = c(1, 2),
              STATS = oracle, FUN = "-", na.rm = na.rm)
}

#' Metric for Bias
#'
#' Computes the Bias
#'
#' @param estimate Estimated values from the model.
#' @param oracle   Known values used to generate the model.
#' @param na.rm    A `logical` indicating if missing values (including NaN) should
#'                 be removed. Default: `FALSE`
#'
#' @return
#' A `numeric` value for each parameter comparison.
#'
#' @seealso
#' [base::norm()]
#'
#' @section Recovery Use:
#' The bias may be used to understand the difference between
#' \eqn{\hat\theta} matrix and the oracle \eqn{\theta} matrix.
#'
#' @details
#' The bias measures the difference between expected value of an estimated
#' parameter and its true value.
#'
#' The metric is computed under:
#' \deqn{\operatorname{Bias}(\hat \theta, \theta) = E\left[\hat\theta\right] - \theta}
#'
#' @export
#' @examples
#' # Construct data
#' estimate = matrix(c(1,1,2,4,3,6), nrow = 2, ncol = 3)
#' truth = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)
#'
#' # Compute the bias
#' metric_bias(estimate, truth)
metric_bias = function(estimate, oracle, na.rm = FALSE) {
  mean(estimate, na.rm = na.rm) - oracle
}
