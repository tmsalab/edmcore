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


#' Metric for Frobenius Norm
#'
#' Computes the Frobenius norm of matrix entries
#'
#' @param estimate Estimated values from the model.
#' @param oracle   Known values used to generate the model.
#' @param na.rm    A `logical` indicating if missing values (including NaN) should
#'                 be removed. Default: `FALSE`
#'
#' @return
#' A single `numeric` value.
#'
#' @seealso
#' [base::norm()]
#'
#' @section Recovery Use:
#' The Frobenius norm is best used to understand differences between
#' the estimated \eqn{\hat\theta} matrix and the oracle \eqn{\theta} matrix.
#'
#' @details
#' The Frobenius norm is an extension of the Euclidean norm to \eqn{\mathcal{K}^{n\times n}}.
#'
#' The metric is computed under:
#' \deqn{\|A\|_{\rm F} = \left(\sum_{i=1}^m \sum_{j=1}^n |a_{ij}|^2\right)^{\frac{1}{2}}}
#'
#' @export
#' @examples
#' # Construct data
#' estimate = matrix(c(1,1,2,4,3,6), nrow = 2, ncol = 3)
#' truth = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)
#'
#' # Compute the frobenius norm
#' metric_frobenius_norm(estimate, truth)
metric_frobenius_norm = function(estimate, oracle, na.rm = FALSE) {

  # Pass computation off to R's norm function
  base::norm(center_values(estimate = estimate, oracle = oracle, na.rm = na.rm),
             type = "F")
}

#' Metric for Element-Wise Accuracy
#'
#' Computes the element-wise accuracy of matrices.
#'
#' @param estimate Estimated values from the model.
#' @param oracle   Known values used to generate the model.
#' @param na.rm    A `logical` indicating if missing values (including NaN) should
#'                 be removed. Default: `FALSE`
#'
#' @return
#' A single `numeric` value between 0 and 1.
#'
#' @seealso
#' [base::norm()]
#'
#' @section Recovery Use:
#' The element-wise recovery metric is best used to understand differences
#' between dichotomous matrices such as the \eqn{\boldsymbol{Q}} and
#' \eqn{\boldsymbol{\Delta}} matrices.
#'
#' @details
#' The element-wise metric is also known as accuracy or
#' the proportion of estimated values that are equivalent to the same elements
#' in the oracle value.
#'
#' The metric is computed under:
#' \deqn{\frac{1}{JK}\sum _{j=1}^J\sum _{k=1}^K\mathcal I(\hat{\theta}_{jk}=\theta_{jk})}
#'
#' @export
#' @examples
#' # Construct data
#' estimate = matrix(c(1,1,2,4,3,6), nrow = 2, ncol = 3)
#' truth = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)
#'
#' # Compute the frobenius norm
#' metric_element_wise(estimate, truth)
metric_element_wise = function(estimate, oracle, na.rm = FALSE) {
  mean(estimate == oracle)
}


#' Metric for Matrix-Wise Accuracy
#'
#' Computes the matrix-wise accuracy.
#'
#' @param estimate Estimated values from the model.
#' @param oracle   Known values used to generate the model.
#' @param na.rm    A `logical` indicating if missing values (including NaN) should
#'                 be removed. Default: `FALSE`
#'
#' @return
#' A single `numeric` value between 0 and 1.
#'
#' @seealso
#' [base::norm()]
#'
#' @section Recovery Use:
#' The element-wise recovery metric is best used to understand differences
#' between dichotomous matrices such as the \eqn{\boldsymbol{Q}} and
#' \eqn{\boldsymbol{\Delta}} matrices.
#'
#' @details
#' The matrix-wise metric is a variant of accuracy that holistically looks
#' at the entire estimated matrix against the entire
#'
#' The metric is computed under:
#' \deqn{I(\hat{\theta}=\theta)}
#'
#' @export
#' @examples
#' # Construct data
#' estimate = matrix(c(1,1,2,4,3,6), nrow = 2, ncol = 3)
#' truth = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)
#'
#' # Compute the frobenius norm
#' metric_matrix_wise(estimate, truth)
metric_matrix_wise = function(estimate, oracle, na.rm = FALSE) {
  1 * isTRUE(all.equal(estimate, oracle))
}

