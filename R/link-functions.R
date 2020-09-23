#' Link Functions
#'
#' Functions that rescale the data according to a given distribution.
#'
#' @param x A `matrix` or `vector` containing values to be transformed.
#'
#' @return
#' A `matrix` or `vector` with values transformed under the link function.
#'
#' @seealso [`theta_to_beta`]
#' @rdname link-functions
#' @export
link_logit = function(x) {
  log(x / (1 - x))
}

#' @rdname link-functions
#' @export
link_probit = function(x) {
  qnorm(x)
}

#' Translate Theta values into Beta values
#'
#' Converts theta values underneath a link function to beta values.
#'
#' @param theta         Matrix of theta values
#' @param k             Number of Attributes
#' @param link_function Link function to use
#'
#' @return
#' Beta matrix after translation.
#'
#' @seealso [`link_logit`] and [`link_probit`]
#' @export
#' @rdname theta-to-beta
theta_to_beta = function(theta, k, link_function) {
  link_function(theta) %*% solve(t(GenerateAtable(2 ^ k, k, 2, k)$Atable))
}

#' @export
#' @rdname theta-to-beta
theta_probit_to_beta = function(theta, k) {
  theta_to_beta(theta, k, link_probit)
}

#' @export
#' @rdname theta-to-beta
theta_logit_to_beta = function(theta, k) {
  theta_to_beta(theta, k, link_probit)
}

