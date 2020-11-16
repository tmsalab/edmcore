#' Link Functions
#'
#' Functions that rescale the data according to a given distribution.
#'
#' @param x A `matrix` or `vector` containing quantiles to be transformed.
#' @param p A `matrix` or `vector` containing probabilities to be transformed.
#'
#' @return
#' A `matrix` or `vector` with values transformed under the link function.
#'
#' @seealso [`theta_to_beta`]
#' @rdname link-functions
#' @export
link_logit = function(p) {
  # log(p / (1 - p))
  stats::qlogis(p)
}

#' @rdname link-functions
#' @export
link_logit_inv = function(x) {
  # exp(p) / (1 + exp(x))
  stats::plogis(x)
}

#' @rdname link-functions
#' @export
link_probit = function(p) {
  stats::qnorm(p)
}

#' @rdname link-functions
#' @export
link_probit_inv = function(x) {
  stats::pnorm(x)
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

