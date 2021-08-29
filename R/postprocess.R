
# Latent Class Membership ----

#' Correlation Among Attribute Membership
#'
#' Given the estimated latent class membership probabilities, we compute the
#' correlation across classes.
#'
#' @param pi_hat Estimated latent class membership probabilities
#' @param k      Number of attributes
#' @param order  Number of interactions in the model
#'
#' @examples
#'
#' # Sample calculation when K = 3 and number of classses is 2^3 = 8
#' pi_hat = c(0.10, 0.10, 0.01, 0.13, 0.03, 0.26, 0.05, 0.31)
#' k = 3
#'
#' # Compute attribute correlation matrix
#' postprocess_attribute_correlation(pi_hat, k)
postprocess_attribute_correlation = function(pi_hat, k, order = k) {

  # Obtain the attribute profiles that follow the bijection
  alpha_profiles =
    t(GenerateAtable(
      nClass = 2 ^ k,
      k,
      M = 2,
      order = order
    )$DtoQtable)

  # Obtain the expected value (mean) and standard deviation
  # using formulas for the binomial distribution
  E_alpha = pi_hat %*% alpha_profiles      # p * n
  sd_alpha = sqrt(E_alpha * (1 - E_alpha)) # p * n * q

  # Construct a correlation matrix
  cor_alpha = matrix(0, k, k)

  # Compute under the assumption of an upper triangular matrix
  for (k_1 in seq_len(k - 1)) {
    for (k_2 in (k_1 + 1):k) {
      cov_tmp =
        pi_hat %*% (alpha_profiles[, k_1] * alpha_profiles[, k_2]) - E_alpha[k_1] * E_alpha[k_2]
      cor_alpha[k_1, k_2] =
        cov_tmp / (sd_alpha[k_1] * sd_alpha[k_2])
    }
  }
  # Modify to incorporate
  cor_alpha = cor_alpha + t(cor_alpha) + diag(k)

  # Return correlation matrix
  cor_alpha
}
