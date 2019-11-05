#' Generate all attribute permutations possible for each trait.
#'
#' Create all possible orderings for trait columns within the Q matrix.
#'
#' @param K Number of Traits
#'
#' @return
#' K! x K
#'
#' @details
#' This is an internal function wrapper.
#'
#' @noRd
attribute_permutations = function(K) {
  gtools::permutations(K, K)
}


#' Matrix Subscript Position of Highest Matrix Entry
#'
#' Locates the highest matrix entry and returns its location in matrix
#' subscript.
#'
#' @param x A `matrix` with dimensions \eqn{M \times N}{M x N}.
#'
#' @return
#' The location of the highest value reported as a `vector` with entries given
#' as:
#'
#' - `row`
#' - `column`
#'
#' @seealso [`base::which.max()`], [`base::arrayInd()`]
#' @noRd
matrix_max_index = function(x) {
  # Obtain the location of the maximum value in the matrix
  highest_index = which.max(x)

  # Release the subscript
  arrayInd(highest_index, .dim = dim(x))
}


#' Element-wise Accuracy Table for Q Matrix Estimation
#'
#' Given the oracle Q matrix, compare the estimated Q matrix to
#' determine how accurate the estimation routine is.
#'
#' @param oracle    Set containing the actual strategies used to
#'                  generate the data. Dimensions: \eqn{J x K x S}.
#' @param estimated Set containing the estimated strategies values
#'                  on generated data. Dimensions: \eqn{J x K x S}.
#' @param decision  Threshold for applying a classification to estimated
#'                  entry. Default: `0.5`
#'
#' @return
#' A `list` containing:
#'
#' - `element_wise_recovery`:
#'    Average number of entries matching between the Q matrix and estimated
#'    matrix.
#' - `permutated_q`:
#'    Ordered Q matrix by permutation under decision rule of 0.5.
#' - `decision_rule`:
#'    Threshold cut-off applied to the estimated Q matrix to ensure dichotomous
#'    entries.
#'
#' @export
#'
#' @examples
#' # Create a Q matrix
#' Q_oracle =
#'   rbind(c(0, 1),
#'         c(1, 0),
#'         c(0, 1),
#'         c(1, 0),
#'         c(0, 1),
#'         c(0, 1),
#'         c(1, 0))
#'
#' # Simulate a random Q matrix
#' Q_est = matrix(runif(7*2), ncol = 2)
#'
#' #
#' q_recovery_element(Q_oracle, Q_est)
q_recovery_element = function(oracle, estimated, decision = 0.5) {

  if (dim(oracle) != dim(estimated)) {
    msg = sprintf(
      "Dimensions of `oracle` (%s, %s) must match `estimated` (%s, %s).",
        nrow(oracle), ncol(oracle),
        nrow(estimated), ncol(estimated)
      )
    stop(msg, call. = FALSE)
  }

  # Data Dimensions
  K = ncol(oracle)

  # Create permutation tables
  attribute_permutation_table = attribute_permutations(K)

  # Setup a vector to store results from matching columns
  col_match = rep(
    NA, nrow(attribute_permutation_table),
  )

  # Computed averaged Q matrices and perform binary classification
  Q_est_binary = 1 * (estimated > decision)

  # Oracle Q matrices used to generate the data
  Q_oracle = oracle

  # Iterate on trait permutations
  for (r in seq_len(nrow(attribute_permutation_table))) {
    # Obtain the mean difference between the oracle and binary estimate matrix
    col_match[r] =
      mean(Q_oracle == Q_est_binary[, attribute_permutation_table[r,]])
  }

  # Obtain the desired column permutation index
  permutation_idx = which.max(col_match)

  # Retrieve the column ordering
  permutation_column_order = attribute_permutation_table[permutation_idx, ]

  # Retrieve and set the best pairing
  permutated_q = Q_est_binary[, permutation_column_order]

  # Retrieve the element-wise mean
  element_wise_recovery = col_match[r]

  out = structure(
    list(element_wise_recovery = element_wise_recovery,
         permutated_q = permutated_q,
         decision_rule = decision),
    class = c("q_element_recovery", "list")
  )

  return(invisible(out))
}
