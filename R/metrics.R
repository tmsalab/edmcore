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
attribute_permutation_table = function(K) {
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

#' Reorder an object to match a target as closely as possible
#'
#' Rearranges columns within a matrix until the closest permutation
#' to the target is found.
#'
#' @param current Present `matrix` object.
#' @param target  Desired `matrix` object.
#'
#' @return
#' A `permutated_matrix` object that contains:
#'
#' - `permutate_matrix`:
#'    The `current` matrix permutated to match the `target`
#' - `permutate_order`:
#'    The permutation order applied to the columns of `current`.
#' - `permutate_mean`:
#'    The mean difference between the permutation current and target.
#'    Bound between 0 and 1 with a mean value close to 1 indicating an
#'    exact permutation between `current` and `target` is possible.
#' - `permutate_target`:
#'    Reference object used for permutation.
#'
#' @export
permutate_binary_matrix = function(current, target) {

  if (all(dim(current) != dim(target))) {
    msg = sprintf(
      "Dimensions of `current` (%s, %s) must match `target` (%s, %s).",
      nrow(current), ncol(current),
      nrow(target), ncol(target)
    )
    stop(msg, call. = FALSE)
  }

  # Data Dimensions
  K = ncol(current)

  # Create permutation tables
  permutation_table = attribute_permutation_table(K)

  n_permutations = nrow(permutation_table)

  # Setup a vector to store results from matching columns
  col_distances = rep(NA, n_permutations)

  # Iterate on trait permutations
  for (r in seq_len(n_permutations)) {
    # Retrieve permutation indices
    col_indices = permutation_table[r,]

    # Obtain the mean difference between the current and target matrix
    col_distances[r] = mean(current[, col_indices] == target)
  }

  # Obtain the desired column permutation index
  permutation_idx = which.max(col_distances)

  # Release the best permutation order
  permutation_order = permutation_table[permutation_idx,]

  out = structure(list(
    permutate_matrix = current[, permutation_order],
    permutate_order  = permutation_order,
    permutate_mean   = col_distances[permutation_idx],
    permutate_target = target
  ), class = "permutation_matrix")

  return(out)
}

#' Custom print handle for displaying the permutation ordering
#'
#' @param x   A `permutation_matrix` object.
#' @param ... Not used.
#'
#' @export
print.permutation_matrix = function(x, ...) {
  cat("Results of permutating binary matrix...")
  cat("Column Permutation:", x$permutation_order, "\n")
  cat("Average Permutation Element Match:", x$permutation_mean, "\n\n")
  cat("Permutated Matrix:\n")
  print(x$permutate_matrix)
}

#' Element-wise Accuracy Table for Q Matrix Estimation
#'
#' Given the oracle Q matrix, compare the estimated Q matrix to
#' determine how accurate the estimation routine is.
#'
#' @param estimated Set containing the estimated strategies values
#'                  on generated data. Dimensions: \eqn{J x K x S}.
#' @param oracle    Set containing the actual strategies used to
#'                  generate the data. Dimensions: \eqn{J x K x S}.
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
#' # Obtain the recovery metric for elementwise matches
#' recovery_element(Q_est, Q_oracle)
recovery_element = function(estimated, oracle, decision = 0.5) {

  # Perform the binary classification
  Q_est_binary = 1 * (estimated > decision)

  # Rename oracle
  Q_oracle = oracle

  permutation = permutate_binary_matrix(Q_est_binary, Q_oracle)

  return(permutation)
}
