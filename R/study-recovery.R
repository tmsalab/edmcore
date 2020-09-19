#' Reorganize parameters in Theta Matrix
#'
#' @param estimated_theta  Estimated Theta Matrix
#' @param oracle_theta     Known Theta Matrix
#' @param K                Number of Attributes
#'
#' @return
#' A vector containing the permutation order
#' @export
permutate_order = function(estimated_theta, oracle_theta, K = 3) {

  # Create permutation table
  permutation_table = gtools::permutations(K, K)
  n_permutations = nrow(permutation_table)

  # Setup a vector to store results from matching columns
  distance_by_column = rep(NA, n_permutations)

  # Construct a bijection vector into M = 2
  vv = attribute_bijection(K)

  # Iterate on trait permutations
  for (r in seq_len(n_permutations)) {
    # Retrieve permutation indices
    col_indices = permutation_table[r,]

    beta_order = c(permuteAtableIndices(
      nClass = 2 ^ K,
      K,
      order = K,
      vv,
      col_indices - 1 # remove 1 from indices
    )) + 1

    # Obtain the mean difference between the current and target matrix
    distance_by_column[r] = mean(abs(estimated_theta[, beta_order] - oracle_theta))
  }

  # Obtain the desired column permutation index
  permutation_idx = which.min(distance_by_column)

  # Determine the best permutation order from 1 to K
  permutation_order = permutation_table[permutation_idx,]

  beta_order = c(permuteAtableIndices(
    nClass = 2 ^ K,
    K,
    order = K,
    vv,
    permutation_order - 1 # remove 1 from indices
  )) + 1

  # Define a permutation search
  structure(
    list(
      "m" = 2,
      "k" = K,
      "permutation_order" = permutation_order,
      "beta_order" = beta_order
       ),
    class = "permutate_order")
}


