#' Permutation Attribute Table
#'
#' Generate a table containing attribute permutations.
#'
#' @param k Supply the number of attributes to change
#' @param m Number of Responses. Default: 2.
#' @examples
#'
#' # Example of generating all attribute level swap permutations
#' k = 3
#' nClass = 2 ^ k
#'
#' # Total number of label swaps for unstructured mixture
#' factorial(nClass)
#'
#' # Total number of label swaps for structured mixture
#' factorial(k) * (2 ^ k)
#'
#' # Create the permutation table
#' permutation_table = permutate_attribute_level_table(k = 3, m = 2)
#'
#' # Loop over columns to find the positions of the equivalent attribute level swaps
#' @export
permutate_attribute_level_table = function(k, m = 2) {

  nClass = 2 ^ k

  # Create a 2^K by K table of attribute classes
  all_binary_attribute_classes = matrix(0, nClass, k)
  for (cc in 1:(nClass - 1)) {
    all_binary_attribute_classes[cc + 1, ] =
      t(attribute_inv_gen_bijection(k, m, cc))
  }

  # Establish a vector to map between binary classes and integers
  vv <- 2 ^ {
    (k:1) - 1
  }

  # Creating a 2^K by 2^K matrix
  # Each row corresponds with focal attribute classes
  # (e.g., true data generating arrangement)
  # Each column indicates whether an attribute level is swapped
  # 0 = not swapped, 1 = swapped
  # e.g., 000 = no attributes levels swapped;
  # 011 = 2nd and 3rd attribute levels swapped

  level_swap_table = matrix(0, nClass, nClass)
  for (perms in 0:(nClass - 1)) {
    # Identify which attribute levels are swapped
    binary_levels_swapped = attribute_inv_gen_bijection(k, m, perms)
    attributes_with_swapped_levels = which(binary_levels_swapped == 1)

    # Create a temporary table of swapped attributes modify
    temp_attributes <- all_binary_attribute_classes
    temp_attributes[, attributes_with_swapped_levels] =
      1 - all_binary_attribute_classes[, attributes_with_swapped_levels]
    level_swap_table[, perms + 1] <- temp_attributes %*% vv
  }

  # Release table
  level_swap_table
}

#' Reorganize parameters in Theta Matrix
#'
#' @param estimated_theta  Estimated Theta Matrix
#' @param oracle_theta     Known Theta Matrix
#' @param k                Number of Attributes
#'
#' @return
#' A vector containing the permutation order
#' @export
permutate_theta_order = function(estimated_theta, oracle_theta, k = 3) {

  # Create permutation table
  permutation_table = gtools::permutations(k, k)
  n_permutations = nrow(permutation_table)

  # Setup a vector to store results from matching columns
  distance_by_column = rep(NA, n_permutations)

  # Construct a bijection vector into M = 2
  vv = attribute_bijection(k)

  # Iterate on trait permutations
  for (r in seq_len(n_permutations)) {
    # Retrieve permutation indices
    col_indices = permutation_table[r,]

    beta_order = c(permuteAtableIndices(
      nClass = 2 ^ k,
      k,
      order = k,
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
    nClass = 2 ^ k,
    k,
    order = k,
    vv,
    permutation_order - 1 # remove 1 from indices
  )) + 1

  # Define a permutation search
  structure(
    list(
      "m" = 2,
      "k" = k,
      "permutation_order" = permutation_order,
      "beta_order" = beta_order
    ),
    class = "permutate_order")
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
