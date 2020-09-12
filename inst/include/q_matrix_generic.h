#ifndef EDMCORE_Q_MATRIX_GENERIC_H
#define EDMCORE_Q_MATRIX_GENERIC_H

namespace edmcore {

inline bool check_generic_complete(const arma::mat& Q) {

  // Define variables
  unsigned int
    K = Q.n_cols,
    nClass = static_cast<unsigned int>(pow(2.0, static_cast<double>(K)));

  // Check if Hall's marriage condition is satisfied
  // Check if for any subset S of [K], if |N(S)| >= |S|

  // Get all the 2^K subsets of [K]
  // Effectively, DtoQ table without one row.
  arma::mat all_subset = binary_q_ideal(K);

  // Compute the number of neighbors of each subset of [K]
  for(unsigned int kk = 0; kk < nClass - 1; ++kk) {

    // Subset
    arma::rowvec kk_subset = all_subset.row(kk);

    // 1:(2^K-1)
    // extract columns of Q corresponding to the kk-th subset of attributes
    arma::mat Q_cols = Q.cols(find(kk_subset != 0)); // all_subset(kk,:) ~= 0

    // take the union of items that require the attributes in kk-th subset
    // sum max(Q_cols, [],  2)
    double num_neighbor = arma::sum(arma::max(Q_cols, 1));

    // Escape the loop early if any summed value is not true.
    // Why? The product would be 0 since 1 * 1 * 0 * 1 = 0
    if ( !( num_neighbor >= arma::sum(kk_subset) ) ) {
      return false;
    }
  }

  // If we complete the loop with all true statements, this means we have
  // generic completeness.
  return true;
}


inline bool is_q_generic_identified(const arma::mat& Q) {

  // Unpack Q dimensions
  unsigned int
    J = Q.n_rows,
    K = Q.n_cols;

  // Obtain column sum
  arma::rowvec attr_count = sum(Q, 0);

  // Checks if some attributes is required by <=2 item, then not generic ID  by Theorem 3
  // Verify if Q is generically complete; if not, then not generically identifiable.
  // Ensure J >= 2K+1; if not, conditions not satisfied
  if ( any(attr_count <= 2) || !check_generic_complete(Q) || (J < 2*K + 1)) {
    return false;
  }

  // Create a matrix containing all combinations from j N k
  // Change the matrix to be 0-based indices instead of 1-based indices.
  arma::umat all_submat_K1 = combination_matrix(J, K) - 1;

  // Construct a linear sequence from 0 to J - 1
  arma::urowvec J_item_index = seq_linear_increase<arma::urowvec>(0, J - 1);

  for (unsigned int ii = 0;  ii < all_submat_K1.n_rows; ++ii) {
    arma::urowvec submat_K1_indices = all_submat_K1.row(ii);

    // ii is the index for a K*K submat
    arma::mat Q1 = Q.rows(submat_K1_indices);
    bool Q1_GC = check_generic_complete(Q1);

    // if Q1 is not Gen. Complete, then go to next Q1
    if (Q1_GC == false) {
      continue;
    }

    // if Q1 is Gen. Complete, search for Q2
    arma::urowvec ind_remain = set_diff(J_item_index, submat_K1_indices);

    arma::umat all_submat_K2 = combination_matrix_from_vector(ind_remain, K);

    for(unsigned int jj = 0; jj < all_submat_K2.n_rows; ++jj) {
      arma::urowvec submat_K2_indices = all_submat_K2.row(jj);

      arma::mat Q2 = Q.rows(submat_K2_indices);

      bool Q2_GC = check_generic_complete(Q2);

      if(Q2_GC == false) {
        continue;
      }

      // if two Gen. Complete Q1 and Q2 are found, then find indices
      // of the remaining Q' part, and check Condition E
      arma::urowvec ind_last = set_diff(ind_remain, submat_K2_indices);

      if(!arma::all(arma::sum(Q.rows(ind_last), 0) >= 1)) {
        continue;
      }

      return true;
    }
  }

  // Not identifiable.
  return false;
}

}

#endif
