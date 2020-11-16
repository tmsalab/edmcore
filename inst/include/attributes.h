#ifndef EDMCORE_ATTRIBUTES
#define EDMCORE_ATTRIBUTES

namespace edmcore {

// Binary Attribute bijections  ----

template <typename T>
inline T attribute_bijection(unsigned int K)
{
  T vv(K);
  for (unsigned int k = 0; k < K; k++) {
    vv(k) = static_cast<typename T::elem_type>(std::pow(2.0, static_cast<double>(K - k) - 1.0));
  }
  return vv;
}


template <typename T>
inline T attribute_inv_bijection(unsigned int K, double CL)
{
  T alpha(K);
  for (unsigned int k = 0; k < K; k++) {

    double twopow = std::pow(2.0, static_cast<double>(K - k) - 1.0);
    alpha(k) = static_cast<typename T::elem_type>(twopow <= CL);
    CL = CL - twopow * alpha(k);
  }
  return alpha;
}

// General Attribute bijections ----
template <typename T>
inline T attribute_gen_bijection(unsigned int K, unsigned int M)
{
  T vv(K);
  for (unsigned int k = 0; k < K; k++) {
    vv(k) = static_cast<typename T::elem_type>(
      pow(static_cast<double>(M), static_cast<double>(K - k) - 1.0)
    );
  }
  return vv;
}

template <typename T>
inline T attribute_inv_gen_bijection(unsigned int K, unsigned int M, unsigned int CL)
{
  T alpha(K);
  for (unsigned int k = 0; k < K; k++) {
    unsigned int Mpow = static_cast<unsigned int>(pow(M, K - k - 1.0));
    unsigned int ak = 0.;
    while (((ak + 1) * Mpow <= CL) && (ak < M)) {
      ak += 1;
    }
    //Rcpp::Rcout << "k: " << k << "cl: " << CL << " with: " << ak << std::endl;
    alpha(k) = static_cast<typename T::elem_type>(ak);
    CL = CL - Mpow * alpha(k);
  }
  return alpha;
}



inline arma::umat permutate_attribute_level_table(unsigned int K, unsigned int M = 2) {
  unsigned int nClass =
    static_cast<unsigned int>(pow(2.0, static_cast<double>(K)));

  // Create a 2^K by K table of attribute classes
  arma::umat all_binary_attribute_classes = arma::zeros<arma::umat>(nClass, K);

  for (unsigned int cc = 1; cc < nClass; ++cc) {
    all_binary_attribute_classes.row(cc) =
      attribute_inv_gen_bijection<arma::urowvec>(K, M, cc);
  }

  // Establish a vector to map between binary classes and integers
  arma::uvec vv(K);
  for(unsigned int i = 0; i < K; ++i) {
    vv(i) = static_cast<unsigned int>(pow(2.0, static_cast<double>(i)));
  }

  std::reverse(vv.begin(), vv.end());

  // Creating a 2^K by 2^K matrix
  // Each row corresponds with focal attribute classes
  // (e.g., true data generating arrangement)
  // Each column indicates whether an attribute level is swapped
  // 0 = not swapped, 1 = swapped
  // e.g., 000 = no attributes levels swapped;
  // 011 = 2nd and 3rd attribute levels swapped
  arma::umat level_swap_table = arma::zeros<arma::umat>(nClass, nClass);

  for (unsigned int perms = 0; perms < nClass; ++perms) {
    //  Identify which attribute levels are swapped
    arma::uvec binary_levels_swapped = attribute_inv_gen_bijection<arma::uvec>(K, M, perms);
    arma::uvec attributes_with_swapped_levels = arma::find(binary_levels_swapped == 1);

    // Create a temporary table of swapped attributes modify
    arma::umat temp_attributes = all_binary_attribute_classes;

    temp_attributes.cols(attributes_with_swapped_levels) =
      1 - all_binary_attribute_classes.cols(attributes_with_swapped_levels);

    level_swap_table.col(perms) = temp_attributes * vv;

  }

  // Release table
  return level_swap_table;
}

inline arma::uvec permuteAtableIndices(unsigned int nClass, unsigned int K,
                                       unsigned int M, unsigned int order,
                                       const arma::vec &vv, const arma::uvec &perm)
{
  arma::vec vvperm(K);
  arma::vec fullorigperm = arma::linspace(0, nClass - 1, nClass);
  for (unsigned int k = 0; k < K; k++) {
    vvperm(k) = vv(perm(k));
  }
  arma::vec model(nClass);
  arma::vec fullpermindices(nClass);
  for (unsigned int cr = 0; cr < nClass; cr++) {
    arma::vec alpha_r = attribute_inv_gen_bijection<arma::vec>(K, M, cr);
    double nof0s = 0.;
    for (unsigned int k = 0; k < K; k++) {
      nof0s += 1. * (alpha_r(k) == 0);
    }
    model(cr) = 1. * (nof0s > double(K - order) - 1);
    arma::vec alpha_perm(K);
    fullpermindices(cr) = arma::accu(alpha_r % vvperm);
  }
  arma::uvec finalcols = find(model == 1);
  arma::vec origperm = fullorigperm(finalcols);
  arma::vec reducedpermindices = fullpermindices(finalcols);
  arma::uvec permindices(origperm.n_elem);
  for (unsigned int p = 0; p < origperm.n_elem; p++) {
    double origval = origperm(p);
    for (unsigned int pp = 0; pp < origperm.n_elem; pp++) {
      if (origval == reducedpermindices(pp)) {
        permindices(p) = pp;
      }
    }
  }
  return permindices;
}

inline arma::vec mad_row_and_column(arma::mat estimate, arma::mat truth, arma::umat col_perms,
                                    unsigned int K, unsigned int M, unsigned int order) {

  unsigned int nClass = static_cast<unsigned int>(
    pow(static_cast<double>(M), static_cast<double>(K))
  );

  // Under non-monotonicity
  // Find if the betas need to swap or not
  // col_perms = gtools::permutations(K,K) // passed in for now

  // Create permutation tables
  arma::umat permutation_table = permutate_attribute_level_table(K, M);

  arma::vec vv = attribute_gen_bijection<arma::vec>(K, M);

  unsigned int n_permutations = permutation_table.n_rows;

  arma::mat col_match = arma::mat(col_perms.n_rows, n_permutations);

  for(unsigned int r = 0; r < col_perms.n_rows; ++r) {

    for (unsigned int c = 0; c < n_permutations; ++c) {

      arma::uvec perm = col_perms.row(r).t()  - 1;
      // Obtain a new order based on the column swaps
      arma::uvec border = permuteAtableIndices(nClass,
                                               K,
                                               M,
                                               order,
                                               vv,
                                               perm
      );

      // Obtain the permuted order based on the attribute level swap
      arma::urowvec col_indices = permutation_table.row(c); // + 1

      // Compute the MAD value
      arma::uvec bordern = border.elem(col_indices);
      col_match(r, c) = arma::mean(
        arma::mean(arma::abs(estimate.cols(bordern) - truth))
      );
    }
  }

  // Obtain where the lowest element is in the matrix
  arma::uword minimum_value_element_id = col_match.index_min();
  arma::umat colperm_ind = ind2sub(size(col_match), minimum_value_element_id);
  unsigned int column_swap_ind = colperm_ind(0); // row location
  unsigned int level_swap_ind = colperm_ind(1); // column location

  arma::vec out(3);
  out(0) = col_match(minimum_value_element_id); // mean absolute deviation at the lowest point
  out(1) = column_swap_ind;
  out(2) = level_swap_ind;

  return out;


  // // Convert from element location to matrix subscript
  // arma::umat colperm_ind = ind2sub(size(col_match), minimum_value_element_id);
  // unsigned int column_swap_ind = colperm_ind(0); // row location
  // unsigned int level_swap_ind = colperm_ind(1); // column location
  //
  // // Retrieve permutation index for Q matrix
  // // for beta and theta matrices you need to use two swap approach below
  // arma::uvec permind_column = col_perms.row(column_swap_ind).t()  - 1;
  //
  // // Obtain the attribute class permutation from the column swaps
  // arma::uvec border_class_column =
  //   permuteAtableIndices(nClass,
  //                        K,
  //                        M,
  //                        order,
  //                        vv,
  //                        permind_column
  //   );
  //
  // // Obtain the attribute class permutation from the best-level swaps
  // arma::urowvec border_class_level =  permutation_table.row(level_swap_ind);
  //
  // // MAD for beta
  // arma::mat estimate_swap_class_level = estimate.cols(border_class_level);
  // double mean_abs_deviation = arma::mean(
  //   arma::mean(
  //     arma::abs(estimate_swap_class_level.cols(border_class_column) - truth)
  //   )
  // );

}

inline double mad_iteration_row_and_column(arma::cube estimate,
                                           arma::mat truth,
                                           arma::umat col_perms,
                                           unsigned int K, unsigned int M, unsigned int order) {

  unsigned int nClass = static_cast<unsigned int>(
    pow(static_cast<double>(M), static_cast<double>(K))
  );

  // Under non-monotonicity
  // Find if the betas need to swap or not
  // col_perms = gtools::permutations(K,K) // passed in for now

  // Create permutation tables
  arma::umat permutation_table = permutate_attribute_level_table(K, M);

  arma::vec vv = attribute_gen_bijection<arma::vec>(K, M);

  unsigned int n_iterations = estimate.n_slices;
  unsigned int n_permutations = permutation_table.n_rows;
  unsigned int n_column_swaps = col_perms.n_rows;

  arma::mat col_match = arma::mat(n_column_swaps, n_permutations);
  arma::vec mad_value(n_iterations);

  arma::mat running_sum_swapped_estimate = arma::zeros<arma::mat>(estimate.n_rows, estimate.n_cols);

  for(unsigned int i = 0; i < n_iterations; ++i) {
    arma::mat estimate_slice = estimate.slice(i);

    for(unsigned int r = 0; r < n_column_swaps; ++r) {

      for (unsigned int c = 0; c < n_permutations; ++c) {

        arma::uvec perm = col_perms.row(r).t()  - 1;

        // Obtain a new order based on the column swaps
        arma::uvec border = permuteAtableIndices(nClass,
                                                 K,
                                                 M,
                                                 order,
                                                 vv,
                                                 perm
        );

        // Obtain the permuted order based on the attribute level swap
        arma::urowvec col_indices = permutation_table.row(c); // + 1

        // Re-order the border by column swaps
        arma::uvec bordern = border.elem(col_indices);

        // Compute the Mean Absolute Deviation
        col_match(r, c) = arma::mean(
          arma::mean(arma::abs(estimate_slice.cols(bordern) - truth))
        );
      }
    }

    // Obtain where the lowest element is in the matrix
    arma::uword minimum_value_element_id = col_match.index_min();
    arma::umat colperm_ind = ind2sub(size(col_match), minimum_value_element_id);
    unsigned int column_swap_ind = colperm_ind(0); // row location
    unsigned int level_swap_ind = colperm_ind(1); // column location

    // Retrieve permutation index for Q matrix
    // for beta and theta matrices you need to use two swap approach below
    arma::uvec permind_column = col_perms.row(column_swap_ind).t()  - 1;

    // Obtain the attribute class permutation from the column swaps
    arma::uvec border_class_column =
      permuteAtableIndices(nClass,
                           K,
                           M,
                           order,
                           vv,
                           permind_column
      );

    // Obtain the attribute class permutation from the best-level swaps
    arma::urowvec border_class_level =  permutation_table.row(level_swap_ind);

    // Need to swap once temporarily
    arma::mat estimate_swap_class_level = estimate_slice.cols(border_class_level);

    running_sum_swapped_estimate += estimate_swap_class_level.cols(border_class_column);
  }

  arma::mat mean_swap_estimate =
    1.0/static_cast<double>(n_iterations) * running_sum_swapped_estimate;

  double mean_abs_deviation_iteration = arma::mean(
    arma::mean(
      arma::abs(mean_swap_estimate - truth)
    )
  );

  return mean_abs_deviation_iteration;
}

}// end edmcore namespace

#endif
