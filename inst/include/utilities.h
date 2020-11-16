#ifndef EDMCORE_UTILITIES_H
#define EDMCORE_UTILITIES_H

namespace edmcore {

// Comparison ----

template <class T> bool is_needle_in_haystack(T x, unsigned int needle)
{
  arma::uvec m = arma::find(x == needle);

  if (m.n_elem > 0) {
    return true;
  }

  return false;
}


// Sequence construction ----

// Short hand sequence generators with a fixed length.out based on integer
// distance from start to end.
template <class T>
inline T seq_linear_decrease(unsigned int start, unsigned int end) {
  if (start < end) {
    Rcpp::stop("start must be greater than end.");
  }
  return arma::linspace<T>(start, end, start - end + 1);
}

template <class T>
inline T seq_linear_increase(unsigned int start, unsigned int end) {
  if (end < start) {
    Rcpp::stop("end must be greater than start.");
  }
  return arma::linspace<T>(start, end, end - start + 1);
}


// Set differencing ----

inline arma::urowvec set_diff(arma::urowvec x, arma::urowvec y) {
  // Sort and load
  std::vector<unsigned int> x_sorted = arma::conv_to< std::vector<unsigned int> >::from(arma::sort(x));
  std::vector<unsigned int> y_sorted = arma::conv_to< std::vector<unsigned int> >::from(arma::sort(y));
  std::vector<unsigned int> out;

  std::set_difference(x_sorted.begin(), x_sorted.end(),
                      y_sorted.begin(), y_sorted.end(),
                      std::inserter(out, out.end()));

  return arma::conv_to<arma::urowvec>::from(out);
}

// Combination computation ----

// Fast computation of nCk / total number of combination
inline unsigned int n_choose_k( unsigned int n, unsigned int k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

// Heavily modified version from Rosetta code:
// http://rosettacode.org/wiki/Combinations#C.2B.2B
// Enables a generation of all combinations under n_choose_k (nCk) that
// is equivalent to MATLAB's nchoosek(1:N, k) function.
inline arma::umat combination_matrix_from_vector(arma::urowvec x, unsigned int k) {

  unsigned int n = x.size();

  // Pre-allocate space
  arma::umat all_combs(n_choose_k(n, k), k);

  std::string bitmask(k, 1); // K leading 1's
  bitmask.resize(n, 0); // N-K trailing 0's

  unsigned int row_index = 0;
  // Permutate the bitmask
  do {
    // Set the column index to 0
    unsigned int col_index = 0;

    // Process all possible integer pairings.
    // Note: We're interested in i = 1 up to N integers.
    // Need to index shift back for C++ usage
    for (unsigned int i = 0; i < n; ++i)
    {
      // If we hit a 1, then stored the numeric value.
      // Increase column position
      if (bitmask[i]) {
        all_combs(row_index, col_index) = x[i];
        ++col_index;
      }
    }

    // Increase row index as we move to the next permutation
    ++row_index;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  // Return the final amount of combinations
  return all_combs;
}

// Helper function to automatically generate a sequence
inline arma::umat combination_matrix(unsigned int n, unsigned int k) {
  arma::urowvec x = seq_linear_increase<arma::urowvec>(1, n);
  return combination_matrix_from_vector(x, k);
}


//' Generate a binary matrix
//'
//' @param k    dimension of the matrix to compute
//'
//' @details
//' This function mimics the `all_subset = binary(1:(2^K-1), K);` output
//'
//' As an example, consider `k = 3`, we would get back:
//'
//' ```
//' 0 0 1
//' 0 1 0
//' 0 1 1
//' 1 0 0
//' 1 0 1
//' 1 1 0
//' 1 1 1
//' ```
//'
//' Effectively, this is the DtoQ matrix without the first row of zero padding.
// [[Rcpp::export]]
inline arma::mat binary_q_ideal(unsigned int k) {

  // If K drops below 2, this algorithm fails without a special definition.
  // The same holds true for the MATLAB code that was ported.
  if (k == 1) {
    return arma::ones<arma::mat>(1, 1);
  }

  unsigned int nClass =
    static_cast<unsigned int>( pow(2.0, static_cast<double>(k)) - 1 );

  // Create a sequence
  arma::vec x = seq_linear_increase<arma::vec>(1, nClass);

  // Recreate MATLAB approach:
  // x = 1:(2^k - 1)
  // divs = floor(bsxfun(@rdivide, x, base.^((k-1):-1:0)));

  // Construct a divisions matrix
  arma::mat divs(nClass, k);

  // Fill each column with sequence
  divs.each_col() = x;

  // Construct a division sequence
  arma::rowvec division_b(k);
  for(int i = k; i > 0; --i) {
    division_b(k - i) = pow(2.0, static_cast<double>(i - 1));
  }

  // Perform an inplace division on rows
  divs.each_row() /= division_b;

  // Element-wise round down
  divs = floor(divs);

  // We access the 2 column to K columns.
  // We subtract out base times the 1 column to K - 1 columns.
  divs.cols(1, k - 1) -= 2.0*divs.cols(0, k - 2);

  return divs;
}

}

#endif
