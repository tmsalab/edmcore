#include "attributes.h"
#include <edmcore.h>

// [[Rcpp::export]]
unsigned int n_choose_k(unsigned int n, unsigned int k)
{
    return edmcore::n_choose_k(n, k);
}

// [[Rcpp::export]]
arma::umat combination_matrix(unsigned int n, unsigned int k)
{
    return edmcore::combination_matrix(n, k);
}

// [[Rcpp::export]]
arma::umat combination_matrix_from_vector(arma::urowvec x, unsigned int k)
{
    return edmcore::combination_matrix_from_vector(x, k);
}

// [[Rcpp::export]]
arma::urowvec set_diff_cpp(arma::urowvec x, arma::urowvec y)
{
    return edmcore::set_diff(x, y);
}

// [[Rcpp::export]]
bool is_q_generic_identified(const arma::mat& Q)
{
    return edmcore::is_q_generic_identified(Q);
}

// [[Rcpp::export]]
bool is_q_generic_complete(const arma::mat& Q)
{
    return edmcore::check_generic_complete(Q);
}

// [[Rcpp::export]]
arma::mat binary_q_ideal(unsigned int k)
{
    return edmcore::binary_q_ideal(k);
}

// [[Rcpp::export]]
arma::rowvec seq_linear_increase(unsigned int start, unsigned int end)
{
    return edmcore::seq_linear_increase<arma::rowvec>(start, end);
}


// [[Rcpp::export]]
arma::rowvec seq_linear_decrease(unsigned int start, unsigned int end)
{
    return edmcore::seq_linear_decrease<arma::rowvec>(start, end);
}

//' Verify Q Matrix is Strictly Identifiable
//'
//' Performs a check to see if Q is strictly identified or not.
//'
//' @param Q The Q matrix to be checked with dimensions \eqn{K \times J}{K x J}.
//'
//' @return
//' A `bool` with value either: false or true
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
bool is_strict_q_identified(const arma::mat Q)
{
    return edmcore::is_strict_q_identified(Q);
}

//' Permute indices
//'
//' Constructs and permutes indices within a vector
//'
//' @param nClass Number of classes given by \eqn{2^K}.
//' @param K      Number of Attributes
//' @param order  Degree of Order
//' @param vv     Bijection vector
//' @param perm   Permutations
//' @noRd
//'
//' @return
//' Permuted table at specified indices.
// [[Rcpp::export]]
arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K,
                               unsigned int order, const arma::vec &vv,
                               const arma::vec &perm)
{
    arma::vec vvperm(K);
    arma::vec fullorigperm = arma::linspace(0, nClass - 1, nClass);
    for (unsigned int k = 0; k < K; ++k) {
        vvperm(k) = vv(perm(k));
    }
    arma::vec model(nClass);
    arma::vec fullpermindices(nClass);
    for (unsigned int cr = 0; cr < nClass; ++cr) {
        arma::vec alpha_r = attribute_inv_bijection(K, cr);
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
    arma::vec permindices(origperm.n_elem);
    for (unsigned int p = 0; p < origperm.n_elem; ++p) {
        double origval = origperm(p);
        for (unsigned int pp = 0; pp < origperm.n_elem; ++pp) {
            if (origval == reducedpermindices(pp)) {
                permindices(p) = pp;
            }
        }
    }
    return permindices;
}
