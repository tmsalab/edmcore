#include <RcppArmadillo.h>

//' Constructs Unique Attribute Pattern Map
//'
//' Computes the powers of 2 from \eqn{0} up to \eqn{K - 1} for
//' \eqn{K}-dimensional attribute pattern.
//'
//' @param K  Number of Attributes.
//'
//' @return
//' A \code{vec} with length \eqn{K} detailing the power's of 2.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [ecdmcore::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute bijection ----
//' biject = attribute_bijection(3)
// [[Rcpp::export]]
arma::vec attribute_bijection(unsigned int K)
{
  arma::vec vv(K);
  for (unsigned int i = 0; i < K; ++i) {
    vv(i) = std::pow(2.0, static_cast<double>(K - i) - 1.0);
  }
  return vv;
}

//' Perform an Inverse Bijection of an Integer to Attribute Pattern
//'
//' Convert an integer between \eqn{0} and \eqn{2^{K-1}} to
//' \eqn{K}-dimensional attribute pattern.
//'
//' @param CL An `integer` between \eqn{0} and \eqn{2^{K-1}}
//' @inheritParams attribute_bijection
//'
//' @return
//' A \eqn{K}-dimensional vector with an attribute pattern corresponding
//' to `CL`.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [ecdmcore::attribute_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute inversion bijection ----
//' inv_biject1 = attribute_inv_bijection(5, 1)
//' inv_biject2 = attribute_inv_bijection(5, 2)
// [[Rcpp::export]]
arma::vec attribute_inv_bijection(unsigned int K, double CL)
{
  arma::vec alpha(K);

  for (unsigned int k = 0; k < K; ++k) {
    double twopow = std::pow(2.0, static_cast<double>(K - k) - 1.0);
    alpha(k) = (twopow <= CL);
    CL = CL - twopow * alpha(k);
  }

  return alpha;
}

//' Simulate all the Latent Attribute Profile \eqn{\mathbf{\alpha}_c} in
//' Matrix form
//'
//' Generate the \eqn{\mathbf{\alpha}_c = (\alpha_{c1}, \ldots, \alpha_{cK})'}
//' attribute profile matrix for members of class \eqn{c} such that
//' \eqn{\alpha_{ck}} ' is 1 if members of class \eqn{c} possess skill \eqn{k}
//' and zero otherwise.
//'
//' @param K Number of Attributes
//'
//' @return
//' A \eqn{2^K} by \eqn{K} `matrix` of latent classes
//' corresponding to entry \eqn{c} of \eqn{pi} based upon
//' mastery and nonmastery of the \eqn{K} skills.
//'
//' @author
//' James Joseph Balamuta and Steven Andrew Culpepper
//'
//' @seealso
//' [simcdm::sim_subject_attributes()] and [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Simulate Attribute Class Matrix ----
//'
//' # Define number of attributes
//' K = 3
//'
//' # Generate an Latent Attribute Profile (Alpha) Matrix
//' alphas = attribute_classes(K)
// [[Rcpp::export]]
arma::mat attribute_classes(int K)
{
  // Modified version of ETAMatrix

  // Calculate number of classes
  unsigned int nClass =
    static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

  // Create alpha matrix
  arma::mat alpha_matrix(nClass, K);

  // Fill alpha matrix with classes under an inverse bijection
  for (unsigned int cc = 0; cc < nClass; cc++) {
    alpha_matrix.row(cc) = attribute_inv_bijection(K, cc).t();
  }

  // Release result
  return alpha_matrix;
}
