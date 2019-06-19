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
