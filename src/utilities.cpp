#include <RcppArmadillo.h>

//' Verify Q Matrix is Identifiable
//'
//' Performs a check to see if Q is identifable or not.
//'
//' @param Q The Q matrix to be checked with dimensions \eqn{K \times J}{K x J}.
//'
//' @return A `bool` with value either: false or true
//' @export
// bool check_identifiability(const arma::mat Q)
// {
//   return ecdmcore::check_identifiability(Q);
// }
