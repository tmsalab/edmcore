#include <RcppArmadillo.h>
#include <edmcore.h>

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
//' [attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute bijection ----
//' biject = attribute_bijection(3)
// [[Rcpp::export]]
arma::uvec attribute_bijection(unsigned int K)
{
  return edmcore::attribute_bijection<arma::uvec>(K);
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
//' [attribute_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute inversion bijection ----
//' inv_biject1 = attribute_inv_bijection(5, 1)
//' inv_biject2 = attribute_inv_bijection(5, 2)
// [[Rcpp::export]]
arma::uvec attribute_inv_bijection(unsigned int K, double CL)
{
  return edmcore::attribute_inv_bijection<arma::uvec>(K, CL);
}

//' Generate a vector to map polytomous vector to integers
//'
//' Converts class into a bijection to integers
//'
//' @param M      Number of Responses
//' @param K      Number of Attributes
//'
//' @return
//'
//' Return a matrix containing the class table
//'
//' @export
// [[Rcpp::export]]
arma::uvec attribute_gen_bijection(unsigned int K, unsigned int M)
{
  return edmcore::attribute_gen_bijection<arma::uvec>(K, M);
}

//' Create the K inverse bijection of attribute vectors
//'
//' Converts the class into a bijection.
//'
//' @param K      Number of Attributes
//' @param M      Number of Options.
//' @param CL     Class Number from 0 to (2^K)-1.
//'
//' @return
//'
//' Return a matrix containing the class table
//' @export
// [[Rcpp::export]]
arma::uvec attribute_inv_gen_bijection(unsigned int K, unsigned int M, double CL) {
  return edmcore::attribute_inv_gen_bijection<arma::uvec>(K, M, CL);
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
//' [simcdm::sim_subject_attributes()] and [attribute_inv_bijection()]
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
    alpha_matrix.row(cc) = edmcore::attribute_inv_bijection<arma::vec>(K, cc).t();
  }

  // Release result
  return alpha_matrix;
}

//' Generate tables to store the results during iterations
//'
//' Generate tables to store the results during iterations
//'
//' @param nClass   Number of Attribute Classes
//' @param M        Number of Responses
//' @param K        Number of Attributes
//' @param order    Order is 1, main-effects, 2 main-effects + interactions.
//'               Highest level of interactions you want.
//'
//' @return
//'
//' Return a list containing the tables for different parameters
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List GenerateAtable(unsigned int nClass, unsigned int K, unsigned int M,
                          unsigned int order)
{

  arma::mat FullAtable(nClass, nClass);
  arma::mat FullLBtable(nClass, nClass);
  arma::mat FullDtoQtable(K, nClass);
  arma::mat Fulladjtable(nClass, nClass);
  arma::vec model(nClass);

  // Setup storage for main effect indices
  arma::uvec model_main_effects(nClass);

  // Construct a vv bijection
  arma::uvec vv_bijection = edmcore::attribute_gen_bijection<arma::uvec>(K, M);

  for (unsigned int cr = 0; cr < nClass; cr++) {
    arma::vec alpha_r = edmcore::attribute_inv_gen_bijection<arma::vec>(K, M, cr);

    double nof0s = 0.;
    for (unsigned int k = 0; k < K; k++) {
      nof0s += 1. * (alpha_r(k) == 0);
      FullDtoQtable(k, cr) = 1. * (alpha_r(k) > 0);
    }

    // Checking if the class is in the bijection vector
    // If integer is corresponding with main effect.
    model_main_effects(cr) = 1*edmcore::is_needle_in_haystack(vv_bijection, cr);

    // Flags coefficients with no greater order
    model(cr) = 1. * (nof0s > double(K - order) - 1);
    for (unsigned int cc = 0; cc < nClass; cc++) {
      arma::vec alpha_c = edmcore::attribute_inv_gen_bijection<arma::vec>(K, M, cc);
      double mindiff = arma::min(alpha_r - alpha_c);
      FullAtable(cr, cc) = 1. * (mindiff > -1);
      double maxdiff = arma::accu(abs(alpha_c - alpha_r));
      FullLBtable(cr, cc) = 1. * (maxdiff == 1) * (mindiff < 0) +
        2. * (mindiff > -1) * (maxdiff != 0);
      Fulladjtable(cr, cc) = 1. * (maxdiff == 1);
    }
  }


  // Columns that should be retained under the appropriate model order
  arma::uvec finalcols = find(model == 1);

  // Columns corresponding to the main effects
  // and is present within the model
  arma::uvec maineffectcols = find(model_main_effects == 1 && model == 1);

  arma::mat Atable = FullAtable.cols(finalcols);
  arma::mat LBtable = FullLBtable.cols(finalcols);
  arma::mat DtoQtable = FullDtoQtable.cols(finalcols);
  arma::mat adjtable = Fulladjtable.submat(finalcols, finalcols);

  return Rcpp::List::create(
    Rcpp::Named("Atable", Atable), Rcpp::Named("LBtable", LBtable),
    Rcpp::Named("finalcols", finalcols),
    Rcpp::Named("DtoQtable", DtoQtable), Rcpp::Named("adjtable", adjtable),
    Rcpp::Named("maineffectcols", maineffectcols)
  );
}
