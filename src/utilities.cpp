#include <ecdmcore.h>
#include "attributes.h"

//' Verify Q Matrix is Identifiable
//'
//' Performs a check to see if Q is identifable or not.
//'
//' @param Q The Q matrix to be checked with dimensions \eqn{K \times J}{K x J}.
//'
//' @return
//' A `bool` with value either: false or true
//'
//' @export
// [[Rcpp::export]]
bool check_identifiability(const arma::mat Q)
{
  return ecdmcore::check_identifiability(Q);
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
                               unsigned int order, const arma::vec& vv,
                               const arma::vec& perm){
  arma::vec vvperm(K);
  arma::vec fullorigperm=arma::linspace(0,nClass-1,nClass);
  for(unsigned int k=0;k<K;++k){
    vvperm(k)=vv(perm(k));
  }
  arma::vec model(nClass);
  arma::vec fullpermindices(nClass);
  for(unsigned int cr=0;cr<nClass;++cr){
    arma::vec alpha_r = attribute_inv_bijection(K,cr);
    double nof0s=0.;
    for(unsigned int k=0;k<K;k++){
      nof0s+=1.*(alpha_r(k)==0);
    }
    model(cr)=1.*(nof0s>double(K-order)-1);
    arma::vec alpha_perm(K);
    fullpermindices(cr)=arma::accu(alpha_r%vvperm);
  }
  arma::uvec finalcols = find(model == 1);
  arma::vec origperm=fullorigperm(finalcols);
  arma::vec reducedpermindices=fullpermindices(finalcols);
  arma::vec permindices(origperm.n_elem);
  for(unsigned int p=0;p<origperm.n_elem;++p){
    double origval=origperm(p);
    for(unsigned int pp=0;pp<origperm.n_elem;++pp){
      if(origval==reducedpermindices(pp)){
        permindices(p)=pp;
      }
    }
  }
  return permindices;
}
