#ifndef EDMCORE_PERMUTATE_A_TABLE
#define EDMCORE_PERMUTATE_A_TABLE
inline arma::vec inv_gen_bijectionvector(unsigned int K, unsigned int M, double CL)
{
  arma::vec alpha(K);
  for (unsigned int k = 0; k < K; k++) {
    double Mpow = pow(M, K - k - 1);
    double ak = 0.;
    while (((ak + 1.) * Mpow <= CL) & (ak < M)) {
      ak += 1.;
    }
    alpha(k) = ak;
    CL = CL - Mpow * alpha(k);
  }
  return alpha;
}

inline arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K,
                               unsigned int M, unsigned int order,
                               const arma::vec &vv, const arma::vec &perm)
{
    arma::vec vvperm(K);
    arma::vec fullorigperm = arma::linspace(0, nClass - 1, nClass);
    for (unsigned int k = 0; k < K; k++) {
        vvperm(k) = vv(perm(k));
    }
    arma::vec model(nClass);
    arma::vec fullpermindices(nClass);
    for (unsigned int cr = 0; cr < nClass; cr++) {
        arma::vec alpha_r = inv_gen_bijectionvector(K, M, cr);
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

#endif
