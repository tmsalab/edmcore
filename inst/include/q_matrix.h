#ifndef EDMCORE_Q_MATRIX_H
#define EDMCORE_Q_MATRIX_H

namespace edmcore {

inline bool is_strict_q_identified(const arma::mat Q) {
    unsigned int K = Q.n_cols;
    unsigned int J = Q.n_rows;

    arma::mat ones_zero_on_diag = -1 * arma::ones<arma::mat>(K, K);
    arma::vec zeros_K = arma::zeros<arma::vec>(K);
    ones_zero_on_diag.diag() = zeros_K;

    arma::vec c_sum = (arma::sum(Q, 0)).t();
    arma::vec r_sum = arma::sum(Q, 1);
    arma::mat I_check = Q * ones_zero_on_diag;
    arma::mat I_count = arma::zeros<arma::mat>(J, K);
    I_count.elem(arma::find(I_check > -1)).fill(1.0);
    arma::vec n_ek = (arma::sum(I_count, 0)).t();

    double min_c = (arma::min(c_sum) > 2);
    double min_r = (arma::min(r_sum) > 0);
    double min_ek = (arma::min(n_ek) > 1);

    return (min_c + min_r + min_ek > 2);
}

}

#endif
