#ifndef EDMCORE_ATTRIBUTES_H
#define EDMCORE_ATTRIBUTES_H

#include <RcppArmadillo.h>

arma::vec attribute_bijection(unsigned int K);

arma::vec attribute_inv_bijection(unsigned int K, double CL);

arma::mat attribute_classes(int K);

#endif
