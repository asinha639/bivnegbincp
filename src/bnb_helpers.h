#ifndef BIVNEGBINCP_BNB_HELPERS_H
#define BIVNEGBINCP_BNB_HELPERS_H

#include <RcppArmadillo.h>

double bnb_log_sum_exp(const Rcpp::NumericVector& x);
Rcpp::NumericVector bnb_prefix_sum(const Rcpp::NumericVector& x);
double bnb_lgamma(double x);

#endif
