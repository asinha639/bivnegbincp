#include "bnb_helpers.h"

double bnb_lgamma(double x) {
  return R::lgammafn(x);
}

double bnb_log_sum_exp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  if (n == 0) {
    return R_NegInf;
  }

  double max_value = R_NegInf;
  for (int i = 0; i < n; ++i) {
    if (Rcpp::NumericVector::is_na(x[i])) {
      return NA_REAL;
    }
    if (x[i] > max_value) {
      max_value = x[i];
    }
  }

  if (!R_FINITE(max_value)) {
    return max_value;
  }

  double sum_exp = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_exp += std::exp(x[i] - max_value);
  }

  return max_value + std::log(sum_exp);
}

Rcpp::NumericVector bnb_prefix_sum(const Rcpp::NumericVector& x) {
  const int n = x.size();
  Rcpp::NumericVector out(n);
  double running = 0.0;

  for (int i = 0; i < n; ++i) {
    running += x[i];
    out[i] = running;
  }

  return out;
}
