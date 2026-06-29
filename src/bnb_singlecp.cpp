#include "bnb_helpers.h"

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericVector;

static double log_product_constant(double sx1, double a) {
  return bnb_lgamma(sx1 + 1.0) +
    bnb_lgamma(a + 2.0) -
    bnb_lgamma(a + sx1 + 2.0);
}

static double log_segment_constant(double sx1,
                                   double sx2,
                                   double segment_length,
                                   double kappa_const) {
  const double a = segment_length * kappa_const;

  return -std::log(a + 1.0) +
    bnb_lgamma(sx1 + a + 2.0) +
    bnb_lgamma(sx2 + 1.0) -
    bnb_lgamma(sx2 + sx1 + a + 3.0) +
    log_product_constant(sx1, a);
}

static double log_product_dirichlet_first(double sx1,
                                          double a,
                                          double v1,
                                          double v3) {
  return bnb_lgamma(sx1 + v1) +
    bnb_lgamma(a + v3 + 1.0) -
    bnb_lgamma(a + sx1 + v1 + v3);
}

static double log_product_dirichlet_second(double sx1,
                                           double a,
                                           double v1,
                                           double v3) {
  return bnb_lgamma(sx1 + v1) +
    bnb_lgamma(a + v3 - v1 + 1.0) -
    bnb_lgamma(a + sx1 + v3);
}

static double log_segment_dirichlet_first(double sx1,
                                          double sx2,
                                          double segment_length,
                                          double kappa_const,
                                          double v1,
                                          double v2,
                                          double v3) {
  const double a = segment_length * kappa_const;

  return -std::log(a + v3) +
    bnb_lgamma(sx1 + a + v1 + v3) +
    bnb_lgamma(sx2 + v2) -
    bnb_lgamma(sx2 + sx1 + a + v1 + v2 + v3) +
    log_product_dirichlet_first(sx1, a, v1, v3);
}

static double log_segment_dirichlet_second(double sx1,
                                           double sx2,
                                           double segment_length,
                                           double kappa_const,
                                           double v1,
                                           double v2,
                                           double v3) {
  const double a = segment_length * kappa_const;

  return -std::log(a + v3) +
    bnb_lgamma(sx1 + a + v1 + v3) +
    bnb_lgamma(sx2 + v2) -
    bnb_lgamma(sx2 + sx1 + a + v1 + v2 + v3) +
    log_product_dirichlet_second(sx1, a, v1, v3);
}

static List finalize_singlecp(const NumericVector& log_post,
                              int n) {
  const int n_tau = log_post.size();
  NumericVector post_prob(n_tau);
  const double normalizer = bnb_log_sum_exp(log_post);

  if (!R_FINITE(normalizer)) {
    return List::create(
      Rcpp::Named("post_prob") = post_prob,
      Rcpp::Named("max_post_prob") = NA_REAL,
      Rcpp::Named("cp_index") = NA_INTEGER,
      Rcpp::Named("cp_status") = NA_INTEGER
    );
  }

  double max_post = R_NegInf;
  int cp_index = NA_INTEGER;

  for (int i = 0; i < n_tau; ++i) {
    post_prob[i] = std::exp(log_post[i] - normalizer);
    if (post_prob[i] > max_post) {
      max_post = post_prob[i];
      cp_index = i + 1;
    }
  }

  int cp_status = 1;
  if (cp_index == 1 || cp_index == 2 || cp_index == n || cp_index == (n - 1)) {
    cp_status = 0;
  }

  return List::create(
    Rcpp::Named("post_prob") = post_prob,
    Rcpp::Named("max_post_prob") = max_post,
    Rcpp::Named("cp_index") = cp_index,
    Rcpp::Named("cp_status") = cp_status
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List bnb_singlecp_constant_cpp(NumericVector x1,
                               NumericVector x2,
                               double kappa_const) {
  const int n = x1.size();
  const NumericVector cs_x1 = bnb_prefix_sum(x1);
  const NumericVector cs_x2 = bnb_prefix_sum(x2);
  const double total_x1 = cs_x1[n - 1];
  const double total_x2 = cs_x2[n - 1];
  NumericVector log_post(n - 1);

  for (int tau = 1; tau <= n - 1; ++tau) {
    const double sx11 = cs_x1[tau - 1];
    const double sx21 = cs_x2[tau - 1];
    const double sx12 = total_x1 - sx11;
    const double sx22 = total_x2 - sx21;

    log_post[tau - 1] =
      log_segment_constant(sx11, sx21, tau, kappa_const) +
      log_segment_constant(sx12, sx22, n - tau, kappa_const);
  }

  return finalize_singlecp(log_post, n);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List bnb_singlecp_dirichlet_cpp(NumericVector x1,
                                NumericVector x2,
                                double v1,
                                double v2,
                                double v3,
                                double kappa_const) {
  const int n = x1.size();
  const NumericVector cs_x1 = bnb_prefix_sum(x1);
  const NumericVector cs_x2 = bnb_prefix_sum(x2);
  const double total_x1 = cs_x1[n - 1];
  const double total_x2 = cs_x2[n - 1];
  NumericVector log_post(n - 1);

  for (int tau = 1; tau <= n - 1; ++tau) {
    const double sx11 = cs_x1[tau - 1];
    const double sx21 = cs_x2[tau - 1];
    const double sx12 = total_x1 - sx11;
    const double sx22 = total_x2 - sx21;

    log_post[tau - 1] =
      log_segment_dirichlet_first(sx11, sx21, tau, kappa_const, v1, v2, v3) +
      log_segment_dirichlet_second(sx12, sx22, n - tau, kappa_const, v1, v2, v3);
  }

  return finalize_singlecp(log_post, n);
}
