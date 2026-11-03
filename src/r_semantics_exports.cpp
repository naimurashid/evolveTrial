// r_semantics_exports.cpp -- thin R-callable wrappers over src/r_semantics.h, used by tests/testthat/test-r-semantics.R.
// All are exported with rng = false: none of them touches R's random number stream, and they must not add an RNGScope.
#include <Rcpp.h>
#include <limits>
#include "r_semantics.h"

// [[Rcpp::export(rng = false)]]
int rsem_ld_digits() { return std::numeric_limits<long double>::digits; }

// [[Rcpp::export(rng = false)]]
double rsem_sum(Rcpp::NumericVector x) { return evt_rsem::sum_ld(REAL(x), XLENGTH(x)); }

// [[Rcpp::export(rng = false)]]
double rsem_mean(Rcpp::NumericVector x) { return evt_rsem::mean_ld(REAL(x), XLENGTH(x)); }

// [[Rcpp::export(rng = false)]]
double rsem_median(Rcpp::NumericVector x) {
  return evt_rsem::median_r(std::vector<double>(x.begin(), x.end()));
}

// [[Rcpp::export(rng = false)]]
double rsem_quantile7(Rcpp::NumericVector x, double p) {
  for (R_xlen_t i = 0; i < x.size(); i++)
    if (ISNAN(x[i])) Rcpp::stop("missing values and NaN's not allowed");
  return evt_rsem::quantile7(std::vector<double>(x.begin(), x.end()), p);
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector rsem_find_interval(Rcpp::NumericVector x, Rcpp::NumericVector breaks, bool rightmost_closed,
                                       bool all_inside, bool left_open) {
  Rcpp::IntegerVector out(x.size());
  const int nb = (int) breaks.size();
  for (R_xlen_t i = 0; i < x.size(); i++)
    out[i] = evt_rsem::find_interval(REAL(breaks), nb, x[i], rightmost_closed, all_inside, left_open);
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector rsem_mean_logical(Rcpp::IntegerVector count, Rcpp::IntegerVector n) {
  if (count.size() != n.size()) Rcpp::stop("count and n must have the same length");
  Rcpp::NumericVector out(count.size());
  for (R_xlen_t i = 0; i < count.size(); i++) out[i] = evt_rsem::mean_logical(count[i], n[i]);
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector rsem_seq_probs(int K) {
  std::vector<double> v = evt_rsem::seq_probs(K);
  return Rcpp::NumericVector(v.begin(), v.end());
}

// Returns NULL where rebalance_cutpoints_from_events returns NULL.
// [[Rcpp::export(rng = false)]]
SEXP rsem_rebalance_cutpoints(Rcpp::NumericVector event_times, double max_follow_up, int num_intervals) {
  std::vector<double> q = evt_rsem::rebalance_cutpoints(
      std::vector<double>(event_times.begin(), event_times.end()), max_follow_up, num_intervals);
  if (q.empty()) return R_NilValue;
  return Rcpp::NumericVector(q.begin(), q.end());
}
