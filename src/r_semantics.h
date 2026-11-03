// r_semantics.h -- floating-point helpers that reproduce R's own arithmetic bit for bit (engine C++ port, step 2b).
//
// Each helper is a transcription of the R 4.6.0 source named in the comment above it (R-4.6.0.tar.gz from CRAN; the
// copies read for this file are under dev/perf/build/rsrc/R-4.6.0/src). The order of operations and the types of the
// intermediates (long double where R uses LDOUBLE) are kept exactly. Do not "tidy" an expression: a reordering that is
// mathematically equal changes the last bit.
//
// Scope (plan 2b.6): sum, mean, median, type 7 quantile, findInterval, mean of a logical vector, seq(0, 1, length.out),
// and the cutpoint rebalancing function. Slice metrics and the posterior are port work, not helpers.
//
// Platform: bit identity with R relies on `long double` being the 80-bit x87 type (x86-64 Linux), and on the build not
// contracting a * b + c into an FMA (src/Makevars: -ffp-contract=off). rsem_ld_digits() lets a test assert the former.
#ifndef EVT_R_SEMANTICS_H
#define EVT_R_SEMANTICS_H

#include <R.h>              // R_NaReal, R_FINITE, R_xlen_t, Rboolean via R_ext
#include <R_ext/Arith.h>
#include <R_ext/Boolean.h>
#include <R_ext/Utils.h>    // findInterval2
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace evt_rsem {

// R: src/main/summary.c, rsum() (called from do_summary() for sum(x) of a double vector, na.rm = FALSE).
//   LDOUBLE s = 0.0; for each k: s += x[k]; then Inf / -Inf if beyond +-DBL_MAX, else (double) s.
// do_summary() then forms zcum.r += tmp with zcum.r a double that starts at 0, so a result of -0 becomes +0.
inline double sum_ld(const double* x, R_xlen_t n) {
  long double s = 0.0;
  for (R_xlen_t k = 0; k < n; k++) s += x[k];
  double value;
  if (s > DBL_MAX) value = R_PosInf;
  else if (s < -DBL_MAX) value = R_NegInf;
  else value = (double) s;
  double zcum = 0.0;
  zcum += value;
  return zcum;
}

// R: src/main/summary.c, real_mean() (called from do_summary() for mean(x) of a double vector, i.e. mean.default).
// A long double sum divided by n, then (when both the sum and the mean are finite) a second pass that adds the mean
// of the residuals. The two fall-back branches for an overflowing sum are kept as written.
inline double mean_ld(const double* x, R_xlen_t n) {
  long double s = 0.0;
  for (R_xlen_t k = 0; k < n; k++) s += x[k];
  bool finite_s = R_FINITE((double) s) != 0;
  if (finite_s) {
    s /= n;
  } else {  // infinite s, maybe just overflowed; try to use smaller terms
    s = 0.;
    for (R_xlen_t k = 0; k < n; k++) s += x[k] / n;
  }
  if (finite_s && R_FINITE((double) s)) {
    long double t = 0.0;
    for (R_xlen_t k = 0; k < n; k++) t += (x[k] - s);
    s += t / n;
  } else if (R_FINITE((double) s)) {  // was infinite: more careful
    long double t = 0.0;
    for (R_xlen_t k = 0; k < n; k++) t += (x[k] - s) / n;
    s += t;
  }
  return (double) s;
}

// R: src/main/summary.c, logical_mean() (mean() of a LGLSXP): LDOUBLE s accumulates the 0/1 values (exact), the result
// is (double) (s / n) with the division in long double. `count` is the number of TRUE values among n non-NA values.
inline double mean_logical(int count, R_xlen_t n) {
  long double s = 0.0;
  s += count;
  return (double) (s / n);
}

// R: stats:::median.default (src/library/stats/R/median.R), na.rm = FALSE.
//   any NA -> NA; n == 0 -> NA; half <- (n + 1L) %/% 2L;
//   odd n: sort(x, partial = half)[half]; even n: mean(sort(x, partial = half + 0L:1L)[half + 0L:1L]).
// The mean is mean.default on two doubles, i.e. mean_ld above. Takes the vector by value because it is reordered.
inline double median_r(std::vector<double> x) {
  const R_xlen_t n = (R_xlen_t) x.size();
  for (double v : x) if (std::isnan(v)) return NA_REAL;
  if (n == 0) return NA_REAL;
  std::sort(x.begin(), x.end());
  const R_xlen_t half = (n + 1) / 2;       // 1-based
  if (n % 2 == 1) return x[half - 1];
  double two[2] = { x[half - 1], x[half] };
  return mean_ld(two, 2);
}

// R: stats:::quantile.default (src/library/stats/R/quantile.R), type = 7, fuzz = 0, one probability, no NA in x.
//   probs <- pmax(0, pmin(1, probs)); index <- 1 + max(n - 1, 0) * probs; lo <- floor(index * (1 + fuzz));
//   hi <- ceiling(index); qs <- x[lo] (x sorted); i <- index > lo & x[hi] != qs;
//   h <- (index - lo)[i]; qs[i] <- (1 - h) * qs[i] + h * x[hi[i]].
// `xs` must be sorted ascending. n == 0 gives NA, as x[1] of an empty vector does.
inline double quantile7_sorted(const std::vector<double>& xs, double p) {
  const double n = (double) xs.size();
  if (xs.empty()) return NA_REAL;
  p = std::max(0.0, std::min(1.0, p));
  const double index = 1 + std::max(n - 1, 0.0) * p;
  const double lo = std::floor(index * (1 + 0.0));
  const double hi = std::ceil(index);
  double qs = xs[(size_t) lo - 1];
  const double xhi = xs[(size_t) hi - 1];
  if (index > lo && xhi != qs) {
    const double h = index - lo;
    qs = (1 - h) * qs + h * xhi;
  }
  return qs;
}
inline double quantile7(std::vector<double> x, double p) {
  std::sort(x.begin(), x.end());
  return quantile7_sorted(x, p);
}

// R: src/appl/interv.c, findInterval2() -- the routine behind findInterval() (src/main/util.c, do_findinterval()).
// R's own routine is called rather than copied (it is declared in R_ext/Utils.h), so the result is R's by construction.
// do_findinterval() starts with ilo = 1 and passes the previous result as the hint; the hint never changes the result
// for a sorted `breaks`, so ilo = 1 is used for every point. NaN gives NA_INTEGER (do_findinterval, checkNA = TRUE).
// The engine uses rightmost_closed = FALSE, all_inside = FALSE, left_open = FALSE (state_management.R:27-31).
inline int find_interval(const double* breaks, int nb, double x, bool rightmost_closed, bool all_inside, bool left_open) {
  if (std::isnan(x)) return NA_INTEGER;
  int mflag;
  return findInterval2(const_cast<double*>(breaks), nb, x, (Rboolean) rightmost_closed, (Rboolean) all_inside,
                       (Rboolean) left_open, 1, &mflag);
}

// R: base:::seq.default (src/library/base/R/seq.R), the branch for seq(from = 0, to = 1, length.out = K + 1):
//   length.out > 2 and from != to:  c(from, from + seq_len(length.out - 2L) * (del / n1), to)  with del = to - from, n1 = length.out - 1L
//   length.out == 2:                c(from, to)
//   length.out == 1:                from
// so the interior points are 0 + i * (1 / K) (a product with a precomputed reciprocal), not i / K.
inline std::vector<double> seq_probs(int K) {
  if (K < 0) throw std::invalid_argument("seq_probs: K must be non-negative");
  const double from = 0.0, to = 1.0;
  if (K == 0) return std::vector<double>(1, from);
  if (K == 1) return std::vector<double>{ from, to };
  const double del = to - from;
  const double n1 = (double) K;   // length.out - 1L
  std::vector<double> out(K + 1);
  out[0] = from;
  for (int i = 1; i <= K - 1; i++) out[i] = from + (double) i * (del / n1);
  out[K] = to;
  return out;
}

// R: evolveTrial:::rebalance_cutpoints_from_events (R/simulation_driver.R:1331-1398), built from the helpers above.
//   if (length(event_times) < max(2, K)) NULL;  event_times <- pmin(pmax(event_times, 0), max_follow_up);
//   probs <- seq(0, 1, length.out = K + 1);  q <- quantile(event_times, probs, na.rm = TRUE, type = 7);
//   q[1] <- 0; q[K + 1] <- max_follow_up;
//   for i in 2..(K + 1): if (!is.finite(q[i]) || q[i] <= q[i - 1]) q[i] <- min(max_follow_up, q[i - 1] + 1e-6);
//   if (any(diff(q) <= 0)) NULL else q.
// Returns an empty vector where R returns NULL. After the loop a NaN difference would make R's `if (any(...))` fail with
// an error (NA condition); that cannot happen with finite max_follow_up, and is reported as an exception here.
inline std::vector<double> rebalance_cutpoints(const std::vector<double>& event_times, double max_follow_up, int K) {
  const double need = std::max(2.0, (double) K);
  if ((double) event_times.size() < need) return std::vector<double>();
  std::vector<double> x;
  x.reserve(event_times.size());
  for (double v : event_times) {
    if (std::isnan(v)) continue;                       // na.rm = TRUE (pmax/pmin keep NA, quantile drops it)
    double w = v;
    if (0.0 > w) w = 0.0;                              // pmax(v, 0)        (summary.c do_pmin: replace when tmp > ra)
    if (max_follow_up < w) w = max_follow_up;          // pmin(w, max_follow_up) (replace when tmp < ra)
    x.push_back(w);
  }
  std::sort(x.begin(), x.end());
  const std::vector<double> probs = seq_probs(K);
  std::vector<double> q(probs.size());
  for (size_t i = 0; i < probs.size(); i++) q[i] = quantile7_sorted(x, probs[i]);
  q[0] = 0.0;
  q[q.size() - 1] = max_follow_up;
  for (size_t i = 1; i < q.size(); i++) {
    if (!R_FINITE(q[i]) || q[i] <= q[i - 1]) {
      const double b = q[i - 1] + 1e-6;
      q[i] = (max_follow_up < b) ? max_follow_up : b;  // min(max_follow_up, q[i - 1] + 1e-6)
    }
  }
  for (size_t i = 0; i + 1 < q.size(); i++) {
    const double d = q[i + 1] - q[i];
    if (std::isnan(d)) throw std::runtime_error("rebalance_cutpoints: NaN cutpoint difference");
    if (d <= 0) return std::vector<double>();
  }
  return q;
}

}  // namespace evt_rsem

#endif  // EVT_R_SEMANTICS_H
