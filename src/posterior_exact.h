// posterior_exact.h -- exact posterior probabilities for the phase 2 chunk engine (engine speed-up phase 2, step P2-1;
// plan dev/plans/2026-11-04-phase2-plan.md Section 1 item 1, report dev/perf/phase2/P2-1_REPORT.md).
//
// Plain C++ (no R objects), callable from the phase 2 chunk engine. Two quantities replace the Monte Carlo posterior
// probabilities of the phase 1 engine (src/simulate_chunk.cpp):
//
//  1. Single-arm (HC). With lambda_k ~ Gamma(shape alpha_k + E_k, rate beta_k + PT_k) independent (the posterior drawn by
//     draw_posterior_hazard_samples_cpp, src/posterior_sampling.cpp:49-63) and med the piecewise-exponential median of
//     calculate_median_survival_piecewise_impl (src/posterior_sampling.cpp:98-133; the last interval extends past the
//     last cutpoint), med < thr exactly when S = sum_k w_k lambda_k > log 2, w_k = length of interval k below thr.
//     Pr(med < thr) = Pr(S > log 2) and Pr(med > thr) = 1 - Pr(med < thr).
//  2. Between-arm (vs-reference, PH). logHR ~ N(bmean, bsd) as drawn by R::rnorm(bmean, bsd) in ph_sample_core
//     (src/simulate_chunk.cpp:448); Pr(logHR >= m) and Pr(logHR < m) by R::pnorm, with R's rnorm() edge semantics.
//
// Nothing here draws random numbers.

#ifndef EVT_POSTERIOR_EXACT_H
#define EVT_POSTERIOR_EXACT_H

#include <vector>

namespace evt_exact {

// Diagnostics of one tail evaluation (filled when a pointer is passed).
struct TailInfo {
  int method = 0;        // 0 closed by an edge rule, 1 pgamma (one distinct scale), 2 contour integral
  int side = 0;          // contour only: +1 vertex right of the pole (computes the CDF), -1 left (computes the tail)
  int n_terms = 0;       // gamma terms after dropping zero terms and merging equal scales
  int n_eval = 0;        // integrand evaluations
  int n_panels = 0;      // quadrature panels at exit
  double err_est = 0.0;  // quadrature error estimate (sum over panels of |K21 - G10|, divided by pi)
  double round_est = 0.0;// rounding estimate: DBL_EPSILON * integral of |integrand| / pi
  bool converged = true; // false if the panel budget ran out before err_est met the tolerance
};

// Pr(sum_k scale_k * G_k > c), G_k ~ Gamma(shape_k, 1) independent, shape_k >= 0, scale_k >= 0, c > 0 finite.
// Terms with shape 0 or scale 0 are identically 0 and are dropped. Returns NaN for invalid input (negative or NaN shape or
// scale, c not positive and finite).
double gamma_sum_upper_tail(const std::vector<double>& shape, const std::vector<double>& scale, double c,
                            TailInfo* info = nullptr);

// status: 0 ok; 1 ok but at least one posterior rate beta_k + PT_k was <= 0 and was replaced by 1e-6, as
// draw_posterior_hazard_samples_cpp does (it also emits an R warning; the caller decides whether to warn);
// 2 invalid input (any NaN, negative shape, NaN interval length): both probabilities NaN.
struct MedianProbs {
  double below = 0.0;  // Pr(med < thr)
  double above = 0.0;  // Pr(med > thr)
  int status = 0;
  TailInfo info;
};

// E, PT, alpha, beta, lens: K-vectors (events, person-time, prior shape, prior rate, interval lengths diff(cutpoints)).
MedianProbs median_probs(const int* E, const double* PT, const double* alpha, const double* beta, const double* lens,
                         int K, double thr);

// Pr(logHR >= m) (ge) and Pr(logHR < m) (lt) for logHR ~ R::rnorm(bmean, bsd). Mirrors R's rnorm(): NaN mean, NaN or
// infinite or negative sd give NaN draws, so both probabilities are NaN (phase 1: mean() is NA, the engine stops);
// sd == 0 or an infinite mean give the constant bmean, so the probabilities are the indicators. NaN m gives NaN.
struct LogHRProbs {
  double ge = 0.0;
  double lt = 0.0;
};
LogHRProbs loghr_probs(double bmean, double bsd, double m);

}  // namespace evt_exact

#endif  // EVT_POSTERIOR_EXACT_H
