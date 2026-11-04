// posterior_exact.cpp -- exact posterior probabilities for the phase 2 chunk engine (step P2-1). See posterior_exact.h for
// the contract and dev/perf/phase2/P2-1_REPORT.md for the derivation, the error control and the verification (gate G1).
//
// Single-arm tail Pr(S > c), S = sum_k theta_k G_k, G_k ~ Gamma(a_k, 1):
//   * no term left (all shapes or scales zero): 0; one distinct scale (terms with equal scale merged, shapes added):
//     R::pgamma upper tail.
//   * two or more distinct scales: Laplace inversion on a deformed Bromwich contour. With c scaled to 1 (tau = theta / c)
//     and L(s) = E exp(-s T) = prod_k (1 + tau_k s)^(-a_k), analytic off the cuts (-inf, -1/tau_k],
//         Pr(T <= 1) = (1 / 2 pi i) int_G exp(s) L(s) / s ds
//     for any contour G from -i inf to +i inf passing right of 0 and of every cut, and, when G passes right of every cut
//     but left of the pole at 0, the same integral equals Pr(T <= 1) - 1 = -Pr(T > 1). The vertex gamma of G is the
//     real saddle point of exp(s) L(s) / s on the side (s > 0 or -1/tau_max < s < 0) whose saddle value is smaller, so
//     the smaller of the two probabilities is computed directly. G is the hyperbola
//         s(y) = gamma + i y - beta_c (sqrt(y^2 + b^2) - b),   y in R,  beta_c = 1/2,
//     vertical at the vertex (the steepest-descent direction) and turning into the left half plane, where exp(s) decays
//     exponentially, so the integrand decays like exp(-y / 2) for every total shape (the vertical Gil-Pelaez line decays
//     only like y^-(1 + sum a), which is the hard case for small shapes). By conjugate symmetry the integral is
//     (1 / pi) Im int_0^inf exp(s) L(s) / s s'(y) dy. The upper limit Y is set by a rigorous bound on the remaining tail,
//     and [0, Y] is integrated by globally adaptive Gauss-Kronrod (10/21 points) to an absolute error estimate of 1e-12.
//
// Nothing in this file draws random numbers.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>   // before Rcpp.h
#include <Rmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <vector>
#include "posterior_exact.h"

// The one definition is src/posterior_sampling.cpp:196 (declared here as in src/simulate_chunk.cpp:33-39).
Rcpp::List ph_beta_mode_var_cpp(const arma::ivec& E_C, const arma::vec& PT_C, const arma::ivec& E_T, const arma::vec& PT_T,
                                const arma::vec& alpha_prior, const arma::vec& beta_prior, double mu, double sigma,
                                double tol, int max_iter);

namespace evt_exact {
namespace {

typedef std::complex<double> cd;

const double kBend = 0.5;          // beta_c: asymptotic slope of the contour, Re s ~ gamma - beta_c * |y|
const double kTolAbs = 1e-12;      // target absolute error of the probability (quadrature estimate)
const double kTailTol = 1e-16;     // bound on the neglected integral beyond Y, as a probability
const int kMaxPanels = 2000;

// QUADPACK qk21 abscissae and weights (Kronrod 21 points; Gauss 10 points at the odd-indexed abscissae).
const double xgk[11] = {0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
                        0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
                        0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
                        0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
                        0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
                        0.000000000000000000000000000000000};
const double wgk[11] = {0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
                        0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
                        0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
                        0.123491976262065851077548888495640, 0.134709217311473325928054001771707,
                        0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
                        0.149445554002916905664936468389821};
const double wg[5] = {0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
                      0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
                      0.295524224714752870173892994651338};

// The scaled problem: T = sum tau_k G_k compared with 1. q_k = 1 + tau_k * gamma is stored separately so that it keeps
// full relative accuracy when gamma is close to -1/tau_max.
struct Contour {
  int n;
  const double* a;
  const double* tau;
  std::vector<double> q;   // 1 + tau_k * gamma
  double gamma;
  double b;                // bend scale (local Gaussian width at the vertex)
};

struct PanelOut { double k21, g10, abs21; };

inline double integrand(const Contour& C, double y, double* absval) {
  const double r = std::sqrt(y * y + C.b * C.b);
  const double dre = -kBend * (y * y) / (r + C.b);           // -beta_c (sqrt(y^2 + b^2) - b), no cancellation
  const cd delta(dre, y);                                      // s - gamma
  const cd s = C.gamma + delta;
  double lre = 0.0, lim = 0.0;                                 // log L(s) = -sum a_k log(q_k + tau_k delta)
  for (int k = 0; k < C.n; k++) {
    const cd z(C.q[k] + C.tau[k] * delta.real(), C.tau[k] * delta.imag());
    const cd lz = std::log(z);
    lre -= C.a[k] * lz.real();
    lim -= C.a[k] * lz.imag();
  }
  const cd e = std::exp(cd(s.real() + lre, s.imag() + lim));
  const cd ds(-kBend * y / r, 1.0);
  const cd f = e / s * ds;
  if (absval) *absval = std::abs(f);
  return f.imag();
}

PanelOut gk21(const Contour& C, double lo, double hi, int* n_eval) {
  const double c = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
  double av;
  double fc = integrand(C, c, &av);
  double k21 = fc * wgk[10], g10 = 0.0, a21 = av * wgk[10];
  for (int j = 0; j < 10; j++) {
    double a1, a2;
    const double f1 = integrand(C, c - h * xgk[j], &a1);
    const double f2 = integrand(C, c + h * xgk[j], &a2);
    k21 += wgk[j] * (f1 + f2);
    a21 += wgk[j] * (a1 + a2);
    if (j % 2 == 1) g10 += wg[j / 2] * (f1 + f2);
  }
  *n_eval += 21;
  return PanelOut{k21 * h, g10 * h, a21 * h};
}

// g(s) = s - sum a log(1 + tau s) - log|s| and its first two derivatives, with q_k = 1 + tau_k s given.
struct GVal { double g, g1, g2; };
GVal gfun(int n, const double* a, const double* tau, const std::vector<double>& q, double s) {
  GVal r{s - std::log(std::fabs(s)), 1.0 - 1.0 / s, 1.0 / (s * s)};
  for (int k = 0; k < n; k++) {
    const double u = tau[k] / q[k];
    r.g -= a[k] * std::log(q[k]);
    r.g1 -= a[k] * u;
    r.g2 += a[k] * u * u;
  }
  return r;
}

// Right saddle on (1, A + 2]: g' is increasing, g'(1) <= 0 < g'(A + 2).
double saddle_right(int n, const double* a, const double* tau, std::vector<double>& q, GVal* out) {
  double A = 0.0;
  for (int k = 0; k < n; k++) A += a[k];
  double lo = 1.0, hi = A + 2.0, s = 0.5 * (lo + hi);
  GVal v{};
  for (int it = 0; it < 200; it++) {
    for (int k = 0; k < n; k++) q[k] = 1.0 + tau[k] * s;
    v = gfun(n, a, tau, q, s);
    if (v.g1 > 0) hi = s; else lo = s;
    if (hi - lo <= 1e-12 * hi) break;
    double sn = s - v.g1 / v.g2;                         // Newton step, kept inside the bracket
    if (!(sn > lo && sn < hi)) sn = 0.5 * (lo + hi);
    if (std::fabs(sn - s) <= 1e-13 * std::fabs(s)) { s = sn; break; }
    s = sn;
  }
  for (int k = 0; k < n; k++) q[k] = 1.0 + tau[k] * s;
  *out = gfun(n, a, tau, q, s);
  return s;
}

// Left saddle on (-1/tau_max, 0), parametrized by v = 1 + tau_max s in (0, 1) so that q_k = (1 - r_k) + r_k v with
// r_k = tau_k / tau_max keeps full accuracy near the branch point. g' is increasing in s (and in v).
double saddle_left(int n, const double* a, const double* tau, int kmax, std::vector<double>& q, GVal* out) {
  const double tmax = tau[kmax];
  auto setq = [&](double v) {
    for (int k = 0; k < n; k++) {
      const double rk = tau[k] / tmax;
      q[k] = (k == kmax) ? v : (1.0 - rk) + rk * v;
    }
  };
  double lo = 0.0, hi = 1.0, v = 0.5;   // root in (lo, hi); g'(v -> 0) = -inf, g'(v -> 1) = +inf
  GVal gv{};
  for (int it = 0; it < 400; it++) {
    const double s = -(1.0 - v) / tmax;
    setq(v);
    gv = gfun(n, a, tau, q, s);
    if (gv.g1 > 0) hi = v; else lo = v;
    if (hi - lo <= 1e-12 * hi) break;
    // Newton in v: dg'/dv = g'' * ds/dv = g'' / tmax
    double vn = v - gv.g1 * tmax / gv.g2;
    if (!(vn > lo && vn < hi)) {
      vn = (lo > 0 && hi / lo > 4.0) ? std::sqrt(lo * hi) : 0.5 * (lo + hi);   // geometric bisection across scales
      if (lo == 0.0) vn = (hi > 1e-300) ? hi * 1e-3 : 0.5 * hi;
    }
    if (std::fabs(vn - v) <= 1e-13 * v) { v = vn; break; }
    v = vn;
  }
  setq(v);
  const double s = -(1.0 - v) / tmax;
  *out = gfun(n, a, tau, q, s);
  return s;
}

// log of the bound on |integrand| * (remaining length) beyond Y (see the file header): for y >= Y,
// |exp(s)| <= exp(gamma + beta_c b - beta_c y), |s| >= y, |1 + tau s| >= tau y, |s'| <= sqrt(1 + beta_c^2), and
// prod max(1, (tau y)^-a) / y is non-increasing, so int_Y^inf |integrand| dy <= that / beta_c.
double log_tail_bound(const Contour& C, double Y) {
  double lb = C.gamma + kBend * C.b - kBend * Y - std::log(kBend) + 0.5 * std::log(1.0 + kBend * kBend) - std::log(Y);
  for (int k = 0; k < C.n; k++) {
    const double t = -C.a[k] * std::log(C.tau[k] * Y);
    if (t > 0) lb += t;
  }
  return lb;
}

double contour_tail(int n, const double* a, const double* tau, TailInfo* info) {
  std::vector<double> qr(n), ql(n);
  int kmax = 0;
  for (int k = 1; k < n; k++) if (tau[k] > tau[kmax]) kmax = k;
  GVal gr{}, gl{};
  const double sr = saddle_right(n, a, tau, qr, &gr);
  const double sl = saddle_left(n, a, tau, kmax, ql, &gl);
  const bool right = gr.g <= gl.g;

  Contour C;
  C.n = n;
  C.a = a;
  C.tau = tau;
  C.q = right ? qr : ql;
  C.gamma = right ? sr : sl;
  const GVal& gv = right ? gr : gl;
  C.b = 1.0 / std::sqrt(gv.g2);
  // distance from the vertex to the nearest singular point: the pole at 0, or the branch point -1/tau_max on the left
  double d = std::fabs(C.gamma);
  if (!right) d = std::min(d, C.q[kmax] / tau[kmax]);
  const double h0 = std::min(C.b, d);

  // upper limit: the smallest Y = h0 * 2^j beyond which the neglected integral is below kTailTol * pi
  const double logtol = std::log(kTailTol * M_PI);
  std::vector<double> brk;
  brk.push_back(0.0);
  double Y = h0;
  brk.push_back(Y);
  for (int j = 0; j < 2000 && log_tail_bound(C, Y) > logtol; j++) {
    Y *= 2.0;
    brk.push_back(Y);
  }

  struct Panel { double lo, hi, val, err, absv; };
  std::vector<Panel> P;
  int n_eval = 0;
  double total = 0.0, err = 0.0, absint = 0.0;
  for (size_t i = 0; i + 1 < brk.size(); i++) {
    PanelOut o = gk21(C, brk[i], brk[i + 1], &n_eval);
    P.push_back(Panel{brk[i], brk[i + 1], o.k21, std::fabs(o.k21 - o.g10), o.abs21});
  }
  bool converged = false;
  for (;;) {
    total = 0.0;
    err = 0.0;
    absint = 0.0;
    size_t worst = 0;
    for (size_t i = 0; i < P.size(); i++) {
      total += P[i].val;
      err += P[i].err;
      absint += P[i].absv;
      if (P[i].err > P[worst].err) worst = i;
    }
    if (err / M_PI <= kTolAbs) { converged = true; break; }
    if ((int) P.size() >= kMaxPanels) break;
    const Panel w = P[worst];
    const double mid = 0.5 * (w.lo + w.hi);
    PanelOut o1 = gk21(C, w.lo, mid, &n_eval);
    PanelOut o2 = gk21(C, mid, w.hi, &n_eval);
    P[worst] = Panel{w.lo, mid, o1.k21, std::fabs(o1.k21 - o1.g10), o1.abs21};
    P.push_back(Panel{mid, w.hi, o2.k21, std::fabs(o2.k21 - o2.g10), o2.abs21});
  }
  const double I = total / M_PI;     // right: Pr(T <= 1); left: -Pr(T > 1)
  double tail = right ? 1.0 - I : -I;
  if (tail < 0.0) tail = 0.0;
  if (tail > 1.0) tail = 1.0;
  if (info) {
    info->method = 2;
    info->side = right ? 1 : -1;
    info->n_eval = n_eval;
    info->n_panels = (int) P.size();
    info->err_est = err / M_PI;
    info->round_est = DBL_EPSILON * absint / M_PI;
    info->converged = converged;
  }
  return tail;
}

}  // namespace

double gamma_sum_upper_tail(const std::vector<double>& shape, const std::vector<double>& scale, double c,
                            TailInfo* info) {
  if (info) *info = TailInfo();
  if (shape.size() != scale.size() || !(c > 0) || !R_FINITE(c)) return R_NaN;
  // drop identically-zero terms, merge equal scales (Gamma(a1, t) + Gamma(a2, t) = Gamma(a1 + a2, t))
  std::vector<std::pair<double, double>> ts;   // (scale, shape)
  for (size_t k = 0; k < shape.size(); k++) {
    const double a = shape[k], t = scale[k];
    if (ISNAN(a) || ISNAN(t) || a < 0 || t < 0 || !R_FINITE(a) || !R_FINITE(t)) return R_NaN;
    if (a == 0 || t == 0) continue;
    ts.push_back(std::make_pair(t, a));
  }
  std::sort(ts.begin(), ts.end());
  std::vector<double> A, T;
  for (size_t i = 0; i < ts.size(); i++) {
    if (!T.empty() && ts[i].first == T.back()) {
      A.back() += ts[i].second;
    } else {
      T.push_back(ts[i].first);
      A.push_back(ts[i].second);
    }
  }
  const int n = (int) T.size();
  if (info) info->n_terms = n;
  if (n == 0) {                    // S = 0 < c
    if (info) info->method = 0;
    return 0.0;
  }
  if (n == 1) {
    if (info) info->method = 1;
    return R::pgamma(c, A[0], T[0], /*lower_tail=*/0, /*log_p=*/0);
  }
  for (int k = 0; k < n; k++) T[k] /= c;   // scale c to 1
  return contour_tail(n, A.data(), T.data(), info);
}

MedianProbs median_probs(const int* E, const double* PT, const double* alpha, const double* beta, const double* lens,
                         int K, double thr) {
  MedianProbs out;
  if (K < 1) {
    out.below = out.above = R_NaN;
    out.status = 2;
    return out;
  }
  // posterior parameters exactly as draw_posterior_hazard_samples_cpp (src/posterior_sampling.cpp:49-58)
  std::vector<double> shape(K), rate(K);
  for (int k = 0; k < K; k++) {
    if (E[k] == NA_INTEGER || !R_FINITE(PT[k]) || !R_FINITE(alpha[k]) || !R_FINITE(beta[k]) || !R_FINITE(lens[k]) ||
        lens[k] < 0) {
      out.status = 2;
      break;
    }
    shape[k] = alpha[k] + E[k];
    if (shape[k] < 0) {             // R::rgamma gives NaN, every median is NaN
      out.status = 2;
      break;
    }
    double r = beta[k] + PT[k];
    if (r <= 0) {
      r = 1e-6;
      out.status = 1;
    }
    rate[k] = r;
  }
  if (out.status == 2) {
    out.below = out.above = R_NaN;
    return out;
  }
  // NaN threshold: every comparison med < thr and med > thr is FALSE, so the phase 1 tally is 0 for both.
  if (ISNAN(thr)) {
    out.below = out.above = 0.0;
    out.info.method = 0;
    return out;
  }
  // weights: the part of interval k = [start_k, start_k + lens_k) below thr; the last interval is open-ended.
  // start_k by repeated addition, as current_time += delta_t_j in calculate_median_survival_piecewise_impl.
  std::vector<double> sh, sc;
  double start = 0.0;
  bool infinite = false;
  for (int k = 0; k < K; k++) {
    const double end = start + lens[k];
    double w = 0.0;
    if (thr > start) {
      if (k < K - 1) w = ((thr < end) ? thr : end) - start;
      else w = thr - start;          // +inf when thr = +inf
    }
    if (w > 0 && shape[k] > 0) {
      if (!R_FINITE(w)) {
        infinite = true;             // last hazard > 0 almost surely, cumulative hazard at +inf is infinite
      } else {
        sh.push_back(shape[k]);
        sc.push_back(w / rate[k]);   // w_k lambda_k = (w_k / rate_k) Gamma(shape_k, 1)
      }
    }
    start = end;
  }
  double below;
  if (infinite) {
    below = 1.0;
    out.info.method = 0;
  } else {
    below = gamma_sum_upper_tail(sh, sc, M_LN2, &out.info);
  }
  out.below = below;
  out.above = 1.0 - below;
  return out;
}

LogHRProbs loghr_probs(double bmean, double bsd, double m) {
  LogHRProbs out;
  // R 4.6.0 src/nmath/rnorm.c: if (ISNAN(mu) || !R_FINITE(sigma) || sigma < 0.) ML_WARN_return_NAN;
  //                            if (sigma == 0. || !R_FINITE(mu)) return mu; else mu + sigma * norm_rand()
  if (ISNAN(bmean) || !R_FINITE(bsd) || bsd < 0 || ISNAN(m)) {
    out.ge = out.lt = R_NaN;
    return out;
  }
  if (bsd == 0 || !R_FINITE(bmean)) {
    out.ge = (bmean >= m) ? 1.0 : 0.0;
    out.lt = (bmean < m) ? 1.0 : 0.0;
    return out;
  }
  out.ge = R::pnorm(m, bmean, bsd, /*lower_tail=*/0, /*log_p=*/0);   // Pr(logHR > m) = Pr(logHR >= m)
  out.lt = R::pnorm(m, bmean, bsd, /*lower_tail=*/1, /*log_p=*/0);
  return out;
}

}  // namespace evt_exact

// ---- R test hooks ---------------------------------------------------------------------------------------------------

namespace {
Rcpp::List info_list(const evt_exact::TailInfo& i) {
  return Rcpp::List::create(Rcpp::Named("method") = i.method, Rcpp::Named("side") = i.side,
                            Rcpp::Named("n_terms") = i.n_terms, Rcpp::Named("n_eval") = i.n_eval,
                            Rcpp::Named("n_panels") = i.n_panels, Rcpp::Named("err_est") = i.err_est,
                            Rcpp::Named("round_est") = i.round_est, Rcpp::Named("converged") = i.converged);
}
void check_len(int K, Rcpp::IntegerVector E, Rcpp::NumericVector PT, Rcpp::NumericVector alpha,
               Rcpp::NumericVector beta, Rcpp::NumericVector lens) {
  if (E.size() != K || PT.size() != K || alpha.size() != K || beta.size() != K || lens.size() != K || K < 1) {
    Rcpp::stop("E, PT, alpha, beta and lens must have the same positive length");
  }
}
}  // namespace

// [[Rcpp::export(name = ".exact_gamma_sum_tail", rng = false)]]
Rcpp::List exact_gamma_sum_tail(Rcpp::NumericVector shape, Rcpp::NumericVector scale, double c) {
  // Test hook: Pr(sum scale_k G_k > c), G_k ~ Gamma(shape_k, 1), with diagnostics.
  evt_exact::TailInfo info;
  const double p = evt_exact::gamma_sum_upper_tail(std::vector<double>(shape.begin(), shape.end()),
                                                   std::vector<double>(scale.begin(), scale.end()), c, &info);
  Rcpp::List out = info_list(info);
  out.push_front(p, "tail");
  return out;
}

// [[Rcpp::export(name = ".exact_pr_median", rng = false)]]
Rcpp::List exact_pr_median(Rcpp::IntegerVector E, Rcpp::NumericVector PT, Rcpp::NumericVector alpha,
                           Rcpp::NumericVector beta, Rcpp::NumericVector lens, double thr) {
  // Test hook: Pr(med < thr), Pr(med > thr) and diagnostics for one threshold.
  const int K = E.size();
  check_len(K, E, PT, alpha, beta, lens);
  evt_exact::MedianProbs m = evt_exact::median_probs(E.begin(), PT.begin(), alpha.begin(), beta.begin(), lens.begin(), K,
                                                     thr);
  Rcpp::List out = info_list(m.info);
  out.push_front(m.status, "status");
  out.push_front(m.above, "above");
  out.push_front(m.below, "below");
  return out;
}

// [[Rcpp::export(name = ".exact_pr_median_below", rng = false)]]
Rcpp::NumericVector exact_pr_median_below(Rcpp::IntegerVector E, Rcpp::NumericVector PT, Rcpp::NumericVector alpha,
                                          Rcpp::NumericVector beta, Rcpp::NumericVector lens, Rcpp::NumericVector thr) {
  // Test hook: Pr(med < thr) for each element of thr (the exact counterpart of mean(med < thr) in interim_hc and the
  // final HC futility probability).
  const int K = E.size();
  check_len(K, E, PT, alpha, beta, lens);
  Rcpp::NumericVector out(thr.size());
  for (R_xlen_t i = 0; i < thr.size(); i++) {
    out[i] = evt_exact::median_probs(E.begin(), PT.begin(), alpha.begin(), beta.begin(), lens.begin(), K, thr[i]).below;
  }
  return out;
}

// [[Rcpp::export(name = ".exact_pr_median_above", rng = false)]]
Rcpp::NumericVector exact_pr_median_above(Rcpp::IntegerVector E, Rcpp::NumericVector PT, Rcpp::NumericVector alpha,
                                          Rcpp::NumericVector beta, Rcpp::NumericVector lens, Rcpp::NumericVector thr) {
  // Test hook: Pr(med > thr) (the exact counterpart of mean(med > success_thr) in the final HC analysis).
  const int K = E.size();
  check_len(K, E, PT, alpha, beta, lens);
  Rcpp::NumericVector out(thr.size());
  for (R_xlen_t i = 0; i < thr.size(); i++) {
    out[i] = evt_exact::median_probs(E.begin(), PT.begin(), alpha.begin(), beta.begin(), lens.begin(), K, thr[i]).above;
  }
  return out;
}

// [[Rcpp::export(name = ".exact_loghr_probs", rng = false)]]
Rcpp::NumericMatrix exact_loghr_probs(double bmean, double bsd, Rcpp::NumericVector m) {
  // Test hook: columns ge = Pr(logHR >= m), lt = Pr(logHR < m) for logHR ~ rnorm(bmean, bsd), one row per m.
  Rcpp::NumericMatrix out(m.size(), 2);
  for (R_xlen_t i = 0; i < m.size(); i++) {
    evt_exact::LogHRProbs p = evt_exact::loghr_probs(bmean, bsd, m[i]);
    out(i, 0) = p.ge;
    out(i, 1) = p.lt;
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("ge", "lt");
  return out;
}

// [[Rcpp::export(name = ".exact_pr_loghr", rng = false)]]
Rcpp::List exact_pr_loghr(Rcpp::IntegerVector E_C, Rcpp::NumericVector PT_C, Rcpp::IntegerVector E_T,
                          Rcpp::NumericVector PT_T, Rcpp::NumericVector alpha, Rcpp::NumericVector beta,
                          Rcpp::NumericVector interval_lengths, double mu, double sigma, double tol, int max_iter,
                          Rcpp::NumericVector margins) {
  // Test hook with the inputs of .chunk_ph_sample (minus num_samples, plus the log-HR thresholds m): the Laplace mean and
  // sd from ph_beta_mode_var_cpp, called exactly as ph_sample_core does (src/simulate_chunk.cpp:432-445), and the
  // closed-form Pr(logHR >= m), Pr(logHR < m) for each m. interval_lengths is accepted for signature parity and is not
  // read (the log-HR posterior does not depend on it).
  const int K = alpha.size();
  if (E_C.size() != K || PT_C.size() != K || E_T.size() != K || PT_T.size() != K || beta.size() != K ||
      interval_lengths.size() != K) {
    Rcpp::stop(".exact_pr_loghr: all vectors must have the same length");
  }
  arma::ivec aEC(K), aET(K);
  arma::vec aPTC(K), aPTT(K), aA(K), aB(K);
  for (int k = 0; k < K; k++) {
    aEC(k) = E_C[k];
    aET(k) = E_T[k];
    aPTC(k) = PT_C[k];
    aPTT(k) = PT_T[k];
    aA(k) = alpha[k];
    aB(k) = beta[k];
  }
  Rcpp::List bp = ph_beta_mode_var_cpp(aEC, aPTC, aET, aPTT, aA, aB, mu, sigma, tol, max_iter);
  const double bmean = Rcpp::as<double>(bp["mean"]);
  const double bsd = Rcpp::as<double>(bp["sd"]);
  Rcpp::NumericVector ge(margins.size()), lt(margins.size());
  for (R_xlen_t i = 0; i < margins.size(); i++) {
    evt_exact::LogHRProbs p = evt_exact::loghr_probs(bmean, bsd, margins[i]);
    ge[i] = p.ge;
    lt[i] = p.lt;
  }
  return Rcpp::List::create(Rcpp::Named("mean") = bmean, Rcpp::Named("sd") = bsd, Rcpp::Named("pr_ge") = ge,
                            Rcpp::Named("pr_lt") = lt);
}
