// simulate_chunk_exact.cpp -- phase 2 chunk engine (engine speed-up phase 2, step P2-2; plan
// dev/plans/2026-11-04-phase2-plan.md Section 1, report dev/perf/phase2/P2-2_REPORT.md).
//
// The trial logic is the phase 1 engine's (src/simulate_chunk.cpp, the bit-identical port of simulate_chunk_r) copied
// line for line: arrival loop, the looks_stop_at_final phantom-arrival rules, caps, gates, calendar-beat look schedule,
// rebalanced cutpoints, look slicing, post-enrollment looks, final analysis, NA and error behavior. src/simulate_chunk.cpp
// itself is not touched. The loop is written once, as a template over two policies:
//
//   Rng  -- where patient-level random numbers come from;
//   Post -- how a posterior probability is obtained.
//
// Two instantiations are compiled:
//
//   .simulate_chunk_exact(sim_indices, cfg, scenario_seed): CounterRng x ExactPost, the phase 2 engine. Every posterior
//     probability comes from the exact functions of src/posterior_exact.cpp (P2-1), nothing is drawn for posteriors.
//     Patient-level draws come from a per-trial xoshiro256++ generator seeded by (scenario_seed, trial index) through
//     splitmix64, so a trial's result depends only on the scenario seed and its index. R's generator is not used.
//   .simulate_chunk_exact_replay(sim_indices, cfg): RStreamRng x McPost, a test hook. It draws from R's stream in the
//     phase 1 order and uses the phase 1 Monte Carlo posteriors, so on the same seed it must return exactly what
//     .simulate_chunk_cpp returns. That test (tests/testthat/test-simulate-chunk-exact.R) checks that the copied trial
//     logic is the phase 1 logic; only the two policies differ between the replay and the exact engine.
//
// Supported set: exactly the phase 1 set (the same packed cfg, R/simulate_chunk_cpp.R .pack_chunk_cfg). unpack() stops
// outside it; there is no fallback.
//
// Random numbers in the exact engine (CounterRng):
//   key   = mix64(seed64 ^ mix64(trial_index + 0x9E3779B97F4A7C15)), seed64 = (uint32 seed[1] << 32) | uint32 seed[2]
//   state = four successive splitmix64 outputs from state key (the seeding recommended for xoshiro256++)
//   U     = ((x >> 12) + 0.5) * 2^-52 from one xoshiro256++ output x: in [2^-53, 1 - 2^-53], never 0 or 1
//   arrival gap      -log(U) / accrual_rate                (rexp(1, rate), inversion)
//   event time       scale * (-log(U))^(1 / shape)         (R's rweibull(shape, scale) formula)
//   arm assignment   phase 1's ProbSampleNoReplace port with rT = U, same probabilities; no draw when one arm is eligible
//   dropout          nothing drawn (random_dropout is FALSE in the supported set; censoring time Inf as in phase 1)
// Per arrival the order is: gap, (look), arm, event time.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Random.h>   // unif_rand
#include <R_ext/Utils.h>    // revsort
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include "r_semantics.h"
#include "posterior_exact.h"

// Existing bodies in src/posterior_sampling.cpp (declared as in src/simulate_chunk.cpp).
arma::mat draw_posterior_hazard_samples_cpp(int num_intervals, const arma::ivec& events_per_interval,
                                            const arma::vec& person_time_per_interval,
                                            const arma::vec& prior_alpha_params, const arma::vec& prior_beta_params,
                                            int num_samples);
arma::vec calculate_median_survival_matrix_cpp(const arma::mat& hazard_samples, const arma::vec& interval_lengths);
Rcpp::List ph_beta_mode_var_cpp(const arma::ivec& E_C, const arma::vec& PT_C, const arma::ivec& E_T, const arma::vec& PT_T,
                                const arma::vec& alpha_prior, const arma::vec& beta_prior, double mu, double sigma,
                                double tol, int max_iter);

namespace {

// ---- counter-based generator ------------------------------------------------------------------------------------------
inline uint64_t mix64(uint64_t z) {   // splitmix64 finalizer (a bijection on 64-bit words)
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
inline uint64_t splitmix_next(uint64_t& x) {
  x += 0x9E3779B97F4A7C15ULL;
  return mix64(x);
}

struct Xoshiro256pp {
  uint64_t s[4];
  static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  inline uint64_t next() {
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
  }
  void seed(uint64_t seed64, uint64_t trial) {
    uint64_t x = mix64(seed64 ^ mix64(trial + 0x9E3779B97F4A7C15ULL));
    for (int i = 0; i < 4; i++) s[i] = splitmix_next(x);
    if ((s[0] | s[1] | s[2] | s[3]) == 0) s[0] = 1;   // the all-zero state is the one fixed point; never reached in practice
  }
  // open interval (0, 1): 52 random bits, centred in their cell
  inline double unif() { return ((double) (next() >> 12) + 0.5) * 0x1.0p-52; }
};

// ---- sample(x, size = 1, prob = p): phase 1's port of ProbSampleNoReplace (src/simulate_chunk.cpp:48-76) ----------------
// with the uniform supplied by draw(). draw_single = false skips the draw when only one entry exists (the result is then
// entry 0 whatever the uniform, so the probabilities are unchanged).
template <class Draw>
int prob_sample(std::vector<double> p, Draw&& draw, bool draw_single) {
  const int n = (int) p.size();
  if (n < 1) Rcpp::stop("invalid first argument");
  double sum = 0.0;
  int npos = 0;
  for (int i = 0; i < n; i++) {
    if (!R_FINITE(p[i])) Rcpp::stop("NA in probability vector");
    if (p[i] < 0.0) Rcpp::stop("negative probability");
    if (p[i] > 0.0) {
      npos++;
      sum += p[i];
    }
  }
  if (npos == 0 || 1 > npos) Rcpp::stop("too few positive probabilities");
  if (n == 1 && !draw_single) return 0;
  for (int i = 0; i < n; i++) p[i] /= sum;
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++) perm[i] = i + 1;
  revsort(p.data(), perm.data(), n);
  double totalmass = 1;
  int n1 = n - 1;
  double rT = totalmass * draw();
  double mass = 0;
  int j;
  for (j = 0; j < n1; j++) {
    mass += p[j];
    if (rT <= mass) break;
  }
  return perm[j] - 1;
}

inline double r_pmax2(double a, double b) { return (b > a) ? b : a; }
inline double r_pmin2(double a, double b) { return (b < a) ? b : a; }

// ---- configuration (copy of src/simulate_chunk.cpp:84-221; same packed list, same checks) --------------------------------
struct Cfg {
  int n_arms;
  Rcpp::CharacterVector arm_names;
  std::vector<double> mtpa;
  std::vector<int> capacity;
  std::vector<double> rand_probs;
  std::vector<double> shape, scale;
  double accrual_rate, beat, max_trial_time, min_follow_up_at_final;
  double max_follow_up_sim, censor_max_time_sim;
  bool random_dropout;
  int K;
  std::vector<double> cut_base, len_base;
  double rebalance_threshold;
  std::vector<double> prior_alpha, prior_beta;
  int draws_interim, draws_final;
  std::vector<double> min_pat, min_events, min_median_fu, min_pt_frac, max_pt_arm;
  double futility_prob_interim;
  std::vector<double> futility_thr_interim, success_thr, futility_thr;
  double final_success_prob, final_futility_prob;
  bool collect_raw_n, return_variance;
  bool vsref = false;
  int ref_idx = -1, trt_idx = -1;
  std::vector<double> vr_min_ev, vr_min_mfu, vr_min_pt, vr_max_pt;
  int draws_vsref_interim = 0;
  double ph_mu = 0, ph_sigma = 1, ph_tol = 0;
  int ph_max_iter = 0;
  double log_margin_interim = 0, log_eff_interim = 0, log_margin_final = 0, log_eff_final = 0;
  double vr_fut_prob_interim = 0, vr_eff_prob_final = 0, vr_fut_prob_final = 0;
};

double num1(const Rcpp::List& cfg, const char* name) {
  Rcpp::NumericVector v = cfg[name];
  if (v.size() != 1) Rcpp::stop(std::string("simulate_chunk_exact: cfg$") + name + " must have length 1");
  return v[0];
}
int int1(const Rcpp::List& cfg, const char* name) {
  Rcpp::IntegerVector v = cfg[name];
  if (v.size() != 1 || v[0] == NA_INTEGER) Rcpp::stop(std::string("simulate_chunk_exact: cfg$") + name + " must be one integer");
  return v[0];
}
bool lgl1(const Rcpp::List& cfg, const char* name) {
  Rcpp::LogicalVector v = cfg[name];
  if (v.size() != 1 || v[0] == NA_LOGICAL) Rcpp::stop(std::string("simulate_chunk_exact: cfg$") + name + " must be TRUE or FALSE");
  return v[0] != 0;
}
std::vector<double> numv(const Rcpp::List& cfg, const char* name, int len) {
  Rcpp::NumericVector v = cfg[name];
  if (len >= 0 && v.size() != len) Rcpp::stop(std::string("simulate_chunk_exact: cfg$") + name + " has the wrong length");
  return std::vector<double>(v.begin(), v.end());
}

Cfg unpack(const Rcpp::List& cfg) {
  if (!cfg.containsElementNamed("supported") || !Rcpp::as<bool>(cfg["supported"])) {
    Rcpp::stop("simulate_chunk_exact: the packed configuration is outside the supported set (no fallback here)");
  }
  Cfg c;
  c.arm_names = cfg["arm_names"];
  c.n_arms = c.arm_names.size();
  const std::string mode = cfg.containsElementNamed("mode") ? Rcpp::as<std::string>(cfg["mode"]) : std::string("hc");
  if (mode == "vsref") {
    c.vsref = true;
    if (c.n_arms != 2) Rcpp::stop("simulate_chunk_exact: the vs-reference configuration needs exactly two arms");
  } else if (mode == "hc") {
    if (c.n_arms != 1) Rcpp::stop("simulate_chunk_exact: the single-arm configuration needs exactly one arm");
  } else {
    Rcpp::stop("simulate_chunk_exact: unknown cfg$mode");
  }
  c.mtpa = numv(cfg, "mtpa", c.n_arms);
  Rcpp::IntegerVector cap = cfg["capacity"];
  if (cap.size() != c.n_arms) Rcpp::stop("simulate_chunk_exact: cfg$capacity has the wrong length");
  c.capacity.assign(cap.begin(), cap.end());
  c.rand_probs = numv(cfg, "rand_probs", c.n_arms);
  c.shape = numv(cfg, "shape", c.n_arms);
  c.scale = numv(cfg, "scale", c.n_arms);
  c.accrual_rate = num1(cfg, "accrual_rate");
  c.beat = num1(cfg, "beat");
  c.max_trial_time = num1(cfg, "max_trial_time");
  c.min_follow_up_at_final = num1(cfg, "min_follow_up_at_final");
  c.max_follow_up_sim = num1(cfg, "max_follow_up_sim");
  c.censor_max_time_sim = num1(cfg, "censor_max_time_sim");
  c.random_dropout = lgl1(cfg, "random_dropout");
  if (c.random_dropout) Rcpp::stop("simulate_chunk_exact: random_dropout = TRUE is outside the supported set");
  if (!lgl1(cfg, "looks_stop_at_final")) Rcpp::stop("simulate_chunk_exact: looks_stop_at_final = FALSE is outside the supported set");
  if (!R_FINITE(c.max_trial_time)) Rcpp::stop("simulate_chunk_exact: infinite max_trial_time is outside the supported set");
  c.K = int1(cfg, "K");
  c.cut_base = numv(cfg, "cut_base", c.K + 1);
  c.len_base = numv(cfg, "len_base", c.K);
  c.rebalance_threshold = num1(cfg, "rebalance_threshold");
  c.prior_alpha = numv(cfg, "prior_alpha", c.K);
  c.prior_beta = numv(cfg, "prior_beta", c.K);
  c.draws_final = int1(cfg, "draws_final");
  c.collect_raw_n = lgl1(cfg, "collect_raw_n");
  c.return_variance = lgl1(cfg, "return_variance");
  if (c.vsref) {
    c.ref_idx = int1(cfg, "ref_idx");
    c.trt_idx = int1(cfg, "trt_idx");
    if (c.ref_idx < 0 || c.ref_idx > 1 || c.trt_idx < 0 || c.trt_idx > 1 || c.ref_idx == c.trt_idx) {
      Rcpp::stop("simulate_chunk_exact: cfg$ref_idx and cfg$trt_idx must be the two distinct arm positions");
    }
    c.vr_min_ev = numv(cfg, "vr_min_ev", 2);
    c.vr_min_mfu = numv(cfg, "vr_min_mfu", 2);
    c.vr_min_pt = numv(cfg, "vr_min_pt", 2);
    c.vr_max_pt = numv(cfg, "vr_max_pt", 2);
    c.draws_vsref_interim = int1(cfg, "draws_vsref_interim");
    if (c.draws_vsref_interim < 1 || c.draws_final < 1) Rcpp::stop("simulate_chunk_exact: posterior draw counts must be >= 1");
    c.ph_mu = num1(cfg, "ph_mu");
    c.ph_sigma = num1(cfg, "ph_sigma");
    c.ph_tol = num1(cfg, "ph_tol");
    c.ph_max_iter = int1(cfg, "ph_max_iter");
    c.log_margin_interim = num1(cfg, "log_margin_interim");
    c.log_eff_interim = num1(cfg, "log_eff_interim");
    c.log_margin_final = num1(cfg, "log_margin_final");
    c.log_eff_final = num1(cfg, "log_eff_final");
    c.vr_fut_prob_interim = num1(cfg, "vr_fut_prob_interim");
    c.vr_eff_prob_final = num1(cfg, "vr_eff_prob_final");
    c.vr_fut_prob_final = num1(cfg, "vr_fut_prob_final");
    return c;
  }
  c.draws_interim = int1(cfg, "draws_interim");
  if (c.draws_interim < 1 || c.draws_final < 1) Rcpp::stop("simulate_chunk_exact: posterior draw counts must be >= 1");
  c.min_pat = numv(cfg, "min_pat", c.n_arms);
  c.min_events = numv(cfg, "min_events", c.n_arms);
  c.min_median_fu = numv(cfg, "min_median_fu", c.n_arms);
  c.min_pt_frac = numv(cfg, "min_pt_frac", c.n_arms);
  c.max_pt_arm = numv(cfg, "max_pt_arm", c.n_arms);
  c.futility_prob_interim = num1(cfg, "futility_prob_interim");
  c.futility_thr_interim = numv(cfg, "futility_thr_interim", c.n_arms);
  c.success_thr = numv(cfg, "success_thr", c.n_arms);
  c.futility_thr = numv(cfg, "futility_thr", c.n_arms);
  c.final_success_prob = num1(cfg, "final_success_prob");
  c.final_futility_prob = num1(cfg, "final_futility_prob");
  return c;
}

// ---- per-trial state (copy of src/simulate_chunk.cpp:224-249) ------------------------------------------------------------
enum Status { RECRUITING = 0, STOPPED_EFFICACY = 1, STOPPED_FUTILITY = 2 };

struct Registry {
  std::vector<double> enroll, te, rc;
  int n = 0;
};

struct TrialState {
  std::vector<int> status, enrolled, stop_eff, stop_fut, final_n_cur;
  std::vector<double> stop_time;
  std::vector<Registry> reg;
  void reset(const Cfg& c) {
    status.assign(c.n_arms, RECRUITING);
    enrolled.assign(c.n_arms, 0);
    stop_eff.assign(c.n_arms, 0);
    stop_fut.assign(c.n_arms, 0);
    final_n_cur.assign(c.n_arms, NA_INTEGER);
    stop_time.assign(c.n_arms, NA_REAL);
    reg.assign(c.n_arms, Registry());
    for (int a = 0; a < c.n_arms; a++) {
      reg[a].enroll.reserve(c.capacity[a]);
      reg[a].te.reserve(c.capacity[a]);
      reg[a].rc.reserve(c.capacity[a]);
    }
  }
};

// ---- look slicing (copy of src/simulate_chunk.cpp:252-324) ---------------------------------------------------------------
struct Slice {
  int n_pat = 0;
  int status_sum = 0;
  std::vector<int> events;
  std::vector<double> pt;
  int events_total = 0;
  double pt_total = 0.0;
  double median_fu = 0.0;
  std::vector<double> obs;
  std::vector<int> st;
  std::vector<double> contrib;
};

void slice_arm(const Registry& r, double cal, double max_fu, const std::vector<double>& cut, int K, Slice& s) {
  s.events.assign(K, 0);
  s.pt.assign(K, 0.0);
  s.obs.clear();
  s.st.clear();
  s.n_pat = 0;
  s.status_sum = 0;
  s.events_total = 0;
  s.pt_total = 0.0;
  s.median_fu = 0.0;
  if (r.n == 0) return;
  for (int i = 0; i < r.n; i++) {
    const double ts = r_pmax2(0.0, cal - r.enroll[i]);
    const double ta = r_pmin2(ts, max_fu);
    const double ob = r_pmin2(r_pmin2(r.te[i], r.rc[i]), ta);
    const int es = (r.te[i] <= r_pmin2(r.rc[i], ta)) ? 1 : 0;
    if (ta > 0) {
      s.obs.push_back(ob);
      s.st.push_back(es);
    }
  }
  s.n_pat = (int) s.obs.size();
  for (int e : s.st) s.status_sum += e;
  if (s.n_pat == 0) {
    s.events_total = 0;
    s.pt_total = evt_rsem::sum_ld(s.pt.data(), K);
    return;
  }
  const int nb = K + 1;
  for (int i = 0; i < s.n_pat; i++) {
    if (s.st[i] != 1) continue;
    int idx = evt_rsem::find_interval(cut.data(), nb, s.obs[i], false, false, false);
    if (idx == NA_INTEGER) Rcpp::stop("simulate_chunk_exact: NA event time in the interval tabulation");
    if (idx > 0) {
      if (idx > K) idx = K;
      s.events[idx - 1]++;
    }
  }
  for (int k = 0; k < K; k++) {
    const double lo = cut[k], hi = cut[k + 1];
    s.contrib.clear();
    for (int i = 0; i < s.n_pat; i++) {
      if (s.obs[i] > lo) {
        double c = r_pmin2(s.obs[i], hi) - lo;
        if (c < 0) c = 0;
        s.contrib.push_back(c);
      }
    }
    if (!s.contrib.empty()) s.pt[k] = evt_rsem::sum_ld(s.contrib.data(), (R_xlen_t) s.contrib.size());
  }
  long long et = 0;
  for (int k = 0; k < K; k++) et += s.events[k];
  s.events_total = (int) et;
  s.pt_total = evt_rsem::sum_ld(s.pt.data(), K);
  s.median_fu = evt_rsem::median_r(s.obs);
}

bool any_recruiting(const TrialState& st) {
  for (int s : st.status) if (s == RECRUITING) return true;
  return false;
}

bool gates_pass_vsref(const Cfg& c, const Slice& sC, const Slice& sT) {   // copy of src/simulate_chunk.cpp:491-501
  const double evC = (double) sC.events_total, evT = (double) sT.events_total;
  const double mfuC = sC.median_fu, mfuT = sT.median_fu;
  const double ptC = sC.pt_total, ptT = sT.pt_total;
  const double maxPT_C = c.vr_max_pt[0], maxPT_T = c.vr_max_pt[1];
  const double fracC = (maxPT_C > 0) ? ptC / maxPT_C : 0;
  const double fracT = (maxPT_T > 0) ? ptT / maxPT_T : 0;
  const bool passC = (evC >= c.vr_min_ev[0]) && (mfuC >= c.vr_min_mfu[0]) && (fracC >= c.vr_min_pt[0]);
  const bool passT = (evT >= c.vr_min_ev[1]) && (mfuT >= c.vr_min_mfu[1]) && (fracT >= c.vr_min_pt[1]);
  return passC && passT;
}

bool trt_can_enroll(const Cfg& c, const TrialState& st) {
  const int t = c.trt_idx;
  return st.status[t] == RECRUITING && (double) st.enrolled[t] < c.mtpa[t];
}

// ---- Rng policies ------------------------------------------------------------------------------------------------------
// RStreamRng: R's stream in the phase 1 order (src/simulate_chunk.cpp:613, 677, 680-682). Replay only.
struct RStreamRng {
  double arrival(double scale) { return R::rexp(scale); }
  int arm(const std::vector<double>& probs) { return prob_sample(probs, []() { return unif_rand(); }, true); }
  double event_time(double shape, double scale) { return R::rweibull(shape, scale); }
  double dropout_time(const Cfg& c) {
    double t = R::runif(0.0, c.censor_max_time_sim);   // drawn, then discarded when random_dropout is FALSE
    if (!c.random_dropout) t = R_PosInf;
    return t;
  }
};

// CounterRng: per-trial xoshiro256++ (header comment). The supported set has random_dropout FALSE (unpack stops
// otherwise), so dropout draws nothing and the censoring time is Inf, as phase 1 stores it.
struct CounterRng {
  Xoshiro256pp g;
  double arrival(double scale) { return -std::log(g.unif()) * scale; }
  int arm(const std::vector<double>& probs) { return prob_sample(probs, [this]() { return g.unif(); }, false); }
  double event_time(double shape, double scale) { return scale * std::pow(-std::log(g.unif()), 1. / shape); }
  double dropout_time(const Cfg&) { return R_PosInf; }
};

// ---- Post policies -----------------------------------------------------------------------------------------------------
// McPost: the phase 1 Monte Carlo posteriors (src/simulate_chunk.cpp:330-354, 397-488). Replay only.
struct Tally { int count; bool any_na; int n; };
Tally tally(const arma::vec& med, double thr, bool lt) {
  Tally t{0, false, (int) med.n_elem};
  for (arma::uword i = 0; i < med.n_elem; i++) {
    const double m = med(i);
    if (std::isnan(m)) { t.any_na = true; continue; }
    if (lt ? (m < thr) : (m > thr)) t.count++;
  }
  return t;
}
double tally_mean(const Tally& t) {
  if (t.any_na) return NA_REAL;
  return evt_rsem::mean_logical(t.count, (R_xlen_t) t.n);
}
double median_piecewise_r(const double* lam, const std::vector<double>& lens) {
  const int K = (int) lens.size();
  const double log2 = std::log(2.0);
  double cum = 0.0, t = 0.0;
  for (int k = 0; k < K; k++) {
    const double lj = lam[k], dt = lens[k];
    const double cend = cum + (lj * dt);
    if (ISNAN(cend)) Rcpp::stop("missing value where TRUE/FALSE needed (piecewise median)");
    if (cend >= log2) {
      if (lj == 0) return t;
      return t + (log2 - cum) / lj;
    }
    cum = cend;
    t = t + dt;
  }
  const double last = lam[K - 1];
  if (ISNAN(last)) Rcpp::stop("missing value where TRUE/FALSE needed (piecewise median)");
  if (last > 0) return t + (log2 - cum) / last;
  return R_PosInf;
}
double logical_mean(const std::vector<double>& x, double thr, bool lt) {
  int count = 0;
  for (double v : x) {
    if (std::isnan(v)) return NA_REAL;
    if (lt ? (v < thr) : (v >= thr)) count++;
  }
  return evt_rsem::mean_logical(count, (R_xlen_t) x.size());
}

struct McPost {
  arma::vec med;
  std::vector<double> logHR, shapes, lam, lamT;
  arma::vec medians(const Cfg& c, const Slice& s, const std::vector<double>& lens, int num_samples) {
    arma::ivec ev(c.K);
    for (int k = 0; k < c.K; k++) ev(k) = s.events[k];
    arma::vec pt(s.pt);
    arma::vec a(c.prior_alpha), b(c.prior_beta);
    arma::mat lm = draw_posterior_hazard_samples_cpp(c.K, ev, pt, a, b, num_samples);
    arma::vec il(lens);
    return calculate_median_survival_matrix_cpp(lm, il);
  }
  double hc_interim(const Cfg& c, const Slice& s, const std::vector<double>& lens, double thr) {
    arma::vec m = medians(c, s, lens, c.draws_interim);
    return tally_mean(tally(m, thr, true));
  }
  void hc_final_begin(const Cfg& c, const Slice& s, const std::vector<double>& lens) {
    med = medians(c, s, lens, c.draws_final);
  }
  double hc_final_above(double thr) { return tally_mean(tally(med, thr, false)); }
  double hc_final_below(double thr) { return tally_mean(tally(med, thr, true)); }
  // ph_sample_core (src/simulate_chunk.cpp:430-473) without keep_medians
  void ph(const Cfg& c, const Slice& sC, const Slice& sT, const std::vector<double>& lens, int n) {
    const int K = c.K;
    arma::ivec aEC(K), aET(K);
    arma::vec aPTC(K), aPTT(K);
    for (int k = 0; k < K; k++) {
      aEC(k) = sC.events[k];
      aET(k) = sT.events[k];
      aPTC(k) = sC.pt[k];
      aPTT(k) = sT.pt[k];
    }
    arma::vec aA(c.prior_alpha), aB(c.prior_beta);
    Rcpp::List bp = ph_beta_mode_var_cpp(aEC, aPTC, aET, aPTT, aA, aB, c.ph_mu, c.ph_sigma, c.ph_tol, c.ph_max_iter);
    const double bmean = Rcpp::as<double>(bp["mean"]);
    const double bsd = Rcpp::as<double>(bp["sd"]);
    logHR.resize(n);
    for (int i = 0; i < n; i++) logHR[i] = R::rnorm(bmean, bsd);
    shapes.resize(K);
    for (int k = 0; k < K; k++) shapes[k] = (c.prior_alpha[k] + (double) sC.events[k]) + (double) sT.events[k];
    lam.resize(K);
    lamT.resize(K);
    for (int i = 0; i < n; i++) {
      const double eb = std::exp(logHR[i]);
      for (int k = 0; k < K; k++) {
        const double rate = (c.prior_beta[k] + sC.pt[k]) + eb * sT.pt[k];
        const double scale = 1 / rate;
        lam[k] = R::rgamma(shapes[k], scale);
      }
      for (int k = 0; k < K; k++) lamT[k] = lam[k] * eb;
      (void) median_piecewise_r(lam.data(), lens);
      (void) median_piecewise_r(lamT.data(), lens);
    }
  }
  double ph_ge(double m) { return logical_mean(logHR, m, false); }   // mean(logHR >= m)
  double ph_lt(double m) { return logical_mean(logHR, m, true); }    // mean(logHR < m)
};

// ExactPost: the P2-1 exact probabilities. Nothing is drawn.
struct ExactPost {
  const Slice* fs = nullptr;
  const std::vector<double>* flens = nullptr;
  const Cfg* fc = nullptr;
  double bmean = 0, bsd = 1;
  bool warned = false;
  double below(const Cfg& c, const Slice& s, const std::vector<double>& lens, double thr, bool want_below) {
    evt_exact::MedianProbs p = evt_exact::median_probs(s.events.data(), s.pt.data(), c.prior_alpha.data(),
                                                       c.prior_beta.data(), lens.data(), c.K, thr);
    if (p.status == 1 && !warned) {   // phase 1 warns from draw_posterior_hazard_samples_cpp; unreachable in the set
      warned = true;
      Rcpp::warning("Posterior beta parameter is non-positive. Setting to 1e-6.");
    }
    return want_below ? p.below : p.above;   // status 2: NaN, handled by the callers as phase 1 handles an NA tally
  }
  double hc_interim(const Cfg& c, const Slice& s, const std::vector<double>& lens, double thr) {
    return below(c, s, lens, thr, true);   // Pr(med < thr)
  }
  void hc_final_begin(const Cfg& c, const Slice& s, const std::vector<double>& lens) {
    fc = &c;
    fs = &s;
    flens = &lens;
  }
  double hc_final_above(double thr) { return below(*fc, *fs, *flens, thr, false); }   // Pr(med > thr)
  double hc_final_below(double thr) { return below(*fc, *fs, *flens, thr, true); }    // Pr(med < thr)
  void ph(const Cfg& c, const Slice& sC, const Slice& sT, const std::vector<double>&, int) {
    const int K = c.K;
    arma::ivec aEC(K), aET(K);
    arma::vec aPTC(K), aPTT(K);
    for (int k = 0; k < K; k++) {
      aEC(k) = sC.events[k];
      aET(k) = sT.events[k];
      aPTC(k) = sC.pt[k];
      aPTT(k) = sT.pt[k];
    }
    arma::vec aA(c.prior_alpha), aB(c.prior_beta);
    Rcpp::List bp = ph_beta_mode_var_cpp(aEC, aPTC, aET, aPTT, aA, aB, c.ph_mu, c.ph_sigma, c.ph_tol, c.ph_max_iter);
    bmean = Rcpp::as<double>(bp["mean"]);
    bsd = Rcpp::as<double>(bp["sd"]);
  }
  double ph_ge(double m) { return evt_exact::loghr_probs(bmean, bsd, m).ge; }   // Pr(logHR >= m)
  double ph_lt(double m) { return evt_exact::loghr_probs(bmean, bsd, m).lt; }   // Pr(logHR < m)
};

// ---- interim looks (src/simulate_chunk.cpp:357-385, 505-522), posterior through the policy --------------------------------
template <class Post>
void interim_hc(const Cfg& c, TrialState& st, double cal, const std::vector<double>& cut, const std::vector<double>& lens,
                Slice& sl, Post& post) {
  for (int a = 0; a < c.n_arms; a++) {
    if (st.status[a] != RECRUITING) continue;
    slice_arm(st.reg[a], cal, c.max_follow_up_sim, cut, c.K, sl);
    const double n_pat = (double) sl.n_pat;
    const double events_total = (double) sl.status_sum;
    const double median_fu = (sl.n_pat > 0) ? sl.median_fu : 0.0;
    const double pt_total = evt_rsem::sum_ld(sl.pt.data(), c.K);
    const double max_pt = c.max_pt_arm[a];
    const double pt_frac = (max_pt > 0) ? pt_total / max_pt : 0.0;
    if (n_pat < c.min_pat[a] || events_total < c.min_events[a] || median_fu < c.min_median_fu[a] ||
        pt_frac < c.min_pt_frac[a]) {
      continue;
    }
    const double pr_fut = post.hc_interim(c, sl, lens, c.futility_thr_interim[a]);
    if (ISNAN(pr_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (interim futility probability is NA)");
    if (pr_fut >= c.futility_prob_interim) {
      st.status[a] = STOPPED_FUTILITY;
      st.stop_fut[a] = 1;
      st.final_n_cur[a] = st.enrolled[a];
      st.stop_time[a] = cal;
    }
  }
}

template <class Post>
void interim_vsref(const Cfg& c, TrialState& st, double cal, const std::vector<double>& cut,
                   const std::vector<double>& lens, Slice& sC, Slice& sT, Post& post) {
  const int r = c.ref_idx, t = c.trt_idx;
  if (st.status[t] != RECRUITING) return;
  slice_arm(st.reg[r], cal, c.max_follow_up_sim, cut, c.K, sC);
  slice_arm(st.reg[t], cal, c.max_follow_up_sim, cut, c.K, sT);
  if (!gates_pass_vsref(c, sC, sT)) return;
  post.ph(c, sC, sT, lens, c.draws_vsref_interim);
  const double pr_fut = post.ph_ge(c.log_margin_interim);
  if (ISNAN(pr_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (interim vs-ref futility probability is NA)");
  if (pr_fut >= c.vr_fut_prob_interim) {
    st.status[t] = STOPPED_FUTILITY;
    st.stop_fut[t] = 1;
    st.final_n_cur[t] = st.enrolled[t];
    st.stop_time[t] = cal;
  }
}

// ---- one trial (the body of the loop at src/simulate_chunk.cpp:595-777) ------------------------------------------------
struct Work {
  TrialState st;
  Slice sl, sl_ref, sl_trt;
  std::vector<double> cut, lens, ev_times, probs;
  std::vector<int> elig;
};

struct TrialOut {
  std::vector<double> final_n, final_ev;
  std::vector<int> stop_eff, stop_fut, final_eff, final_fut, final_inc;
  double final_time = 0;
};

template <class Rng, class Post>
void run_trial(const Cfg& c, Rng& rng, Post& post, Work& w, TrialOut& out) {
  const int A = c.n_arms;
  const int K = c.K;
  TrialState& st = w.st;
  std::vector<double>& cut = w.cut;
  std::vector<double>& lens = w.lens;
  const double rexp_scale = 1 / c.accrual_rate;

  cut = c.cut_base;
  lens = c.len_base;
  bool rebalance_done = false;
  double last_enroll_time = 0;
  st.reset(c);
  double current_time = 0.0;
  double next_calendar_look = c.beat;

  auto eligible_any = [&]() {
    for (int a = 0; a < A; a++)
      if (st.status[a] == RECRUITING && (double) st.enrolled[a] < c.mtpa[a]) return true;
    return false;
  };

  // ---- accrual loop ----
  while (eligible_any()) {
    const double interarrival = rng.arrival(rexp_scale);
    current_time = current_time + interarrival;

    if (current_time >= c.max_trial_time) {
      current_time = last_enroll_time;
      break;
    }
    if (c.vsref && !trt_can_enroll(c, st)) {
      current_time = last_enroll_time;
      break;
    }

    if (!rebalance_done) {
      int total_events = 0;
      w.ev_times.clear();
      for (int a = 0; a < A; a++) {
        const Registry& r = st.reg[a];
        int count = 0;
        for (int i = 0; i < r.n; i++) {
          const double ts = r_pmax2(0.0, current_time - r.enroll[i]);
          const double ta = r_pmin2(ts, c.max_follow_up_sim);
          if (!(ta > 0)) continue;
          const double ob = r_pmin2(r_pmin2(r.te[i], r.rc[i]), ta);
          if (r.te[i] <= r_pmin2(r.rc[i], ta)) {
            count++;
            w.ev_times.push_back(ob);
          }
        }
        total_events = total_events + count;
      }
      if ((double) total_events >= c.rebalance_threshold) {
        std::vector<double> new_cuts = evt_rsem::rebalance_cutpoints(w.ev_times, c.max_follow_up_sim, K);
        if (!new_cuts.empty()) {
          cut = new_cuts;
          for (int k = 0; k < K; k++) lens[k] = cut[k + 1] - cut[k];
        }
        rebalance_done = true;
      }
    }

    if (current_time >= next_calendar_look) {
      if (c.vsref) interim_vsref(c, st, current_time, cut, lens, w.sl_ref, w.sl_trt, post);
      else interim_hc(c, st, current_time, cut, lens, w.sl, post);
      next_calendar_look = next_calendar_look + c.beat;
    }

    if (!any_recruiting(st)) break;
    w.elig.clear();
    for (int a = 0; a < A; a++)
      if (st.status[a] == RECRUITING && (double) st.enrolled[a] < c.mtpa[a]) w.elig.push_back(a);
    if (w.elig.empty()) break;
    if (c.vsref && !trt_can_enroll(c, st)) break;

    w.probs.resize(w.elig.size());
    for (size_t i = 0; i < w.elig.size(); i++) w.probs[i] = c.rand_probs[w.elig[i]];
    const double psum = evt_rsem::sum_ld(w.probs.data(), (R_xlen_t) w.probs.size());
    for (size_t i = 0; i < w.probs.size(); i++) w.probs[i] = w.probs[i] / psum;
    const int chosen = w.elig[rng.arm(w.probs)];

    const double t_event_true = rng.event_time(c.shape[chosen], c.scale[chosen]);
    const double t_random_censor = rng.dropout_time(c);

    Registry& r = st.reg[chosen];
    const int idx = r.n + 1;
    if (idx > c.capacity[chosen]) {
      Rcpp::stop("Registry overflow for arm '%s': idx=%d exceeds capacity=%d. This may indicate max_total_patients_per_arm "
                 "is too low or a logic error.", Rcpp::as<std::string>(c.arm_names[chosen]), idx, c.capacity[chosen]);
    }
    r.enroll.push_back(current_time);
    r.te.push_back(t_event_true);
    r.rc.push_back(t_random_censor);
    r.n = idx;
    st.enrolled[chosen] = st.enrolled[chosen] + 1;
    last_enroll_time = current_time;
  }

  // ---- post-enrollment looks, looks_stop_at_final = TRUE ----
  const double t_star = std::min(last_enroll_time + c.min_follow_up_at_final, c.max_trial_time);
  double post_enroll_time = last_enroll_time;
  while (R_FINITE(t_star) && any_recruiting(st)) {
    post_enroll_time = post_enroll_time + c.beat;
    if (post_enroll_time >= t_star) break;
    if (post_enroll_time <= current_time) continue;
    if (c.vsref) interim_vsref(c, st, post_enroll_time, cut, lens, w.sl_ref, w.sl_trt, post);
    else interim_hc(c, st, post_enroll_time, cut, lens, w.sl, post);
  }

  // ---- final analysis time ----
  const bool all_stopped_early = !any_recruiting(st);
  bool any_stop_time = false;
  double max_stop = R_NegInf;
  for (int a = 0; a < A; a++) {
    if (!ISNAN(st.stop_time[a])) {
      any_stop_time = true;
      if (st.stop_time[a] > max_stop) max_stop = st.stop_time[a];
    }
  }
  double final_time;
  if (all_stopped_early && any_stop_time) {
    final_time = max_stop;
  } else {
    final_time = std::min(last_enroll_time + c.min_follow_up_at_final, c.max_trial_time);
  }

  out.final_n.assign(A, 0.0);
  out.final_ev.assign(A, 0.0);
  out.final_eff.assign(A, 0);
  out.final_fut.assign(A, 0);
  out.final_inc.assign(A, 0);
  Slice& sl = w.sl;
  if (c.vsref) {
    slice_arm(st.reg[c.ref_idx], final_time, c.max_follow_up_sim, cut, K, w.sl_ref);
    for (int a = 0; a < A; a++) {
      slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);
      out.final_ev[a] = (double) sl.events_total;
      if (st.status[a] != RECRUITING) continue;
      if (a == c.ref_idx) {
        out.final_inc[a] = 1;
        continue;
      }
      post.ph(c, w.sl_ref, sl, lens, c.draws_final);
      const double p_eff = post.ph_lt(c.log_eff_final);
      const double p_fut = post.ph_ge(c.log_margin_final);
      if (ISNAN(p_eff)) Rcpp::stop("missing value where TRUE/FALSE needed (final vs-ref efficacy probability is NA)");
      if (p_eff >= c.vr_eff_prob_final) {
        out.final_eff[a] = 1;
      } else {
        if (ISNAN(p_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (final vs-ref futility probability is NA)");
        if (p_fut >= c.vr_fut_prob_final) out.final_fut[a] = 1;
        else out.final_inc[a] = 1;
      }
    }
  }
  for (int a = 0; a < A && !c.vsref; a++) {
    if (st.status[a] != RECRUITING) {
      slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);
      out.final_ev[a] = (double) sl.events_total;
      continue;
    }
    slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);
    out.final_ev[a] = (double) sl.events_total;
    post.hc_final_begin(c, sl, lens);
    const double p_eff = post.hc_final_above(c.success_thr[a]);
    if (!ISNAN(p_eff) && p_eff >= c.final_success_prob) {
      out.final_eff[a] = 1;
    } else {
      const double p_fut = post.hc_final_below(c.futility_thr[a]);
      if (!ISNAN(p_fut) && p_fut >= c.final_futility_prob) {
        out.final_fut[a] = 1;
      } else {
        out.final_inc[a] = 1;
      }
    }
  }
  for (int a = 0; a < A; a++) {
    out.final_n[a] = (st.final_n_cur[a] == NA_INTEGER) ? (double) st.enrolled[a] : (double) st.final_n_cur[a];
  }
  out.stop_eff = st.stop_eff;
  out.stop_fut = st.stop_fut;
  out.final_time = final_time;
}

// ---- phase 1 accumulators and return list (src/simulate_chunk.cpp:569-577, 780-836) -------------------------------------
struct Acc {
  int A;
  R_xlen_t n_sims;
  bool collect_raw_n, return_variance;
  std::vector<double> sum_n, sum_ev, sum_se, sum_sf, sum_fe, sum_ff, sum_fi;
  double sum_time = 0;
  std::vector<std::vector<double>> raw_n, raw_ev;
  std::vector<double> raw_time;
  Acc(const Cfg& c, R_xlen_t n)
      : A(c.n_arms), n_sims(n), collect_raw_n(c.collect_raw_n), return_variance(c.return_variance),
        sum_n(A, 0.0), sum_ev(A, 0.0), sum_se(A, 0.0), sum_sf(A, 0.0), sum_fe(A, 0.0), sum_ff(A, 0.0), sum_fi(A, 0.0) {
    if (collect_raw_n) raw_n.assign(A, std::vector<double>(n, 0.0));
    if (return_variance) {
      raw_ev.assign(A, std::vector<double>(n, 0.0));
      raw_time.assign(n, 0.0);
    }
  }
  // the phase 1 per-trial accumulation, plain double +
  void add(const double* fn, const double* fev, const int* se, const int* sf, const int* fe, const int* ff, const int* fi,
           double ft) {
    for (int a = 0; a < A; a++) {
      sum_n[a] = sum_n[a] + fn[a];
      sum_ev[a] = sum_ev[a] + fev[a];
      sum_se[a] = sum_se[a] + (double) se[a];
      sum_sf[a] = sum_sf[a] + (double) sf[a];
      sum_fe[a] = sum_fe[a] + (double) fe[a];
      sum_ff[a] = sum_ff[a] + (double) ff[a];
      sum_fi[a] = sum_fi[a] + (double) fi[a];
    }
    sum_time = sum_time + ft;
  }
  void store(R_xlen_t pos, const double* fn, const double* fev, double ft) {
    if (collect_raw_n) for (int a = 0; a < A; a++) raw_n[a][pos] = fn[a];
    if (return_variance) {
      for (int a = 0; a < A; a++) raw_ev[a][pos] = fev[a];
      raw_time[pos] = ft;
    }
  }
  Rcpp::List list(const Rcpp::CharacterVector& arm_names) const {
    auto named = [&](const std::vector<double>& v) {
      Rcpp::NumericVector out(v.begin(), v.end());
      out.attr("names") = arm_names;
      return out;
    };
    auto raw_list = [&](const std::vector<std::vector<double>>& v) -> SEXP {
      Rcpp::List out(A);
      for (int a = 0; a < A; a++) out[a] = Rcpp::NumericVector(v[a].begin(), v[a].end());
      out.attr("names") = arm_names;
      return out;
    };
    SEXP final_n_raw = collect_raw_n ? raw_list(raw_n) : R_NilValue;
    SEXP final_events_raw = return_variance ? raw_list(raw_ev) : R_NilValue;
    SEXP final_time_raw = return_variance ? (SEXP) Rcpp::NumericVector(raw_time.begin(), raw_time.end()) : R_NilValue;
    return Rcpp::List::create(
        Rcpp::Named("sum_final_n") = named(sum_n),
        Rcpp::Named("sum_final_events") = named(sum_ev),
        Rcpp::Named("sum_stop_eff") = named(sum_se),
        Rcpp::Named("sum_stop_fut") = named(sum_sf),
        Rcpp::Named("sum_final_eff") = named(sum_fe),
        Rcpp::Named("sum_final_fut") = named(sum_ff),
        Rcpp::Named("sum_final_inc") = named(sum_fi),
        Rcpp::Named("sum_final_time") = Rcpp::NumericVector::create(sum_time),
        Rcpp::Named("n_sims") = Rcpp::IntegerVector::create((int) n_sims),
        Rcpp::Named("final_n_raw") = final_n_raw,
        Rcpp::Named("final_events_raw") = final_events_raw,
        Rcpp::Named("final_time_raw") = final_time_raw);
  }
};

uint64_t seed64_from(const Rcpp::IntegerVector& s) {
  if (s.size() != 2 || s[0] == NA_INTEGER || s[1] == NA_INTEGER) {
    Rcpp::stop("simulate_chunk_exact: scenario_seed must be two non-NA integers");
  }
  return ((uint64_t) (uint32_t) s[0] << 32) | (uint64_t) (uint32_t) s[1];
}

}  // namespace

// [[Rcpp::export(name = ".simulate_chunk_exact", rng = false)]]
Rcpp::List simulate_chunk_exact(Rcpp::IntegerVector sim_indices, Rcpp::List cfg, Rcpp::IntegerVector scenario_seed) {
  // Phase 2 engine. sim_indices: trial indices (1-based, non-NA); cfg: the phase 1 packed configuration;
  // scenario_seed: two integers (the driver draws them with one sample.int call). Returns per-trial outcomes in the order
  // of sim_indices: sim_indices, the n x A matrices final_n, final_events (double) and stop_eff, stop_fut, final_eff,
  // final_fut, final_inc (integer), and final_time (length n). R's generator is not touched (rng = false).
  const Cfg c = unpack(cfg);
  const uint64_t seed64 = seed64_from(scenario_seed);
  const int A = c.n_arms;
  const R_xlen_t n = sim_indices.size();
  Rcpp::NumericMatrix fn(n, A), fev(n, A);
  Rcpp::IntegerMatrix se(n, A), sf(n, A), fe(n, A), ff(n, A), fi(n, A);
  Rcpp::NumericVector ft(n);
  Work w;
  TrialOut out;
  CounterRng rng;
  ExactPost post;
  for (R_xlen_t t = 0; t < n; t++) {
    const int idx = sim_indices[t];
    if (idx == NA_INTEGER) Rcpp::stop("simulate_chunk_exact: NA trial index");
    rng.g.seed(seed64, (uint64_t) (int64_t) idx);
    run_trial(c, rng, post, w, out);
    for (int a = 0; a < A; a++) {
      fn(t, a) = out.final_n[a];
      fev(t, a) = out.final_ev[a];
      se(t, a) = out.stop_eff[a];
      sf(t, a) = out.stop_fut[a];
      fe(t, a) = out.final_eff[a];
      ff(t, a) = out.final_fut[a];
      fi(t, a) = out.final_inc[a];
    }
    ft[t] = out.final_time;
  }
  Rcpp::List dn = Rcpp::List::create(R_NilValue, c.arm_names);
  for (SEXP m : {(SEXP) fn, (SEXP) fev, (SEXP) se, (SEXP) sf, (SEXP) fe, (SEXP) ff, (SEXP) fi}) Rf_setAttrib(m, R_DimNamesSymbol, dn);
  return Rcpp::List::create(Rcpp::Named("sim_indices") = Rcpp::IntegerVector(sim_indices.begin(), sim_indices.end()),
                            Rcpp::Named("final_n") = fn, Rcpp::Named("final_events") = fev,
                            Rcpp::Named("stop_eff") = se, Rcpp::Named("stop_fut") = sf,
                            Rcpp::Named("final_eff") = fe, Rcpp::Named("final_fut") = ff,
                            Rcpp::Named("final_inc") = fi, Rcpp::Named("final_time") = ft);
}

// [[Rcpp::export(name = ".exact_reduce", rng = false)]]
Rcpp::List exact_reduce(Rcpp::List blocks, Rcpp::List cfg, int n_sims) {
  // Combines .simulate_chunk_exact block results into the phase 1 chunk list for trials 1..n_sims. Each index must occur
  // exactly once over the blocks. Sums are taken in trial-index order with the phase 1 accumulation (plain double +), and
  // the raw vectors hold trial i at position i, so the result does not depend on how the trials were split or ordered.
  const Cfg c = unpack(cfg);
  const int A = c.n_arms;
  if (n_sims < 1) Rcpp::stop("exact_reduce: n_sims must be positive");
  std::vector<int> blk(n_sims, -1), row(n_sims, -1);
  const R_xlen_t nb = blocks.size();
  std::vector<Rcpp::List> bl(nb);
  for (R_xlen_t b = 0; b < nb; b++) {
    bl[b] = Rcpp::as<Rcpp::List>(blocks[b]);
    Rcpp::IntegerVector ix = bl[b]["sim_indices"];
    for (R_xlen_t i = 0; i < ix.size(); i++) {
      const int k = ix[i];
      if (k == NA_INTEGER || k < 1 || k > n_sims) Rcpp::stop("exact_reduce: trial index out of range");
      if (blk[k - 1] >= 0) Rcpp::stop("exact_reduce: trial index %d occurs twice", k);
      blk[k - 1] = (int) b;
      row[k - 1] = (int) i;
    }
  }
  for (int k = 0; k < n_sims; k++) if (blk[k] < 0) Rcpp::stop("exact_reduce: trial index %d missing", k + 1);
  struct V {
    Rcpp::NumericMatrix fn, fev;
    Rcpp::IntegerMatrix se, sf, fe, ff, fi;
    Rcpp::NumericVector ft;
  };
  std::vector<V> v(nb);
  for (R_xlen_t b = 0; b < nb; b++) {
    v[b] = V{bl[b]["final_n"], bl[b]["final_events"], bl[b]["stop_eff"], bl[b]["stop_fut"], bl[b]["final_eff"],
             bl[b]["final_fut"], bl[b]["final_inc"], bl[b]["final_time"]};
    if (v[b].fn.ncol() != A) Rcpp::stop("exact_reduce: block has the wrong number of arms");
  }
  Acc acc(c, n_sims);
  std::vector<double> fn(A), fev(A);
  std::vector<int> se(A), sf(A), fe(A), ff(A), fi(A);
  for (int k = 0; k < n_sims; k++) {
    const V& b = v[blk[k]];
    const int i = row[k];
    for (int a = 0; a < A; a++) {
      fn[a] = b.fn(i, a);
      fev[a] = b.fev(i, a);
      se[a] = b.se(i, a);
      sf[a] = b.sf(i, a);
      fe[a] = b.fe(i, a);
      ff[a] = b.ff(i, a);
      fi[a] = b.fi(i, a);
    }
    acc.add(fn.data(), fev.data(), se.data(), sf.data(), fe.data(), ff.data(), fi.data(), b.ft[i]);
    acc.store(k, fn.data(), fev.data(), b.ft[i]);
  }
  return acc.list(c.arm_names);
}

// [[Rcpp::export(name = ".simulate_chunk_exact_replay")]]
Rcpp::List simulate_chunk_exact_replay(Rcpp::IntegerVector sim_indices, Rcpp::List cfg) {
  // Test hook: the exact engine's trial loop with R's stream in the phase 1 order and the phase 1 Monte Carlo posteriors.
  // On the same .Random.seed it must return exactly what .simulate_chunk_cpp(sim_indices, cfg) returns and leave the
  // same .Random.seed (tests/testthat/test-simulate-chunk-exact.R).
  const Cfg c = unpack(cfg);
  const int A = c.n_arms;
  const R_xlen_t n_sims = sim_indices.size();
  Acc acc(c, n_sims);
  bool has_dup = false;
  {
    std::vector<int> s(sim_indices.begin(), sim_indices.end());
    std::sort(s.begin(), s.end());
    for (size_t i = 1; i < s.size(); i++) if (s[i] == s[i - 1]) { has_dup = true; break; }
  }
  Work w;
  TrialOut out;
  RStreamRng rng;
  McPost post;
  for (R_xlen_t t = 0; t < n_sims; t++) {
    run_trial(c, rng, post, w, out);
    acc.add(out.final_n.data(), out.final_ev.data(), out.stop_eff.data(), out.stop_fut.data(), out.final_eff.data(),
            out.final_fut.data(), out.final_inc.data(), out.final_time);
    const int sim_idx = sim_indices[t];
    if (sim_idx != NA_INTEGER) {
      if (!has_dup) {
        acc.store(t, out.final_n.data(), out.final_ev.data(), out.final_time);
      } else {
        for (R_xlen_t j = 0; j < n_sims; j++)
          if (sim_indices[j] == sim_idx) acc.store(j, out.final_n.data(), out.final_ev.data(), out.final_time);
      }
    }
  }
  (void) A;
  return acc.list(c.arm_names);
}

// [[Rcpp::export(name = ".exact_rng_uniforms", rng = false)]]
Rcpp::NumericVector exact_rng_uniforms(Rcpp::IntegerVector scenario_seed, int trial, int n) {
  // Test hook: the first n uniforms of trial `trial`'s generator.
  Xoshiro256pp g;
  g.seed(seed64_from(scenario_seed), (uint64_t) (int64_t) trial);
  Rcpp::NumericVector u(n);
  for (int i = 0; i < n; i++) u[i] = g.unif();
  return u;
}

// [[Rcpp::export(name = ".exact_rng_raw", rng = false)]]
Rcpp::NumericVector exact_rng_raw(Rcpp::NumericVector state, int n) {
  // Test hook: n xoshiro256++ outputs from an explicit state (four words given as doubles holding 32-bit halves:
  // state = c(hi0, lo0, hi1, lo1, ...)), each output returned as c(hi, lo). For the reference-vector test.
  if (state.size() != 8) Rcpp::stop("state must have 8 entries (hi, lo per word)");
  Xoshiro256pp g;
  for (int i = 0; i < 4; i++) g.s[i] = ((uint64_t) state[2 * i] << 32) | (uint64_t) state[2 * i + 1];
  Rcpp::NumericVector out(2 * n);
  for (int i = 0; i < n; i++) {
    const uint64_t x = g.next();
    out[2 * i] = (double) (x >> 32);
    out[2 * i + 1] = (double) (x & 0xFFFFFFFFULL);
  }
  return out;
}
