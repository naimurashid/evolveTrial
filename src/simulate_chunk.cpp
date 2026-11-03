// simulate_chunk.cpp -- bit-identical C++ port of the single-arm (HC) chunk body of run_simulation_pure (engine C++ port,
// step 3; plan dev/plans/2026-11-03-engine-cpp-plan.md, design Section 2, inventory dev/perf/evidence/chunk_body_inventory.md).
//
// The R original is the closure `simulate_chunk_r` in R/simulation_driver.R (lines cited as SD:<line>), with its callees in
// R/interim_logic.R (IL), R/state_management.R (SM), R/posterior_cpp_dispatchers.R (PD) and R/gate_diagnostics.R (GD).
// Every random draw goes through R's own stream (R:: samplers and unif_rand) in the R body's order, inside the one RNGScope
// of the generated wrapper. Every sum, mean, median, quantile and findInterval goes through src/r_semantics.h. The gamma
// sampler and the piecewise median are the existing C++ functions in src/posterior_sampling.cpp, called directly.
//
// Scope: exactly the configuration set checked by .chunk_cpp_supported() in R/simulate_chunk_cpp.R. Common to both
// paths: looks_stop_at_final = TRUE, calendar-beat looks, rebalancing on, dropout off, finite max_trial_time. Step 3:
// single arm, HC analysis, interim efficacy stopping off. Step 4 (cfg$mode "vsref"): two arms, one reference and one
// experimental, vs-reference analysis on the PH model (the live R sampler sample_vs_ref_medians_ph, PH:105-135, with
// ph_beta_mode_var_cpp called unchanged), HR futility margin set, efficacy margin NULL, interim efficacy stopping off. The R side refuses anything else before this function is called; the checks below that can be
// made cheaply here stop() as a second line, they never fall back.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Random.h>   // unif_rand
#include <R_ext/Utils.h>    // revsort
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "r_semantics.h"

// Existing bodies in src/posterior_sampling.cpp (exported there; the generated RcppExports.cpp declares them the same way).
arma::mat draw_posterior_hazard_samples_cpp(int num_intervals, const arma::ivec& events_per_interval,
                                            const arma::vec& person_time_per_interval,
                                            const arma::vec& prior_alpha_params, const arma::vec& prior_beta_params,
                                            int num_samples);
arma::vec calculate_median_survival_matrix_cpp(const arma::mat& hazard_samples, const arma::vec& interval_lengths);
// The one definition is src/posterior_sampling.cpp:196 (step 4). It is what the R dispatcher ph_beta_mode_var
// (R/posterior_cpp_dispatchers.R:98-118) reaches through RcppExports. src/hybrid_trial_sim.cpp:88 is a declaration of
// the same function, not a second definition. No default arguments here: every argument is passed explicitly.
Rcpp::List ph_beta_mode_var_cpp(const arma::ivec& E_C, const arma::vec& PT_C, const arma::ivec& E_T, const arma::vec& PT_T,
                                const arma::vec& alpha_prior, const arma::vec& beta_prior, double mu, double sigma,
                                double tol, int max_iter);

namespace {

// ---- sample(x, size = 1, prob = p) ----------------------------------------------------------------------------------
// R 4.6.0 base::sample -> sample.int -> .Internal(sample(n, 1, FALSE, prob)) -> src/main/random.c do_sample():
// FixupProb(p, n, k = 1, replace = FALSE) (plain double sum of the positive entries, then p[i] /= sum), then
// ProbSampleNoReplace(n, p, perm, nans = 1, ans): perm = 1..n, revsort(p, perm, n) (descending, not stable),
// rT = totalmass * unif_rand() with totalmass = 1, scan j < n - 1 accumulating mass until rT <= mass.
// Exactly one unif_rand(), also when n == 1 (inventory 4.5, 6.3; plan P2). Returns a 0-based index into p.
int prob_sample_one(std::vector<double> p) {
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
  for (int i = 0; i < n; i++) p[i] /= sum;
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++) perm[i] = i + 1;
  revsort(p.data(), perm.data(), n);
  double totalmass = 1;
  int n1 = n - 1;
  double rT = totalmass * unif_rand();
  double mass = 0;
  int j;
  for (j = 0; j < n1; j++) {
    mass += p[j];
    if (rT <= mass) break;
  }
  return perm[j] - 1;
}

// ---- R's pmin / pmax on two doubles (src/main/summary.c do_pmin, na.rm = FALSE, no NaN on this path) ------------------
// pmax(a, b): start from a, replace when b > current. pmin(a, b): replace when b < current. Order of arguments as in R.
inline double r_pmax2(double a, double b) { return (b > a) ? b : a; }
inline double r_pmin2(double a, double b) { return (b < a) ? b : a; }

// ---- configuration -------------------------------------------------------------------------------------------------
struct Cfg {
  int n_arms;
  Rcpp::CharacterVector arm_names;
  std::vector<double> mtpa;           // max_total_patients_per_arm, positional (is_eligible, SD:538-540)
  std::vector<int> capacity;          // as.integer(ceiling(mtpa * 1.2)) (SM:74)
  std::vector<double> rand_probs;     // randomization_probs[arm_names] (by name, SD:650)
  std::vector<double> shape, scale;   // weibull_*_true_arms[arm] (by name, SD:655-656)
  double accrual_rate, beat, max_trial_time, min_follow_up_at_final;
  double max_follow_up_sim, censor_max_time_sim;
  bool random_dropout;
  int K;                              // num_intervals
  std::vector<double> cut_base, len_base;
  double rebalance_threshold;
  std::vector<double> prior_alpha, prior_beta;
  int draws_interim, draws_final;
  // HC interim gates and rule (IL:182-253), per arm
  std::vector<double> min_pat, min_events, min_median_fu, min_pt_frac, max_pt_arm;
  double futility_prob_interim;
  // thresholds on the median, per arm: the interim futility threshold as IL:57-61 resolves it from args, and the final
  // success and futility thresholds from the frame (SD:770-786). They come from different variables (args has no
  // median_pfs_futility_threshold_arms element), so they are kept apart.
  std::vector<double> futility_thr_interim, success_thr, futility_thr;
  double final_success_prob, final_futility_prob;
  bool collect_raw_n, return_variance;
  // vs-reference path (step 4): mode "vsref", two arms, one reference, one experimental
  bool vsref = false;
  int ref_idx = -1, trt_idx = -1;                     // 0-based positions in arm_names
  std::vector<double> vr_min_ev, vr_min_mfu, vr_min_pt, vr_max_pt;   // gates, {ctrl, trt} (SM:186-240)
  int draws_vsref_interim = 0;                        // args$num_posterior_draws (IL:7)
  double ph_mu = 0, ph_sigma = 1, ph_tol = 0;         // PH:113-114, PD:98-99
  int ph_max_iter = 0;
  double log_margin_interim = 0, log_eff_interim = 0, log_margin_final = 0, log_eff_final = 0;
  double vr_fut_prob_interim = 0, vr_eff_prob_final = 0, vr_fut_prob_final = 0;
};

double num1(const Rcpp::List& cfg, const char* name) {
  Rcpp::NumericVector v = cfg[name];
  if (v.size() != 1) Rcpp::stop(std::string("simulate_chunk_cpp: cfg$") + name + " must have length 1");
  return v[0];
}
int int1(const Rcpp::List& cfg, const char* name) {
  Rcpp::IntegerVector v = cfg[name];
  if (v.size() != 1 || v[0] == NA_INTEGER) Rcpp::stop(std::string("simulate_chunk_cpp: cfg$") + name + " must be one integer");
  return v[0];
}
bool lgl1(const Rcpp::List& cfg, const char* name) {
  Rcpp::LogicalVector v = cfg[name];
  if (v.size() != 1 || v[0] == NA_LOGICAL) Rcpp::stop(std::string("simulate_chunk_cpp: cfg$") + name + " must be TRUE or FALSE");
  return v[0] != 0;
}
std::vector<double> numv(const Rcpp::List& cfg, const char* name, int len) {
  Rcpp::NumericVector v = cfg[name];
  if (len >= 0 && v.size() != len) Rcpp::stop(std::string("simulate_chunk_cpp: cfg$") + name + " has the wrong length");
  return std::vector<double>(v.begin(), v.end());
}

Cfg unpack(const Rcpp::List& cfg) {
  if (!cfg.containsElementNamed("supported") || !Rcpp::as<bool>(cfg["supported"])) {
    Rcpp::stop("simulate_chunk_cpp: the packed configuration is outside the supported set (no fallback here)");
  }
  Cfg c;
  c.arm_names = cfg["arm_names"];
  c.n_arms = c.arm_names.size();
  const std::string mode = cfg.containsElementNamed("mode") ? Rcpp::as<std::string>(cfg["mode"]) : std::string("hc");
  if (mode == "vsref") {
    c.vsref = true;
    if (c.n_arms != 2) Rcpp::stop("simulate_chunk_cpp: the vs-reference configuration needs exactly two arms");
  } else if (mode == "hc") {
    if (c.n_arms != 1) Rcpp::stop("simulate_chunk_cpp: the single-arm configuration needs exactly one arm");
  } else {
    Rcpp::stop("simulate_chunk_cpp: unknown cfg$mode");
  }
  c.mtpa = numv(cfg, "mtpa", c.n_arms);
  Rcpp::IntegerVector cap = cfg["capacity"];
  if (cap.size() != c.n_arms) Rcpp::stop("simulate_chunk_cpp: cfg$capacity has the wrong length");
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
  if (c.random_dropout) Rcpp::stop("simulate_chunk_cpp: random_dropout = TRUE is outside the supported set");
  if (!lgl1(cfg, "looks_stop_at_final")) Rcpp::stop("simulate_chunk_cpp: looks_stop_at_final = FALSE is outside the supported set");
  if (!R_FINITE(c.max_trial_time)) Rcpp::stop("simulate_chunk_cpp: infinite max_trial_time is outside the supported set");
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
      Rcpp::stop("simulate_chunk_cpp: cfg$ref_idx and cfg$trt_idx must be the two distinct arm positions");
    }
    c.vr_min_ev = numv(cfg, "vr_min_ev", 2);
    c.vr_min_mfu = numv(cfg, "vr_min_mfu", 2);
    c.vr_min_pt = numv(cfg, "vr_min_pt", 2);
    c.vr_max_pt = numv(cfg, "vr_max_pt", 2);
    c.draws_vsref_interim = int1(cfg, "draws_vsref_interim");
    if (c.draws_vsref_interim < 1 || c.draws_final < 1) Rcpp::stop("simulate_chunk_cpp: posterior draw counts must be >= 1");
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
  if (c.draws_interim < 1 || c.draws_final < 1) Rcpp::stop("simulate_chunk_cpp: posterior draw counts must be >= 1");
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

// ---- per-trial state (make_state, SM:70-104) -----------------------------------------------------------------------
enum Status { RECRUITING = 0, STOPPED_EFFICACY = 1, STOPPED_FUTILITY = 2 };

struct Registry {                 // the active rows of the preallocated data frame (rows with non-NA id, SM:108-112)
  std::vector<double> enroll, te, rc;
  int n = 0;
};

struct TrialState {
  std::vector<int> status, enrolled, stop_eff, stop_fut, final_n_cur;  // final_n_cur: NA_INTEGER when unset
  std::vector<double> stop_time;                                        // NA_REAL when unset
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

// ---- slice_arm_data_at_time + calculate_interval_metrics_fast (SM:115-149, SM:11-63) -----------------------------
struct Slice {
  int n_pat = 0;                     // nrow(patient_data): kept patients (time_available > 0, SM:136)
  int status_sum = 0;                // sum(patient_data$event_status) (IL:178)
  std::vector<int> events;           // events_per_interval
  std::vector<double> pt;            // person_time_per_interval
  int events_total = 0;              // sum(events_per_interval) (SM:142)
  double pt_total = 0.0;             // sum(person_time_per_interval) (SM:143, long double)
  double median_fu = 0.0;            // stats::median(observed_time) or 0 (SM:144)
  std::vector<double> obs;           // kept observed times, registry order (scratch)
  std::vector<int> st;               // kept event statuses (scratch)
  std::vector<double> contrib;       // scratch
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
  if (r.n == 0) return;  // SM:118-130: zero vectors, events_total 0L, person_time_total 0, median_followup 0
  for (int i = 0; i < r.n; i++) {
    const double ts = r_pmax2(0.0, cal - r.enroll[i]);                  // pmax(0, calendar_time - enroll_time)
    const double ta = r_pmin2(ts, max_fu);                              // pmin(time_since_enroll, max_follow_up)
    const double ob = r_pmin2(r_pmin2(r.te[i], r.rc[i]), ta);           // pmin(te, rc, time_available)
    const int es = (r.te[i] <= r_pmin2(r.rc[i], ta)) ? 1 : 0;           // as.integer(te <= pmin(rc, time_available))
    if (ta > 0) {                                                       // keep <- time_available > 0
      s.obs.push_back(ob);
      s.st.push_back(es);
    }
  }
  s.n_pat = (int) s.obs.size();
  for (int e : s.st) s.status_sum += e;
  if (s.n_pat == 0) {  // calculate_interval_metrics_fast on 0 rows: rep(0L, K), rep(0, K); median_followup 0
    s.events_total = 0;
    s.pt_total = evt_rsem::sum_ld(s.pt.data(), K);
    return;
  }
  // 1) events per interval: findInterval(event_times, cutpoints, left.open = FALSE, rightmost.closed = FALSE),
  //    drop idx 0, pmin(idx, K), tabulate (SM:23-40)
  const int nb = K + 1;
  for (int i = 0; i < s.n_pat; i++) {
    if (s.st[i] != 1) continue;
    int idx = evt_rsem::find_interval(cut.data(), nb, s.obs[i], false, false, false);
    if (idx == NA_INTEGER) Rcpp::stop("simulate_chunk_cpp: NA event time in the interval tabulation");
    if (idx > 0) {
      if (idx > K) idx = K;
      s.events[idx - 1]++;
    }
  }
  // 2) person-time per interval: at_risk <- obs > lo; contrib <- pmin(obs[at_risk], hi) - lo; contrib[contrib < 0] <- 0;
  //    sum(contrib) in long double (SM:42-60)
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

// ---- posterior: draw_posterior_hazard_samples (PD:29-55) + calculate_median_survival_matrix_cpp -----------------------
// Returns the count of draws with med < thr (lt = true) or med > thr (lt = false), and whether any comparison was NA.
struct Tally { int count; bool any_na; int n; };

arma::vec posterior_medians(const Cfg& c, const Slice& s, const std::vector<double>& lens, int num_samples) {
  arma::ivec ev(c.K);
  for (int k = 0; k < c.K; k++) ev(k) = s.events[k];          // as.integer(events_per_interval)
  arma::vec pt(s.pt);                                         // as.numeric(person_time_per_interval)
  arma::vec a(c.prior_alpha), b(c.prior_beta);
  arma::mat lam = draw_posterior_hazard_samples_cpp(c.K, ev, pt, a, b, num_samples);
  arma::vec il(lens);
  return calculate_median_survival_matrix_cpp(lam, il);
}

Tally tally(const arma::vec& med, double thr, bool lt) {
  Tally t{0, false, (int) med.n_elem};
  for (arma::uword i = 0; i < med.n_elem; i++) {
    const double m = med(i);
    if (std::isnan(m)) { t.any_na = true; continue; }
    if (lt ? (m < thr) : (m > thr)) t.count++;
  }
  return t;
}

// mean() of the logical vector (inventory 4.36): NA if any element is NA, else long double count / n.
double tally_mean(const Tally& t) {
  if (t.any_na) return NA_REAL;
  return evt_rsem::mean_logical(t.count, (R_xlen_t) t.n);
}

// ---- interim_check_hc (IL:171-263) ----------------------------------------------------------------------------------
void interim_hc(const Cfg& c, TrialState& st, double cal, const std::vector<double>& cut, const std::vector<double>& lens,
                Slice& sl) {
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
    // calculate_current_probs_hc (IL:43-66): interval_lengths <- diff(args$interval_cutpoints_sim); pr_eff is not read
    // because efficacy_stopping_rule_hc is FALSE in the supported set (IL:240); computing it draws nothing.
    arma::vec med = posterior_medians(c, sl, lens, c.draws_interim);
    const double pr_fut = tally_mean(tally(med, c.futility_thr_interim[a], true));
    // IL:251-253: !is.null(thr) && is.finite(thr) && pr_fut >= thr; thr is finite in the supported set
    if (ISNAN(pr_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (interim futility probability is NA)");
    if (pr_fut >= c.futility_prob_interim) {
      st.status[a] = STOPPED_FUTILITY;
      st.stop_fut[a] = 1;
      st.final_n_cur[a] = st.enrolled[a];
      st.stop_time[a] = cal;
    }
  }
}

bool any_recruiting(const TrialState& st) {
  for (int s : st.status) if (s == RECRUITING) return true;
  return false;
}

// ---- vs-reference path (step 4) -------------------------------------------------------------------------------------

// calculate_median_survival_piecewise (PH:16-57) on one hazard vector. It draws nothing. Its value is not read on the
// supported (PH) path, but R evaluates it for every sample, and its if() conditions error on NaN ("missing value where
// TRUE/FALSE needed"), so the port evaluates it too and stops where R would stop.
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

struct PhScratch {
  std::vector<double> logHR, shapes, lam, lamT, medC, medT;
};

// sample_vs_ref_medians_ph (PH:105-135, the live definition; PD:241-271 is shadowed, inventory 6.2) on raw slice vectors.
// ph_beta_mode_var (PD:98-118) with EVOLVETRIAL_USE_CPP TRUE: ph_beta_mode_var_cpp(as.integer(E_C), as.numeric(PT_C),
// as.integer(E_T), as.numeric(PT_T), as.numeric(alpha), as.numeric(beta), mu, sigma, tol, as.integer(max_iter)).
// Then rnorm(n, mean, sd): n draws R::rnorm(mean, sd) in index order (R's C_rnorm calls the same rnorm(mu, sigma) per
// element). Then per sample, in order: exp(beta_i); rates <- (beta_prior + PT_C) + exp_beta * PT_T; rgamma(K, shapes,
// rate = rates), which is .Call(C_rgamma, K, shapes, scale = 1 / rates), so K draws R::rgamma(shape_k, 1 / rate_k) in
// interval order; then the two medians (no draws). keep_medians stores them (the unit-test hook); otherwise they are
// evaluated for their error conditions only.
void ph_sample_core(const int* EC, const double* PTC, const int* ET, const double* PTT, const std::vector<double>& alpha,
                    const std::vector<double>& bprior, const std::vector<double>& lens, double mu, double sigma,
                    double tol, int max_iter, int n, bool keep_medians, PhScratch& ph) {
  const int K = (int) alpha.size();
  arma::ivec aEC(K), aET(K);
  arma::vec aPTC(K), aPTT(K);
  for (int k = 0; k < K; k++) {
    aEC(k) = EC[k];
    aET(k) = ET[k];
    aPTC(k) = PTC[k];
    aPTT(k) = PTT[k];
  }
  arma::vec aA(alpha), aB(bprior);
  Rcpp::List bp = ph_beta_mode_var_cpp(aEC, aPTC, aET, aPTT, aA, aB, mu, sigma, tol, max_iter);
  const double bmean = Rcpp::as<double>(bp["mean"]);
  const double bsd = Rcpp::as<double>(bp["sd"]);

  ph.logHR.resize(n);
  for (int i = 0; i < n; i++) ph.logHR[i] = R::rnorm(bmean, bsd);

  ph.shapes.resize(K);
  for (int k = 0; k < K; k++) ph.shapes[k] = (alpha[k] + (double) EC[k]) + (double) ET[k];   // alpha_prior + E_C + E_T
  ph.lam.resize(K);
  ph.lamT.resize(K);
  if (keep_medians) {
    ph.medC.assign(n, 0.0);
    ph.medT.assign(n, 0.0);
  }
  for (int i = 0; i < n; i++) {
    const double eb = std::exp(ph.logHR[i]);
    for (int k = 0; k < K; k++) {
      const double rate = (bprior[k] + PTC[k]) + eb * PTT[k];
      const double scale = 1 / rate;
      ph.lam[k] = R::rgamma(ph.shapes[k], scale);
    }
    for (int k = 0; k < K; k++) ph.lamT[k] = ph.lam[k] * eb;   // lambda * exp_beta
    const double mc = median_piecewise_r(ph.lam.data(), lens);
    const double mt = median_piecewise_r(ph.lamT.data(), lens);
    if (keep_medians) {
      ph.medC[i] = mc;
      ph.medT[i] = mt;
    }
  }
}

void ph_sample(const Cfg& c, const Slice& sC, const Slice& sT, const std::vector<double>& lens, int n, PhScratch& ph) {
  ph_sample_core(sC.events.data(), sC.pt.data(), sT.events.data(), sT.pt.data(), c.prior_alpha, c.prior_beta, lens,
                 c.ph_mu, c.ph_sigma, c.ph_tol, c.ph_max_iter, n, false, ph);
}

// mean(logHR < thr) (lt = true) or mean(logHR >= thr) (lt = false), inventory 4.37: NA if any comparison is NA.
double logical_mean(const std::vector<double>& x, double thr, bool lt) {
  int count = 0;
  for (double v : x) {
    if (std::isnan(v)) return NA_REAL;
    if (lt ? (v < thr) : (v >= thr)) count++;
  }
  return evt_rsem::mean_logical(count, (R_xlen_t) x.size());
}

// gates_pass_for_both_arms (SM:153-258) with the thresholds resolved in R once per frame (cfg$vr_*: {ctrl, trt}).
bool gates_pass_vsref(const Cfg& c, const Slice& sC, const Slice& sT) {
  const double evC = (double) sC.events_total, evT = (double) sT.events_total;   // coalesce_num(events_total, 0)
  const double mfuC = sC.median_fu, mfuT = sT.median_fu;                         // median_followup (0 when no rows)
  const double ptC = sC.pt_total, ptT = sT.pt_total;                             // person_time_total
  const double maxPT_C = c.vr_max_pt[0], maxPT_T = c.vr_max_pt[1];
  const double fracC = (maxPT_C > 0) ? ptC / maxPT_C : 0;
  const double fracT = (maxPT_T > 0) ? ptT / maxPT_T : 0;
  const bool passC = (evC >= c.vr_min_ev[0]) && (mfuC >= c.vr_min_mfu[0]) && (fracC >= c.vr_min_pt[0]);
  const bool passT = (evT >= c.vr_min_ev[1]) && (mfuT >= c.vr_min_mfu[1]) && (fracT >= c.vr_min_pt[1]);
  return passC && passT;
}

// interim_check_vs_ref (IL:84-169) on the supported set: one experimental arm, PH model, no target override, interim
// efficacy stopping off (so pr_eff, IL:30, is computed by R but never read: IL:147 short-circuits on the rule).
void interim_vsref(const Cfg& c, TrialState& st, double cal, const std::vector<double>& cut,
                   const std::vector<double>& lens, Slice& sC, Slice& sT, PhScratch& ph) {
  const int r = c.ref_idx, t = c.trt_idx;
  if (st.status[t] != RECRUITING) return;                                       // IL:96-103
  slice_arm(st.reg[r], cal, c.max_follow_up_sim, cut, c.K, sC);                 // IL:111-112
  // IL:116-119: the control-cache precompute is for the independent model only (not reached, V2)
  slice_arm(st.reg[t], cal, c.max_follow_up_sim, cut, c.K, sT);                 // IL:124-125
  if (!gates_pass_vsref(c, sC, sT)) return;                                     // IL:127-133
  ph_sample(c, sC, sT, lens, c.draws_vsref_interim, ph);                        // IL:3-9, PH:93-135
  const double pr_fut = logical_mean(ph.logHR, c.log_margin_interim, false);    // IL:31
  if (ISNAN(pr_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (interim vs-ref futility probability is NA)");
  if (pr_fut >= c.vr_fut_prob_interim) {                                        // IL:158-166
    st.status[t] = STOPPED_FUTILITY;
    st.stop_fut[t] = 1;
    st.final_n_cur[t] = st.enrolled[t];
    st.stop_time[t] = cal;
  }
}

// the experimental arm is recruiting and below its cap (SD:562-566, SD:641-644, by name; names(mtpa) == arm_names)
bool trt_can_enroll(const Cfg& c, const TrialState& st) {
  const int t = c.trt_idx;
  return st.status[t] == RECRUITING && (double) st.enrolled[t] < c.mtpa[t];
}

}  // namespace

// [[Rcpp::export(name = ".chunk_ph_sample")]]
Rcpp::List chunk_ph_sample(Rcpp::IntegerVector E_C, Rcpp::NumericVector PT_C, Rcpp::IntegerVector E_T,
                           Rcpp::NumericVector PT_T, Rcpp::NumericVector alpha, Rcpp::NumericVector beta,
                           Rcpp::NumericVector interval_lengths, double mu, double sigma, double tol, int max_iter,
                           int num_samples) {
  // Test hook for the PH sampler port (plan P5, P15, P16): the same draws, in the same order, as
  // sample_vs_ref_medians_ph (PH:105-135). Returns list(medCtrl, medTrt, logHR) like the R function.
  const int K = alpha.size();
  if (E_C.size() != K || PT_C.size() != K || E_T.size() != K || PT_T.size() != K || beta.size() != K ||
      interval_lengths.size() != K) {
    Rcpp::stop(".chunk_ph_sample: all vectors must have the same length");
  }
  PhScratch ph;
  std::vector<double> a(alpha.begin(), alpha.end()), b(beta.begin(), beta.end());
  std::vector<double> il(interval_lengths.begin(), interval_lengths.end());
  std::vector<double> ptc(PT_C.begin(), PT_C.end()), ptt(PT_T.begin(), PT_T.end());
  std::vector<int> ec(E_C.begin(), E_C.end()), et(E_T.begin(), E_T.end());
  ph_sample_core(ec.data(), ptc.data(), et.data(), ptt.data(), a, b, il, mu, sigma, tol, max_iter, num_samples, true, ph);
  return Rcpp::List::create(Rcpp::Named("medCtrl") = Rcpp::NumericVector(ph.medC.begin(), ph.medC.end()),
                            Rcpp::Named("medTrt") = Rcpp::NumericVector(ph.medT.begin(), ph.medT.end()),
                            Rcpp::Named("logHR") = Rcpp::NumericVector(ph.logHR.begin(), ph.logHR.end()));
}

// [[Rcpp::export(name = ".chunk_sample_one")]]
int chunk_sample_one(Rcpp::NumericVector prob) {
  // Test hook for the sample() port (plan P2): the 1-based index that sample(seq_along(prob), 1, prob = prob) returns.
  return prob_sample_one(std::vector<double>(prob.begin(), prob.end())) + 1;
}

// [[Rcpp::export(name = ".simulate_chunk_cpp")]]
Rcpp::List simulate_chunk_cpp(Rcpp::IntegerVector sim_indices, Rcpp::List cfg) {
  const Cfg c = unpack(cfg);
  const int A = c.n_arms;
  const int K = c.K;
  const R_xlen_t n_sims = sim_indices.size();

  // chunk accumulators (SD:472-498)
  std::vector<double> sum_n(A, 0.0), sum_ev(A, 0.0), sum_se(A, 0.0), sum_sf(A, 0.0), sum_fe(A, 0.0), sum_ff(A, 0.0),
      sum_fi(A, 0.0);
  double sum_time = 0;
  std::vector<std::vector<double>> raw_n, raw_ev;
  std::vector<double> raw_time;
  if (c.collect_raw_n) raw_n.assign(A, std::vector<double>(n_sims, 0.0));
  if (c.return_variance) {
    raw_ev.assign(A, std::vector<double>(n_sims, 0.0));
    raw_time.assign(n_sims, 0.0);
  }
  // which(sim_indices == sim_idx) (SD:867, 875): every position holding the same index; NA matches nothing.
  bool has_dup = false;
  {
    std::vector<int> v(sim_indices.begin(), sim_indices.end());
    std::sort(v.begin(), v.end());
    for (size_t i = 1; i < v.size(); i++) if (v[i] == v[i - 1]) { has_dup = true; break; }
  }

  TrialState st;
  Slice sl, sl_ref, sl_trt;
  PhScratch ph;
  std::vector<double> cut, lens, ev_times;
  std::vector<double> final_n(A), final_ev(A);
  std::vector<int> final_eff(A), final_fut(A), final_inc(A);
  const double rexp_scale = 1 / c.accrual_rate;   // rexp(1, rate) is .Call(C_rexp, n, 1/rate) (stats::rexp)

  for (R_xlen_t t = 0; t < n_sims; t++) {
    // per-trial reset (SD:522-536, inventory 6.16)
    cut = c.cut_base;
    lens = c.len_base;
    bool rebalance_done = false;            // is.null(rebalance_threshold) is FALSE in the supported set
    double last_enroll_time = 0;
    st.reset(c);
    double current_time = 0.0;
    double next_calendar_look = c.beat;

    auto eligible_any = [&]() {
      for (int a = 0; a < A; a++)
        if (st.status[a] == RECRUITING && (double) st.enrolled[a] < c.mtpa[a]) return true;
      return false;
    };

    // ---- accrual loop (SD:542-679) ----
    while (eligible_any()) {
      const double interarrival = R::rexp(rexp_scale);
      current_time = current_time + interarrival;

      // SD:552-555 (looks_stop_at_final): an arrival at or after max_trial_time enrolls nobody
      if (current_time >= c.max_trial_time) {
        current_time = last_enroll_time;
        break;
      }
      // SD:561-570 (looks_stop_at_final, vs-reference): an arrival that would end the loop through the vs-reference
      // break below enrolls nobody and fires no rebalancing and no look
      if (c.vsref && !trt_can_enroll(c, st)) {
        current_time = last_enroll_time;
        break;
      }

      // SD:572-603: rebalancing check on every arrival until done
      if (!rebalance_done) {
        int total_events = 0;
        ev_times.clear();
        for (int a = 0; a < A; a++) {
          // extract_events_fast (SD:501-520)
          const Registry& r = st.reg[a];
          int count = 0;
          for (int i = 0; i < r.n; i++) {
            const double ts = r_pmax2(0.0, current_time - r.enroll[i]);
            const double ta = r_pmin2(ts, c.max_follow_up_sim);
            if (!(ta > 0)) continue;
            const double ob = r_pmin2(r_pmin2(r.te[i], r.rc[i]), ta);
            if (r.te[i] <= r_pmin2(r.rc[i], ta)) {
              count++;
              ev_times.push_back(ob);
            }
          }
          total_events = total_events + count;
        }
        if ((double) total_events >= c.rebalance_threshold) {
          std::vector<double> new_cuts = evt_rsem::rebalance_cutpoints(ev_times, c.max_follow_up_sim, K);
          if (!new_cuts.empty()) {
            cut = new_cuts;
            for (int k = 0; k < K; k++) lens[k] = cut[k + 1] - cut[k];   // diff(new_cuts)
          }
          rebalance_done = true;
        }
      }

      // SD:628-633: calendar-beat look (pt_targets_abs is NULL in the supported set)
      if (current_time >= next_calendar_look) {
        if (c.vsref) interim_vsref(c, st, current_time, cut, lens, sl_ref, sl_trt, ph);
        else interim_hc(c, st, current_time, cut, lens, sl);
        next_calendar_look = next_calendar_look + c.beat;
      }

      if (!any_recruiting(st)) break;                                   // SD:635
      std::vector<int> elig;                                            // SD:637-638
      for (int a = 0; a < A; a++)
        if (st.status[a] == RECRUITING && (double) st.enrolled[a] < c.mtpa[a]) elig.push_back(a);
      if (elig.empty()) break;
      if (c.vsref && !trt_can_enroll(c, st)) break;                     // SD:640-646 (current_time not reset, 6.12)

      // SD:648-652: probs <- randomization_probs[elig_arms]; probs <- probs / sum(probs); sample(elig_arms, 1, prob = probs)
      std::vector<double> probs(elig.size());
      for (size_t i = 0; i < elig.size(); i++) probs[i] = c.rand_probs[elig[i]];
      const double psum = evt_rsem::sum_ld(probs.data(), (R_xlen_t) probs.size());
      for (size_t i = 0; i < probs.size(); i++) probs[i] = probs[i] / psum;
      const int chosen = elig[prob_sample_one(probs)];

      // SD:655-658: rweibull, then the dropout uniform (drawn, then discarded when random_dropout is FALSE)
      const double t_event_true = R::rweibull(c.shape[chosen], c.scale[chosen]);
      double t_random_censor = R::runif(0.0, c.censor_max_time_sim);
      if (!c.random_dropout) t_random_censor = R_PosInf;

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

    // ---- post-enrollment looks, looks_stop_at_final = TRUE (SD:697-716) ----
    const double t_star = std::min(last_enroll_time + c.min_follow_up_at_final, c.max_trial_time);  // finite branch
    double post_enroll_time = last_enroll_time;
    while (R_FINITE(t_star) && any_recruiting(st)) {
      post_enroll_time = post_enroll_time + c.beat;                    // repeated addition (inventory 4.41, plan P17)
      if (post_enroll_time >= t_star) break;
      if (post_enroll_time <= current_time) continue;
      if (c.vsref) interim_vsref(c, st, post_enroll_time, cut, lens, sl_ref, sl_trt, ph);
      else interim_hc(c, st, post_enroll_time, cut, lens, sl);
    }

    // ---- final analysis time (SD:722-730) ----
    const bool all_stopped_early = !any_recruiting(st);
    bool any_stop_time = false;
    double max_stop = R_NegInf;
    for (int a = 0; a < A; a++) {
      if (!ISNAN(st.stop_time[a])) {
        any_stop_time = true;
        if (st.stop_time[a] > max_stop) max_stop = st.stop_time[a];   // max(stop_time, na.rm = TRUE)
      }
    }
    double final_time;
    if (all_stopped_early && any_stop_time) {
      final_time = max_stop;
    } else {
      final_time = std::min(last_enroll_time + c.min_follow_up_at_final, c.max_trial_time);
    }

    std::fill(final_eff.begin(), final_eff.end(), 0);
    std::fill(final_fut.begin(), final_fut.end(), 0);
    std::fill(final_inc.begin(), final_inc.end(), 0);
    if (c.vsref) {
      // ---- final vs-reference analysis (SD:736-742, 744-756, 794-833), PH branch ----
      slice_arm(st.reg[c.ref_idx], final_time, c.max_follow_up_sim, cut, K, sl_ref);   // ref_slice_final
      for (int a = 0; a < A; a++) {
        slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);            // arm slice (or SD:845-853)
        final_ev[a] = (double) sl.events_total;
        if (st.status[a] != RECRUITING) continue;
        if (a == c.ref_idx) {
          final_inc[a] = 1;                                                            // SD:795-796
          continue;
        }
        ph_sample(c, sl_ref, sl, lens, c.draws_final, ph);                            // SD:798-803
        const double p_eff = logical_mean(ph.logHR, c.log_eff_final, true);           // SD:816
        const double p_fut = logical_mean(ph.logHR, c.log_margin_final, false);       // SD:817
        if (ISNAN(p_eff)) Rcpp::stop("missing value where TRUE/FALSE needed (final vs-ref efficacy probability is NA)");
        if (p_eff >= c.vr_eff_prob_final) {
          final_eff[a] = 1;
        } else {
          if (ISNAN(p_fut)) Rcpp::stop("missing value where TRUE/FALSE needed (final vs-ref futility probability is NA)");
          if (p_fut >= c.vr_fut_prob_final) final_fut[a] = 1;
          else final_inc[a] = 1;
        }
      }
    }
    // ---- final HC analysis (SD:744-793) ----
    for (int a = 0; a < A && !c.vsref; a++) {
      if (st.status[a] != RECRUITING) {
        slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);   // SD:845-853, re-sliced for events
        final_ev[a] = (double) sl.events_total;
        continue;
      }
      slice_arm(st.reg[a], final_time, c.max_follow_up_sim, cut, K, sl);
      final_ev[a] = (double) sl.events_total;
      arma::vec med = posterior_medians(c, sl, lens, c.draws_final);
      const double p_eff = tally_mean(tally(med, c.success_thr[a], false));   // is.finite(success_thr) in the set
      if (!ISNAN(p_eff) && p_eff >= c.final_success_prob) {
        final_eff[a] = 1;
      } else {
        const double p_fut = tally_mean(tally(med, c.futility_thr[a], true));
        if (!ISNAN(p_fut) && p_fut >= c.final_futility_prob) {
          final_fut[a] = 1;
        } else {
          final_inc[a] = 1;
        }
      }
    }
    for (int a = 0; a < A; a++) {
      final_n[a] = (st.final_n_cur[a] == NA_INTEGER) ? (double) st.enrolled[a] : (double) st.final_n_cur[a];
    }

    // ---- accumulators in trial order, plain double + (SD:856-863, plan P9) ----
    for (int a = 0; a < A; a++) {
      sum_n[a] = sum_n[a] + final_n[a];
      sum_ev[a] = sum_ev[a] + final_ev[a];
      sum_se[a] = sum_se[a] + (double) st.stop_eff[a];
      sum_sf[a] = sum_sf[a] + (double) st.stop_fut[a];
      sum_fe[a] = sum_fe[a] + (double) final_eff[a];
      sum_ff[a] = sum_ff[a] + (double) final_fut[a];
      sum_fi[a] = sum_fi[a] + (double) final_inc[a];
    }
    sum_time = sum_time + final_time;

    // ---- raw per-trial values at which(sim_indices == sim_idx) (SD:866-880) ----
    const int sim_idx = sim_indices[t];
    if (sim_idx != NA_INTEGER) {
      auto store = [&](R_xlen_t pos) {
        if (c.collect_raw_n) for (int a = 0; a < A; a++) raw_n[a][pos] = final_n[a];
        if (c.return_variance) {
          for (int a = 0; a < A; a++) raw_ev[a][pos] = final_ev[a];
          raw_time[pos] = final_time;
        }
      };
      if (!has_dup) {
        store(t);
      } else {
        for (R_xlen_t j = 0; j < n_sims; j++) if (sim_indices[j] == sim_idx) store(j);
      }
    }
  }

  // ---- return value (SD:885-898, inventory Section 5) ----
  auto named = [&](const std::vector<double>& v) {
    Rcpp::NumericVector out(v.begin(), v.end());
    out.attr("names") = c.arm_names;
    return out;
  };
  auto raw_list = [&](const std::vector<std::vector<double>>& v) -> SEXP {
    Rcpp::List out(A);
    for (int a = 0; a < A; a++) out[a] = Rcpp::NumericVector(v[a].begin(), v[a].end());
    out.attr("names") = c.arm_names;
    return out;
  };
  SEXP final_n_raw = c.collect_raw_n ? raw_list(raw_n) : R_NilValue;
  SEXP final_events_raw = c.return_variance ? raw_list(raw_ev) : R_NilValue;
  SEXP final_time_raw = c.return_variance ? (SEXP) Rcpp::NumericVector(raw_time.begin(), raw_time.end()) : R_NilValue;
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
