// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// =============================================================================
// BOUNDS CHECKING MACROS
// =============================================================================

#define CHECK_VEC_BOUNDS(vec, idx, func_name) \
  if ((idx) < 0 || (idx) >= static_cast<int>((vec).n_elem)) { \
    Rcpp::stop("[%s] Index %d out of bounds for vector of size %d", \
               func_name, static_cast<int>(idx), static_cast<int>((vec).n_elem)); \
  }

#define CHECK_POSITIVE(val, name, func_name) \
  if ((val) <= 0) { \
    Rcpp::stop("[%s] %s must be positive, got %g", func_name, name, static_cast<double>(val)); \
  }

// =============================================================================
// CONSTANTS
// =============================================================================

enum HybridState {
  STATE_SINGLE = 0,
  STATE_CONSIDER_CONVERSION = 1,
  STATE_BETWEEN = 2,
  STATE_STOP = 3
};

// Trial modes for different design types
enum TrialMode {
  MODE_HYBRID = 0,           // Default: SA → conversion → BA
  MODE_SINGLE_ARM = 1,       // SA only, vs historical control
  MODE_BETWEEN_ARM = 2,      // BA only, randomized comparison
  MODE_DUAL_SINGLE_ARM = 3   // Two independent SAs (both arms reported)
};

// Helper function to parse trial mode from string
TrialMode parse_trial_mode(const std::string& mode_str) {
  if (mode_str == "single_arm") return MODE_SINGLE_ARM;
  if (mode_str == "between_arm") return MODE_BETWEEN_ARM;
  if (mode_str == "dual_single_arm") return MODE_DUAL_SINGLE_ARM;
  return MODE_HYBRID;  // Default
}

// External declarations for PP functions from predictive_probability.cpp
double compute_pp_efficacy_sa_cpp(
    const arma::vec& a_arm,
    const arma::vec& b_arm,
    const arma::vec& hist_hazard,
    double hr_threshold,
    int n_add,
    const arma::vec& interval_cutpoints,
    double accrual_rate,
    double followup,
    double eff_threshold,
    int n_outer,
    bool use_antithetic
);

double compute_pp_futility_sa_cpp(
    const arma::vec& a_arm,
    const arma::vec& b_arm,
    const arma::vec& hist_hazard,
    double hr_threshold,
    int n_add,
    const arma::vec& interval_cutpoints,
    double accrual_rate,
    double followup,
    double fut_threshold,
    int n_outer,
    bool use_antithetic
);

// =============================================================================
// HELPER STRUCTURES
// =============================================================================

struct PatientRecord {
  int arm_idx;
  double enrollment_time;
  double event_time;      // Time from enrollment to event (or Inf if censored)
  bool had_event;
  int phase;              // 1 = SA, 2 = BA
};

struct TrialState {
  // State machine
  HybridState current_state;
  double current_time;
  int current_look;

  // Arms (0 = reference, 1+ = experimental)
  int n_arms;
  int n_intervals;
  std::vector<bool> arm_active;
  int reference_arm_idx;

  // Enrollment
  std::vector<int> n_enrolled;
  std::vector<int> n_enrolled_phase1;

  // Patient registry
  std::vector<PatientRecord> patients;

  // Posteriors per arm (n_arms x n_intervals)
  std::vector<arma::vec> posterior_a;
  std::vector<arma::vec> posterior_b;

  // SA phase decisions
  std::vector<bool> sa_efficacy_reached;
  std::vector<bool> sa_futility_reached;
  std::vector<double> sa_posterior_prob;

  // Conversion
  std::string conversion_decision;  // "GO", "NO_GO", "AMBIGUOUS", "PENDING"
  double pp_at_conversion;
  bool conversion_evaluated;  // Flag to ensure PP is only evaluated once

  // BA phase decisions
  std::vector<bool> ba_efficacy_reached;
  std::vector<bool> ba_futility_reached;
  std::vector<double> ba_posterior_prob;

  // Outcome
  std::string trial_outcome;
  std::string stop_reason;
};

// =============================================================================
// FORWARD DECLARATIONS
// =============================================================================

arma::vec simulate_pwe_survival_batch_internal(int n, const arma::vec& lambda,
                                                const arma::vec& interval_cutpoints);
double compute_pwe_median_internal(const arma::vec& lambda,
                                    const arma::vec& interval_cutpoints);
double compute_ba_posterior_internal(const arma::vec& a_exp, const arma::vec& b_exp,
                                      const arma::vec& a_ref, const arma::vec& b_ref,
                                      const arma::vec& interval_cutpoints,
                                      int n_samples);
double compute_p_single_arm_internal(const arma::vec& post_a, const arma::vec& post_b,
                                      const arma::vec& hist_hazard, double hr_threshold,
                                      const arma::vec& interval_cutpoints,
                                      int n_samples);
double compute_pp_predictive_internal(const arma::vec& a_exp, const arma::vec& b_exp,
                                       const arma::vec& a_ref, const arma::vec& b_ref,
                                       int n_add, const arma::vec& interval_cutpoints,
                                       double accrual_rate, double followup,
                                       double eff_ba, double pp_go, double pp_nogo,
                                       int n_outer);

// =============================================================================
// SIMULATION HELPERS
// =============================================================================

// Simulate PWE survival times (vectorized)
arma::vec simulate_pwe_survival_batch_internal(int n, const arma::vec& lambda,
                                                const arma::vec& interval_cutpoints) {
  int n_intervals = lambda.n_elem;

  // BOUNDS CHECK: interval_cutpoints must have n_intervals + 1 elements
  if (static_cast<int>(interval_cutpoints.n_elem) != n_intervals + 1) {
    Rcpp::stop("[simulate_pwe] interval_cutpoints size (%d) must equal lambda size + 1 (%d)",
               static_cast<int>(interval_cutpoints.n_elem), n_intervals + 1);
  }

  // Handle all zero hazards
  bool all_zero = true;
  for (int j = 0; j < n_intervals; j++) {
    if (lambda(j) > 0) { all_zero = false; break; }
  }
  if (all_zero) {
    arma::vec result(n); result.fill(R_PosInf);
    return result;
  }

  // Compute cumulative hazards
  arma::vec int_widths(n_intervals);
  for (int j = 0; j < n_intervals; j++) {
    int_widths(j) = interval_cutpoints(j + 1) - interval_cutpoints(j);
  }
  arma::vec lambda_safe = arma::clamp(lambda, 0.0, R_PosInf);
  arma::vec interval_hazards = lambda_safe % int_widths;
  arma::vec cum_haz(n_intervals + 1);
  cum_haz(0) = 0.0;
  for (int j = 0; j < n_intervals; j++) {
    cum_haz(j + 1) = cum_haz(j) + interval_hazards(j);
  }

  arma::vec surv_times(n);
  for (int i = 0; i < n; i++) {
    double u = R::runif(0.0, 1.0);
    double target = -std::log(u);

    int idx = 0;
    for (int j = 0; j < n_intervals; j++) {
      if (target >= cum_haz(j)) idx = j;
    }
    if (idx >= n_intervals) idx = n_intervals - 1;

    if (target > cum_haz(n_intervals)) {
      if (lambda_safe(n_intervals - 1) <= 0) {
        surv_times(i) = R_PosInf;
      } else {
        double rem = target - cum_haz(n_intervals);
        surv_times(i) = interval_cutpoints(n_intervals) + rem / lambda_safe(n_intervals - 1);
      }
    } else if (lambda_safe(idx) <= 0) {
      surv_times(i) = interval_cutpoints(idx);
    } else {
      double rem = target - cum_haz(idx);
      surv_times(i) = interval_cutpoints(idx) + rem / lambda_safe(idx);
    }
  }
  return surv_times;
}

// Compute interval metrics (events and exposure)
void compute_interval_metrics(const std::vector<PatientRecord>& patients,
                              double analysis_time,
                              const arma::vec& interval_cutpoints,
                              int arm_idx,
                              arma::vec& events_out,
                              arma::vec& exposure_out) {
  // BOUNDS CHECK: Need at least 2 cutpoints for 1 interval
  if (interval_cutpoints.n_elem < 2) {
    Rcpp::stop("[compute_interval_metrics] Need at least 2 interval cutpoints, got %d",
               static_cast<int>(interval_cutpoints.n_elem));
  }

  int n_intervals = interval_cutpoints.n_elem - 1;
  events_out.zeros(n_intervals);
  exposure_out.zeros(n_intervals);

  for (const auto& p : patients) {
    if (p.arm_idx != arm_idx) continue;

    double time_on_study = analysis_time - p.enrollment_time;
    if (time_on_study <= 0) continue;

    double obs_time = p.had_event ? std::min(p.event_time, time_on_study) : time_on_study;

    for (int j = 0; j < n_intervals; j++) {
      double int_start = interval_cutpoints(j);
      double int_end = interval_cutpoints(j + 1);

      // Exposure
      double exp_end = std::min(obs_time, int_end);
      double exp_start = int_start;
      if (exp_end > exp_start) {
        exposure_out(j) += (exp_end - exp_start);
      }

      // Events
      if (p.had_event && p.event_time >= int_start && p.event_time < int_end &&
          p.event_time <= time_on_study) {
        events_out(j) += 1.0;
      }
    }
  }
}

// Compute median survival from PWE hazards (matches R compute_pwe_median)
// Returns time t where S(t) = 0.5, i.e., cumulative hazard = log(2)
double compute_pwe_median_internal(const arma::vec& lambda,
                                    const arma::vec& interval_cutpoints) {
  int n_intervals = lambda.n_elem;

  // Handle all zero/negative hazards
  bool all_zero = true;
  for (int j = 0; j < n_intervals; j++) {
    if (lambda(j) > 0) {
      all_zero = false;
      break;
    }
  }
  if (all_zero) return R_PosInf;

  double cum_haz = 0.0;
  for (int j = 0; j < n_intervals; j++) {
    double int_start = interval_cutpoints(j);
    double int_end = interval_cutpoints(j + 1);
    double int_width = int_end - int_start;

    // Skip zero hazard intervals
    if (lambda(j) <= 0) continue;

    // Cumulative hazard at end of this interval
    double cum_haz_end = cum_haz + lambda(j) * int_width;

    // Check if median is in this interval (cum_haz = log(2))
    if (cum_haz_end >= M_LN2) {
      double remaining = M_LN2 - cum_haz;
      return int_start + remaining / lambda(j);
    }

    cum_haz = cum_haz_end;
  }

  // Median beyond last interval - extrapolate if last interval has positive hazard
  int last_idx = n_intervals - 1;
  if (lambda(last_idx) <= 0) return R_PosInf;

  double remaining = M_LN2 - cum_haz;
  return interval_cutpoints(n_intervals) + remaining / lambda(last_idx);
}

// Compute P(HR < 1) using median-based Monte Carlo (matches R compute_p_between_arm)
// For PWE models: samples per-interval hazards, computes medians, compares
double compute_ba_posterior_internal(const arma::vec& a_exp, const arma::vec& b_exp,
                                      const arma::vec& a_ref, const arma::vec& b_ref,
                                      const arma::vec& interval_cutpoints,
                                      int n_samples = 2000) {
  int n_intervals = a_exp.n_elem;

  // Guard against invalid parameters
  for (int j = 0; j < n_intervals; j++) {
    if (a_exp(j) <= 0 || b_exp(j) <= 0 || a_ref(j) <= 0 || b_ref(j) <= 0) {
      return NA_REAL;
    }
  }

  int count_exp_better = 0;

  for (int i = 0; i < n_samples; i++) {
    // Sample interval-specific hazards from posteriors
    arma::vec lambda_exp(n_intervals);
    arma::vec lambda_ref(n_intervals);

    for (int j = 0; j < n_intervals; j++) {
      lambda_exp(j) = R::rgamma(a_exp(j), 1.0 / b_exp(j));
      lambda_ref(j) = R::rgamma(a_ref(j), 1.0 / b_ref(j));
      // Clamp to avoid numerical issues
      if (lambda_exp(j) < 1e-10) lambda_exp(j) = 1e-10;
      if (lambda_ref(j) < 1e-10) lambda_ref(j) = 1e-10;
    }

    // Compute median survival for each arm
    double median_exp = compute_pwe_median_internal(lambda_exp, interval_cutpoints);
    double median_ref = compute_pwe_median_internal(lambda_ref, interval_cutpoints);

    // Experimental better if longer median survival (lower hazard)
    if (!R_IsNA(median_exp) && !R_IsNA(median_ref) && median_exp > median_ref) {
      count_exp_better++;
    }
  }

  return (double)count_exp_better / n_samples;
}

// Overload for backward compatibility (aggregated F-test for exponential)
double compute_ba_posterior_internal_aggregated(const arma::vec& a_exp, const arma::vec& b_exp,
                                                 const arma::vec& a_ref, const arma::vec& b_ref) {
  double a_exp_sum = arma::sum(a_exp);
  double b_exp_sum = arma::sum(b_exp);
  double a_ref_sum = arma::sum(a_ref);
  double b_ref_sum = arma::sum(b_ref);

  if (a_exp_sum <= 0 || a_ref_sum <= 0 || b_exp_sum <= 0 || b_ref_sum <= 0) {
    return NA_REAL;
  }

  double q = (b_exp_sum / b_ref_sum) * (a_ref_sum / a_exp_sum);
  return R::pf(q, 2.0 * a_exp_sum, 2.0 * a_ref_sum, 1, 0);
}

// Compute P(HR < threshold) for single-arm vs historical using median-based MC
// (matches R compute_p_single_arm for PWE models)
double compute_p_single_arm_internal(const arma::vec& post_a, const arma::vec& post_b,
                                      const arma::vec& hist_hazard, double hr_threshold,
                                      const arma::vec& interval_cutpoints,
                                      int n_samples = 2000) {
  int n_intervals = post_a.n_elem;

  // For single interval (exponential), use closed-form Gamma CDF
  if (n_intervals == 1) {
    if (post_a(0) <= 0 || post_b(0) <= 0 || hist_hazard(0) <= 0) return NA_REAL;
    double threshold = hr_threshold * hist_hazard(0);
    return R::pgamma(threshold, post_a(0), 1.0 / post_b(0), 1, 0);
  }

  // For PWE: Monte Carlo sampling (matches R implementation)
  // Guard against invalid parameters
  for (int j = 0; j < n_intervals; j++) {
    if (post_a(j) <= 0 || post_b(j) <= 0) {
      return NA_REAL;
    }
  }

  // Historical median - use full PWE computation for multi-interval models
  // For single interval, this reduces to log(2) / hist_hazard[0]
  double hist_median = compute_pwe_median_internal(hist_hazard, interval_cutpoints);

  int count_better = 0;

  for (int i = 0; i < n_samples; i++) {
    // Sample interval-specific hazards from posterior
    arma::vec lambda_samples(n_intervals);
    for (int j = 0; j < n_intervals; j++) {
      lambda_samples(j) = R::rgamma(post_a(j), 1.0 / post_b(j));
      if (lambda_samples(j) < 1e-10) lambda_samples(j) = 1e-10;
    }

    // Compute median survival for this sample
    double median_arm = compute_pwe_median_internal(lambda_samples, interval_cutpoints);

    // HR < c means median_arm > hist_median / c (longer survival is better)
    if (!R_IsNA(median_arm) && !std::isinf(median_arm) &&
        median_arm > hist_median / hr_threshold) {
      count_better++;
    }
  }

  return (double)count_better / n_samples;
}

// Simulate future arm data and compute PP (simplified version)
List simulate_future_arm_internal(int n_patients, const arma::vec& lambda,
                                   const arma::vec& interval_cutpoints,
                                   double accrual_rate, double followup) {
  int n_intervals = lambda.n_elem;

  if (n_patients <= 0) {
    return List::create(
      Named("events") = arma::vec(n_intervals, fill::zeros),
      Named("exposure") = arma::vec(n_intervals, fill::zeros)
    );
  }

  double enrollment_duration = (double)n_patients / accrual_rate;
  arma::vec enrollment_times(n_patients);
  for (int i = 0; i < n_patients; i++) {
    enrollment_times(i) = R::runif(0.0, enrollment_duration);
  }
  enrollment_times = sort(enrollment_times);

  double analysis_time = enrollment_duration + followup;
  arma::vec surv_times = simulate_pwe_survival_batch_internal(n_patients, lambda, interval_cutpoints);
  arma::vec time_available = analysis_time - enrollment_times;
  arma::vec observed_times = arma::min(surv_times, time_available);
  arma::uvec had_event = surv_times <= time_available;

  arma::vec events(n_intervals, fill::zeros);
  arma::vec exposure(n_intervals, fill::zeros);

  for (int j = 0; j < n_intervals; j++) {
    double int_start = interval_cutpoints(j);
    double int_end = interval_cutpoints(j + 1);

    for (int i = 0; i < n_patients; i++) {
      double exp_end = std::min(observed_times(i), int_end);
      if (exp_end > int_start) {
        exposure(j) += (exp_end - int_start);
      }
      if (had_event(i) && surv_times(i) >= int_start && surv_times(i) < int_end) {
        events(j) += 1.0;
      }
    }
  }

  return List::create(Named("events") = events, Named("exposure") = exposure);
}

// Compute predictive probability (Monte Carlo)
double compute_pp_predictive_internal(const arma::vec& a_exp, const arma::vec& b_exp,
                                       const arma::vec& a_ref, const arma::vec& b_ref,
                                       int n_add, const arma::vec& interval_cutpoints,
                                       double accrual_rate, double followup,
                                       double eff_ba, double pp_go, double pp_nogo,
                                       int n_outer) {
  int n_intervals = a_exp.n_elem;
  int success_count = 0;
  int checked = 0;

  // Use antithetic variates
  int n_pairs = (n_outer + 1) / 2;

  for (int i = 0; i < n_pairs; i++) {
    arma::vec u_exp(n_intervals), u_ref(n_intervals);
    for (int j = 0; j < n_intervals; j++) {
      u_exp(j) = R::runif(0.0, 1.0);
      u_ref(j) = R::runif(0.0, 1.0);
    }

    for (int anti = 0; anti < 2; anti++) {
      arma::vec lambda_exp(n_intervals), lambda_ref(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        double ue = (anti == 0) ? u_exp(j) : (1.0 - u_exp(j));
        double ur = (anti == 0) ? u_ref(j) : (1.0 - u_ref(j));
        lambda_exp(j) = std::max(1e-6, R::qgamma(ue, a_exp(j), 1.0/b_exp(j), 1, 0));
        lambda_ref(j) = std::max(1e-6, R::qgamma(ur, a_ref(j), 1.0/b_ref(j), 1, 0));
      }

      List fut_exp = simulate_future_arm_internal(n_add, lambda_exp, interval_cutpoints, accrual_rate, followup);
      List fut_ref = simulate_future_arm_internal(n_add, lambda_ref, interval_cutpoints, accrual_rate, followup);

      arma::vec a_exp_final = a_exp + as<arma::vec>(fut_exp["events"]);
      arma::vec b_exp_final = b_exp + as<arma::vec>(fut_exp["exposure"]);
      arma::vec a_ref_final = a_ref + as<arma::vec>(fut_ref["events"]);
      arma::vec b_ref_final = b_ref + as<arma::vec>(fut_ref["exposure"]);

      double p_between = compute_ba_posterior_internal(a_exp_final, b_exp_final, a_ref_final, b_ref_final, interval_cutpoints, 2000);
      if (!ISNA(p_between) && p_between >= eff_ba) success_count++;
      checked++;
    }

    // Early stopping
    if (checked >= 100 && (checked % 50) < 2) {
      double pp = (double)success_count / checked;
      double se = std::sqrt(pp * (1.0 - pp) / checked);
      if (pp - 2*se > pp_go || pp + 2*se < pp_nogo) break;
    }
  }

  return (double)success_count / checked;
}

// =============================================================================
// MAIN SIMULATION FUNCTION
// =============================================================================

//' Simulate a single hybrid trial (C++ implementation)
//'
//' @param theta_list List of design parameters
//' @param base_args_list List of base trial configuration
//' @param scenario_params_list List of scenario parameters
//' @return List with trial results
//' @export
// [[Rcpp::export]]
List simulate_hybrid_trial_cpp(List theta_list, List base_args_list, List scenario_params_list) {

  // Extract parameters
  double eff_sa = as<double>(theta_list["eff_sa"]);
  double fut_sa = as<double>(theta_list["fut_sa"]);
  double eff_ba = as<double>(theta_list["eff_ba"]);
  double fut_ba = as<double>(theta_list["fut_ba"]);
  int nmax_sa = as<int>(theta_list["nmax_sa"]);
  int nmax_ba = as<int>(theta_list["nmax_ba"]);
  double pp_go = as<double>(theta_list["pp_go"]);
  double pp_nogo = as<double>(theta_list["pp_nogo"]);
  int ev_sa = theta_list.containsElementNamed("ev_sa") ? as<int>(theta_list["ev_sa"]) : 10;
  // FIX (2025-12-24): Add ev_ba extraction for BA phase event gate
  int ev_ba = theta_list.containsElementNamed("ev_ba") ? as<int>(theta_list["ev_ba"]) : 15;
  double hr_threshold_sa = theta_list.containsElementNamed("hr_threshold_sa") ?
    as<double>(theta_list["hr_threshold_sa"]) : 1.0;
  int n_outer = theta_list.containsElementNamed("n_outer") ? as<int>(theta_list["n_outer"]) : 200;

  // Parse trial mode and decision methods
  std::string mode_str = theta_list.containsElementNamed("trial_mode") ?
    as<std::string>(theta_list["trial_mode"]) : "hybrid";
  TrialMode trial_mode = parse_trial_mode(mode_str);

  // Parse efficacy/futility decision methods (posterior vs predictive)
  std::string eff_method = theta_list.containsElementNamed("efficacy_method") ?
    as<std::string>(theta_list["efficacy_method"]) : "posterior";
  std::string fut_method = theta_list.containsElementNamed("futility_method") ?
    as<std::string>(theta_list["futility_method"]) : "posterior";
  bool use_pp_efficacy = (eff_method == "predictive");
  bool use_pp_futility = (fut_method == "predictive");

  // Parse predictive probability futility threshold (default 0.5)
  double pp_futility_threshold = theta_list.containsElementNamed("pp_futility_threshold") ?
    as<double>(theta_list["pp_futility_threshold"]) : 0.5;

  // Parse futility action: "drop_arm" (default), "stop_trial", "continue"
  std::string futility_action = theta_list.containsElementNamed("futility_action") ?
    as<std::string>(theta_list["futility_action"]) : "drop_arm";

  // Parse conversion trigger: "any_single_success" (default), "all_single_success", "k_of_K"
  std::string conversion_trigger = theta_list.containsElementNamed("conversion_trigger") ?
    as<std::string>(theta_list["conversion_trigger"]) : "any_single_success";
  int k_required = theta_list.containsElementNamed("k_required") ?
    as<int>(theta_list["k_required"]) : 1;

  arma::vec interval_cutpoints = as<arma::vec>(base_args_list["interval_cutpoints_sim"]);

  // INPUT VALIDATION: Check interval_cutpoints
  if (interval_cutpoints.n_elem < 2) {
    Rcpp::stop("[simulate_hybrid_trial_cpp] interval_cutpoints must have at least 2 elements");
  }

  int n_intervals = interval_cutpoints.n_elem - 1;
  double accrual_rate = as<double>(base_args_list["overall_accrual_rate"]);

  // INPUT VALIDATION: Check accrual_rate
  if (accrual_rate <= 0) {
    Rcpp::stop("[simulate_hybrid_trial_cpp] accrual_rate must be positive, got %g", accrual_rate);
  }
  double followup = base_args_list.containsElementNamed("max_follow_up_sim") ?
    as<double>(base_args_list["max_follow_up_sim"]) : 24.0;
  double look_interval = base_args_list.containsElementNamed("look_interval") ?
    as<double>(base_args_list["look_interval"]) : 3.0;
  double max_trial_time = base_args_list.containsElementNamed("max_trial_time") ?
    as<double>(base_args_list["max_trial_time"]) : 72.0;

  // Get true hazards from scenario
  arma::vec lambda_exp = as<arma::vec>(scenario_params_list["lambda_exp"]);
  arma::vec lambda_ref = as<arma::vec>(scenario_params_list["lambda_ref"]);
  arma::vec lambda_hist = scenario_params_list.containsElementNamed("lambda_hist") ?
    as<arma::vec>(scenario_params_list["lambda_hist"]) : lambda_ref;

  // INPUT VALIDATION: Check lambda vector sizes match n_intervals
  if (static_cast<int>(lambda_exp.n_elem) != n_intervals) {
    Rcpp::stop("[simulate_hybrid_trial_cpp] lambda_exp size (%d) must equal n_intervals (%d)",
               static_cast<int>(lambda_exp.n_elem), n_intervals);
  }
  if (static_cast<int>(lambda_ref.n_elem) != n_intervals) {
    Rcpp::stop("[simulate_hybrid_trial_cpp] lambda_ref size (%d) must equal n_intervals (%d)",
               static_cast<int>(lambda_ref.n_elem), n_intervals);
  }
  if (static_cast<int>(lambda_hist.n_elem) != n_intervals) {
    Rcpp::stop("[simulate_hybrid_trial_cpp] lambda_hist size (%d) must equal n_intervals (%d)",
               static_cast<int>(lambda_hist.n_elem), n_intervals);
  }

  // Prior parameters - check both naming conventions for compatibility
  arma::vec prior_a;
  if (base_args_list.containsElementNamed("prior_alpha_params_model")) {
    prior_a = as<arma::vec>(base_args_list["prior_alpha_params_model"]);
  } else if (base_args_list.containsElementNamed("prior_a")) {
    prior_a = as<arma::vec>(base_args_list["prior_a"]);
  } else {
    prior_a = arma::vec(n_intervals, fill::ones) * 0.5;
  }

  arma::vec prior_b;
  if (base_args_list.containsElementNamed("prior_beta_params_model")) {
    prior_b = as<arma::vec>(base_args_list["prior_beta_params_model"]);
  } else if (base_args_list.containsElementNamed("prior_b")) {
    prior_b = as<arma::vec>(base_args_list["prior_b"]);
  } else {
    prior_b = arma::vec(n_intervals, fill::ones) * 0.5;  // Use same default as alpha
  }

  // Initialize state
  // NOTE: This implementation assumes exactly 2 arms (1 reference + 1 experimental).
  // Key places that would need modification for N>2 arms:
  //   - state.n_arms initialization (here)
  //   - BA posterior computation (lines ~965-966): compares arms [0] vs [1] only
  //   - Event gate (line ~950): uses ev_ba * 2, should be ev_ba * n_arms
  //   - Result extraction (lines ~1027-1067): hardcodes [0] and [1] indices
  //   - Enrollment checks: only tracks n_enrolled[0] and n_enrolled[1]
  // See docs/TWO_ARM_ASSUMPTIONS.md for full audit.
  TrialState state;
  state.current_time = 0.0;
  state.current_look = 0;
  state.n_arms = 2;  // LIMITATION: Hardcoded to 2 arms (reference + 1 experimental)
  state.n_intervals = n_intervals;
  state.arm_active = {true, true};
  state.reference_arm_idx = 0;  // LIMITATION: Reference arm always index 0
  state.n_enrolled = {0, 0};
  state.n_enrolled_phase1 = {0, 0};
  state.posterior_a = {prior_a, prior_a};
  state.posterior_b = {prior_b, prior_b};
  state.sa_efficacy_reached = {false, false};
  state.sa_futility_reached = {false, false};
  state.sa_posterior_prob = {0.0, 0.0};
  state.ba_efficacy_reached = {false, false};
  state.ba_futility_reached = {false, false};
  state.ba_posterior_prob = {0.0, 0.0};
  state.pp_at_conversion = NA_REAL;
  state.conversion_evaluated = false;
  state.trial_outcome = "";
  state.stop_reason = "";

  // Initialize state based on trial mode
  if (trial_mode == MODE_BETWEEN_ARM) {
    // Between-arm mode: Skip SA phase, start directly in BA
    state.current_state = STATE_BETWEEN;
    state.conversion_decision = "GO";  // Implicit conversion
  } else {
    // Hybrid, single_arm, dual_single_arm: Start in SA phase
    state.current_state = STATE_SINGLE;
    state.conversion_decision = "PENDING";
  }

  // Simulation loop
  double enrollment_time = 0.0;
  int patients_per_look = (int)(accrual_rate * look_interval);

  while (state.current_state != STATE_STOP && state.current_time < max_trial_time) {
    state.current_look++;
    state.current_time += look_interval;

    // Enroll patients
    int n_to_enroll = patients_per_look;
    int active_arms = 0;
    for (int a = 0; a < state.n_arms; a++) {
      if (state.arm_active[a]) active_arms++;
    }
    if (active_arms == 0) break;

    int per_arm = n_to_enroll / active_arms;
    int remainder = n_to_enroll % active_arms;  // Distribute remainder to avoid losing patients
    int arm_idx_in_active = 0;
    for (int a = 0; a < state.n_arms; a++) {
      if (!state.arm_active[a]) continue;

      // Add 1 extra to first 'remainder' arms to distribute evenly
      int this_arm_target = per_arm + (arm_idx_in_active < remainder ? 1 : 0);
      arm_idx_in_active++;

      // Check enrollment limits
      int max_n = (state.current_state == STATE_SINGLE) ? nmax_sa : nmax_ba;
      int can_enroll = std::min(this_arm_target, max_n - state.n_enrolled[a]);
      if (can_enroll <= 0) continue;

      // Simulate patient events
      arma::vec& true_lambda = (a == 0) ? lambda_ref : lambda_exp;
      arma::vec event_times = simulate_pwe_survival_batch_internal(can_enroll, true_lambda, interval_cutpoints);

      for (int i = 0; i < can_enroll; i++) {
        PatientRecord p;
        p.arm_idx = a;
        p.enrollment_time = enrollment_time + R::runif(0, look_interval);
        p.event_time = event_times(i);
        p.had_event = !std::isinf(event_times(i));
        p.phase = (state.current_state == STATE_SINGLE) ? 1 : 2;
        state.patients.push_back(p);
        state.n_enrolled[a]++;
        if (state.current_state == STATE_SINGLE) {
          state.n_enrolled_phase1[a]++;
        }
      }
    }
    enrollment_time += look_interval;

    // Update posteriors
    for (int a = 0; a < state.n_arms; a++) {
      if (!state.arm_active[a]) continue;
      arma::vec events, exposure;
      compute_interval_metrics(state.patients, state.current_time, interval_cutpoints, a, events, exposure);
      state.posterior_a[a] = prior_a + events;
      state.posterior_b[a] = prior_b + exposure;
    }

    // State-specific logic
    if (state.current_state == STATE_SINGLE) {
      // Determine which arms to evaluate based on trial mode
      // For dual_single_arm: evaluate both arms (0 and 1)
      // For single_arm and hybrid: evaluate only experimental arm (1)
      // EXCEPTION: For "not_both_futile" trigger, evaluate both arms even in hybrid mode
      int start_arm = (trial_mode == MODE_DUAL_SINGLE_ARM ||
                       conversion_trigger == "not_both_futile") ? 0 : 1;

      for (int arm = start_arm; arm < state.n_arms; arm++) {
        if (!state.arm_active[arm]) continue;
        if (state.sa_efficacy_reached[arm] || state.sa_futility_reached[arm]) continue;

        // Compute posterior probability vs historical (always computed)
        // Uses median-based MC to match R implementation for PWE models
        double p_single = compute_p_single_arm_internal(
          state.posterior_a[arm], state.posterior_b[arm],
          lambda_hist, hr_threshold_sa,
          interval_cutpoints, 2000  // Harmonized with R sample size
        );
        state.sa_posterior_prob[arm] = p_single;

        // Count events for this arm
        arma::vec events, exposure;
        compute_interval_metrics(state.patients, state.current_time, interval_cutpoints, arm, events, exposure);
        int total_events = (int)arma::sum(events);

        if (total_events >= ev_sa) {
          // Compute decision probabilities (posterior or predictive)
          double p_eff = p_single;  // Default to posterior
          double p_fut = p_single;  // Default to posterior

          // Use predictive probability for efficacy if requested
          if (use_pp_efficacy) {
            int n_add = nmax_sa - state.n_enrolled[arm];
            if (n_add > 0) {
              p_eff = compute_pp_efficacy_sa_cpp(
                state.posterior_a[arm], state.posterior_b[arm],
                lambda_hist, hr_threshold_sa, n_add,
                interval_cutpoints, accrual_rate, followup,
                eff_sa, n_outer, true  // use_antithetic
              );
            }
          }

          // Use predictive probability for futility if requested
          if (use_pp_futility) {
            int n_add = nmax_sa - state.n_enrolled[arm];
            if (n_add > 0) {
              p_fut = compute_pp_futility_sa_cpp(
                state.posterior_a[arm], state.posterior_b[arm],
                lambda_hist, hr_threshold_sa, n_add,
                interval_cutpoints, accrual_rate, followup,
                fut_sa, n_outer, true  // use_antithetic
              );
            }
          }

          // Check efficacy (using p_eff, which is posterior or PP)
          if (p_eff >= eff_sa) {
            state.sa_efficacy_reached[arm] = true;

            // Mode-specific behavior after SA efficacy
            if (trial_mode == MODE_SINGLE_ARM) {
              // Single-arm mode: Stop with success (no conversion)
              state.current_state = STATE_STOP;
              state.conversion_decision = "SKIP";
              state.trial_outcome = "sa_efficacy";
              state.stop_reason = "SA efficacy reached (single-arm mode)";
            } else if (trial_mode == MODE_DUAL_SINGLE_ARM) {
              // Dual single-arm: Continue evaluating, don't transition yet
              // Will check if both arms done below
            }
            // Note: For hybrid mode, conversion trigger is checked after the loop
          }
          // Check futility (using p_fut, which is posterior or PP)
          // Note: For PP futility, high p_fut means high probability of eventual futility
          else if (use_pp_futility ? (p_fut >= pp_futility_threshold) : (p_single <= fut_sa)) {
            state.sa_futility_reached[arm] = true;

            // Handle futility action (matches R handle_sa_futility)
            if (futility_action == "stop_trial") {
              state.current_state = STATE_STOP;
              state.trial_outcome = "sa_futility_stop";
              state.stop_reason = "Arm hit SA futility, trial stopped";
            } else if (futility_action == "drop_arm") {
              state.arm_active[arm] = false;
            }
            // "continue" - do nothing, continue monitoring
          }
        }
      }

      // Check conversion trigger for hybrid mode (matches R check_conversion_trigger)
      if (trial_mode == MODE_HYBRID && state.current_state == STATE_SINGLE) {
        int efficacy_count = 0;
        int n_active = 0;
        for (int arm = 1; arm < state.n_arms; arm++) {  // Skip reference arm (0)
          if (state.arm_active[arm]) {
            n_active++;
            if (state.sa_efficacy_reached[arm]) efficacy_count++;
          }
        }

        bool trigger_met = false;
        if (conversion_trigger == "any_single_success") {
          trigger_met = (efficacy_count >= 1);
        } else if (conversion_trigger == "all_single_success") {
          trigger_met = (efficacy_count == n_active && n_active > 0);
        } else if (conversion_trigger == "k_of_K") {
          trigger_met = (efficacy_count >= k_required);
        } else if (conversion_trigger == "not_both_futile") {
          // NEW: "NOT(both SA futile)" conversion rule
          // Wait until BOTH arms have reached SA decision (efficacy, futility, or max)
          // Then convert to BA unless BOTH arms hit SA futility
          bool all_decided = true;
          int futility_count = 0;
          int decided_count = 0;
          for (int arm = 0; arm < state.n_arms; arm++) {  // Include reference arm (0)
            bool arm_decided = state.sa_efficacy_reached[arm] ||
                               state.sa_futility_reached[arm] ||
                               state.n_enrolled[arm] >= nmax_sa;
            if (!arm_decided && state.arm_active[arm]) {
              all_decided = false;
            }
            if (arm_decided) {
              decided_count++;
              if (state.sa_futility_reached[arm]) {
                futility_count++;
              }
            }
          }
          // Only trigger conversion when all arms decided AND NOT both futile
          if (all_decided && decided_count >= 2) {
            bool both_futile = (futility_count == decided_count);
            if (!both_futile) {
              // At least one arm not futile - proceed to BA comparison
              trigger_met = true;
            } else {
              // Both arms futile - stop trial (don't convert)
              state.current_state = STATE_STOP;
              state.conversion_decision = "NO_GO";
              state.trial_outcome = "both_sa_futile";
              state.stop_reason = "Both arms declared SA futile - stopping";
            }
          }
        } else if (conversion_trigger == "exp_not_futile") {
          // NEW: Convert only if Experimental arm (index 1) is NOT futile
          // Ignore Reference arm entirely for conversion decision
          int exp_idx = 1;  // Experimental arm is index 1
          bool exp_decided = state.sa_efficacy_reached[exp_idx] ||
                             state.sa_futility_reached[exp_idx] ||
                             state.n_enrolled[exp_idx] >= nmax_sa;

          if (exp_decided) {
            if (!state.sa_futility_reached[exp_idx]) {
              // Experimental arm not futile - proceed to BA comparison
              trigger_met = true;
            } else {
              // Experimental arm futile - stop trial immediately
              state.current_state = STATE_STOP;
              state.conversion_decision = "NO_GO";
              state.trial_outcome = "exp_sa_futile";
              state.stop_reason = "Experimental arm declared SA futile - stopping";
            }
          }
        }

        if (trigger_met) {
          state.current_state = STATE_CONSIDER_CONVERSION;
        }
      }

      // For dual_single_arm: check if all arms have reached a decision
      if (trial_mode == MODE_DUAL_SINGLE_ARM && state.current_state == STATE_SINGLE) {
        bool all_decided = true;
        for (int arm = 0; arm < state.n_arms; arm++) {
          if (state.arm_active[arm] &&
              !state.sa_efficacy_reached[arm] &&
              !state.sa_futility_reached[arm] &&
              state.n_enrolled[arm] < nmax_sa) {
            all_decided = false;
            break;
          }
        }
        if (all_decided) {
          state.current_state = STATE_STOP;
          state.conversion_decision = "SKIP";
          state.trial_outcome = "dual_sa_complete";
          state.stop_reason = "Dual single-arm phase complete";
        }
      }

      // Check if all arms dropped (only if still in SA phase)
      if (state.current_state == STATE_SINGLE) {
        bool any_active = false;
        int check_start = (trial_mode == MODE_DUAL_SINGLE_ARM ||
                           conversion_trigger == "not_both_futile") ? 0 : 1;
        for (int a = check_start; a < state.n_arms; a++) {
          if (state.arm_active[a]) any_active = true;
        }
        if (!any_active) {
          state.current_state = STATE_STOP;
          state.trial_outcome = "all_arms_futile";
          state.stop_reason = "All arms dropped for futility";
        }
      }

      // Check max N (only if still in SA phase - don't override other transitions)
      if (state.current_state == STATE_SINGLE) {
        bool max_n_reached = true;
        int check_start = (trial_mode == MODE_DUAL_SINGLE_ARM ||
                           conversion_trigger == "not_both_futile") ? 0 : 1;
        for (int a = check_start; a < state.n_arms; a++) {
          if (state.arm_active[a] && state.n_enrolled[a] < nmax_sa) {
            max_n_reached = false;
            break;
          }
        }
        if (max_n_reached) {
          state.current_state = STATE_STOP;
          state.trial_outcome = "max_n_single_phase";
          state.stop_reason = "Max enrollment reached in SA phase";
        }
      }
    }
    else if (state.current_state == STATE_CONSIDER_CONVERSION) {
      // Only evaluate PP once (matching R behavior)
      if (!state.conversion_evaluated) {
        state.conversion_evaluated = true;

        // For "not_both_futile" or "exp_not_futile" triggers, conversion is automatic (no PP check)
        // We already decided to convert when setting STATE_CONSIDER_CONVERSION
        if (conversion_trigger == "not_both_futile" || conversion_trigger == "exp_not_futile") {
          state.conversion_decision = "GO";
          state.current_state = STATE_BETWEEN;
          state.pp_at_conversion = 1.0;  // Marker for automatic conversion
        } else {
          // Standard PP-based conversion decision
          int n_add = nmax_ba - state.n_enrolled[1];
          if (n_add > 0) {
            double pp = compute_pp_predictive_internal(
              state.posterior_a[1], state.posterior_b[1],
              state.posterior_a[0], state.posterior_b[0],
              n_add, interval_cutpoints, accrual_rate, followup,
              eff_ba, pp_go, pp_nogo, n_outer
            );
            state.pp_at_conversion = pp;

            if (pp >= pp_go) {
              // GO: Proceed to between-arm phase
              state.conversion_decision = "GO";
              state.current_state = STATE_BETWEEN;
            } else if (pp < pp_nogo) {
              // NO-GO: PP too low
              state.conversion_decision = "NO_GO";
              state.current_state = STATE_STOP;
              state.trial_outcome = "conversion_nogo";
              state.stop_reason = "PP too low for conversion";
            } else {
              // AMBIGUOUS: PP between pp_nogo and pp_go
              // Match R behavior: immediately stop with AMBIGUOUS_NOGO
              state.conversion_decision = "AMBIGUOUS";
              state.current_state = STATE_STOP;
              state.trial_outcome = "conversion_ambiguous";
              state.stop_reason = "PP in ambiguous region - default to no-go";
            }
          } else {
            state.current_state = STATE_STOP;
            state.trial_outcome = "max_n_reached";
            state.stop_reason = "Max enrollment reached at conversion";
          }
        }
      }
    }
    else if (state.current_state == STATE_BETWEEN) {
      // BA phase - check efficacy/futility
      // FIX (2025-12-24): Add event count check before BA decision (matches R code hybrid_trial.R:422)
      arma::vec events_exp, exp_exp, events_ref, exp_ref;
      compute_interval_metrics(state.patients, state.current_time, interval_cutpoints, 1, events_exp, exp_exp);
      compute_interval_metrics(state.patients, state.current_time, interval_cutpoints, 0, events_ref, exp_ref);
      int total_ba_events = (int)(arma::sum(events_exp) + arma::sum(events_ref));

      // Skip BA decision if insufficient events (ev_ba * 2 = events across both arms)
      // LIMITATION: Hardcoded * 2 for 2 arms; for N arms use ev_ba * state.n_arms
      if (total_ba_events < ev_ba * 2) {
        // Not enough events for reliable BA inference - continue enrollment
        // Check max N only
        if (state.n_enrolled[1] >= nmax_ba) {
          state.current_state = STATE_STOP;
          state.trial_outcome = "max_n_ba";
          state.stop_reason = "Max enrollment reached in BA phase (insufficient events)";
        }
        continue;
      }

      // LIMITATION: BA comparison hardcoded for arms [1] vs [0] (exp vs ref)
      // For N>2 arms, would need pairwise or global comparison strategy
      double p_between = compute_ba_posterior_internal(
        state.posterior_a[1], state.posterior_b[1],
        state.posterior_a[0], state.posterior_b[0],
        interval_cutpoints, 2000  // Harmonized with R sample size
      );
      state.ba_posterior_prob[1] = p_between;  // Only stores result for arm [1]

      if (p_between >= eff_ba) {
        state.ba_efficacy_reached[1] = true;
        state.current_state = STATE_STOP;
        state.trial_outcome = "ba_efficacy";
        state.stop_reason = "BA efficacy reached";
      } else if (p_between <= fut_ba) {
        state.ba_futility_reached[1] = true;
        state.current_state = STATE_STOP;
        state.trial_outcome = "futility_ba";
        state.stop_reason = "BA futility reached";
      }

      // Check max N
      if (state.n_enrolled[1] >= nmax_ba) {
        state.current_state = STATE_STOP;
        if (state.trial_outcome.empty()) {
          state.trial_outcome = "max_n_ba";
          state.stop_reason = "Max enrollment reached in BA phase";
        }
      }
    }
  }

  // Handle case where trial ended with empty outcome
  // This can happen if stuck in CONSIDER_CONVERSION with ambiguous PP
  if (state.trial_outcome.empty()) {
    if (state.current_state == STATE_CONSIDER_CONVERSION) {
      // Ambiguous PP - treat as NO-GO (conservative)
      state.trial_outcome = "ambiguous_no_conversion";
      state.stop_reason = "Max time reached with ambiguous PP";
      state.conversion_decision = "AMBIGUOUS";
    } else if (state.current_state == STATE_SINGLE) {
      state.trial_outcome = "max_time_single_phase";
      state.stop_reason = "Max time reached in SA phase";
    } else if (state.current_state == STATE_BETWEEN) {
      state.trial_outcome = "max_time_ba";
      state.stop_reason = "Max time reached in BA phase";
    } else {
      state.trial_outcome = "unknown";
      state.stop_reason = "Unknown termination";
    }
  }

  // Compile results
  int total_n = 0;
  for (int a = 0; a < state.n_arms; a++) {
    total_n += state.n_enrolled[a];
  }

  // Mode-specific success definition
  bool is_success = false;
  bool is_success_exp = false;  // Experimental arm success
  bool is_success_ref = false;  // Reference arm success (for dual_single_arm)

  if (trial_mode == MODE_SINGLE_ARM) {
    // Single-arm mode: Success = SA efficacy for experimental arm
    is_success_exp = state.sa_efficacy_reached[1];
    is_success = is_success_exp;
  } else if (trial_mode == MODE_BETWEEN_ARM) {
    // Between-arm mode: Success = BA efficacy
    is_success_exp = state.ba_efficacy_reached[1];
    is_success = is_success_exp;
  } else if (trial_mode == MODE_DUAL_SINGLE_ARM) {
    // Dual single-arm: Report both arms separately
    is_success_exp = state.sa_efficacy_reached[1];
    is_success_ref = state.sa_efficacy_reached[0];
    is_success = is_success_exp;  // Default to experimental for backward compat
  } else {
    // Hybrid mode (default): BA efficacy OR SA efficacy with NO-GO
    bool ba_efficacy = (state.trial_outcome == "ba_efficacy");
    bool sa_efficacy_nogo = state.sa_efficacy_reached[1] &&
      (state.conversion_decision == "NO_GO" || state.conversion_decision == "AMBIGUOUS");
    is_success_exp = ba_efficacy || sa_efficacy_nogo;
    is_success = is_success_exp;
  }

  bool converted = (state.conversion_decision == "GO");

  // Convert std::vectors to IntegerVector for R
  IntegerVector n_enrolled_r(state.n_enrolled.begin(), state.n_enrolled.end());
  IntegerVector n_phase1_r(state.n_enrolled_phase1.begin(), state.n_enrolled_phase1.end());

  // Extract bool values explicitly (std::vector<bool> uses bit references)
  bool sa_eff_exp = state.sa_efficacy_reached[1];
  bool sa_fut_exp = state.sa_futility_reached[1];
  bool sa_eff_ref = state.sa_efficacy_reached[0];
  bool sa_fut_ref = state.sa_futility_reached[0];
  bool ba_eff = state.ba_efficacy_reached[1];
  bool ba_fut = state.ba_futility_reached[1];

  return List::create(
    Named("trial_outcome") = state.trial_outcome,
    Named("stop_reason") = state.stop_reason,
    Named("total_n") = total_n,
    Named("n_enrolled") = n_enrolled_r,
    Named("n_enrolled_exp") = state.n_enrolled[1],
    Named("n_enrolled_ref") = state.n_enrolled[0],
    Named("n_phase1") = n_phase1_r,
    Named("converted") = converted,
    Named("pp_at_conversion") = state.pp_at_conversion,
    Named("is_success") = is_success,
    Named("is_success_exp") = is_success_exp,
    Named("is_success_ref") = is_success_ref,
    Named("trial_mode") = mode_str,
    Named("final_look") = state.current_look,
    Named("final_time") = state.current_time,
    Named("sa_efficacy") = sa_eff_exp,
    Named("sa_futility") = sa_fut_exp,
    Named("sa_efficacy_exp") = sa_eff_exp,
    Named("sa_futility_exp") = sa_fut_exp,
    Named("sa_efficacy_ref") = sa_eff_ref,
    Named("sa_futility_ref") = sa_fut_ref,
    Named("ba_efficacy") = ba_eff,
    Named("ba_futility") = ba_fut
  );
}

//' Run multiple hybrid trial simulations (C++ implementation)
//'
//' @param n_sim Number of simulations
//' @param theta_list Design parameters
//' @param base_args_list Base configuration
//' @param scenario_params_list Scenario parameters
//' @return Data frame with simulation results
//' @export
// [[Rcpp::export]]
DataFrame run_hybrid_simulations_cpp(int n_sim, List theta_list,
                                      List base_args_list, List scenario_params_list) {

  std::vector<std::string> outcomes(n_sim);
  std::vector<int> total_n(n_sim);
  std::vector<int> n_enrolled_exp(n_sim), n_enrolled_ref(n_sim);
  std::vector<bool> is_success(n_sim), is_success_exp(n_sim), is_success_ref(n_sim);
  std::vector<bool> converted(n_sim);
  std::vector<double> pp_conversion(n_sim);
  std::vector<bool> sa_eff_exp(n_sim), sa_fut_exp(n_sim);
  std::vector<bool> sa_eff_ref(n_sim), sa_fut_ref(n_sim);
  std::vector<bool> ba_eff(n_sim), ba_fut(n_sim);

  for (int i = 0; i < n_sim; i++) {
    List result = simulate_hybrid_trial_cpp(theta_list, base_args_list, scenario_params_list);

    outcomes[i] = as<std::string>(result["trial_outcome"]);
    total_n[i] = as<int>(result["total_n"]);
    n_enrolled_exp[i] = as<int>(result["n_enrolled_exp"]);
    n_enrolled_ref[i] = as<int>(result["n_enrolled_ref"]);
    is_success[i] = as<bool>(result["is_success"]);
    is_success_exp[i] = as<bool>(result["is_success_exp"]);
    is_success_ref[i] = as<bool>(result["is_success_ref"]);
    converted[i] = as<bool>(result["converted"]);
    pp_conversion[i] = as<double>(result["pp_at_conversion"]);
    sa_eff_exp[i] = as<bool>(result["sa_efficacy_exp"]);
    sa_fut_exp[i] = as<bool>(result["sa_futility_exp"]);
    sa_eff_ref[i] = as<bool>(result["sa_efficacy_ref"]);
    sa_fut_ref[i] = as<bool>(result["sa_futility_ref"]);
    ba_eff[i] = as<bool>(result["ba_efficacy"]);
    ba_fut[i] = as<bool>(result["ba_futility"]);
  }

  return DataFrame::create(
    Named("outcome") = outcomes,
    Named("total_n") = total_n,
    Named("n_enrolled_exp") = n_enrolled_exp,
    Named("n_enrolled_ref") = n_enrolled_ref,
    Named("is_success") = is_success,
    Named("is_success_exp") = is_success_exp,
    Named("is_success_ref") = is_success_ref,
    Named("converted") = converted,
    Named("pp_conversion") = pp_conversion,
    Named("sa_efficacy") = sa_eff_exp,
    Named("sa_futility") = sa_fut_exp,
    Named("sa_efficacy_exp") = sa_eff_exp,
    Named("sa_futility_exp") = sa_fut_exp,
    Named("sa_efficacy_ref") = sa_eff_ref,
    Named("sa_futility_ref") = sa_fut_ref,
    Named("ba_efficacy") = ba_eff,
    Named("ba_futility") = ba_fut
  );
}

//' Compute hybrid trial operating characteristics (C++ implementation)
//'
//' @param n_sim Number of simulations
//' @param theta_list Design parameters
//' @param base_args_list Base configuration
//' @param scenario_params_list Scenario parameters (alternative hypothesis)
//' @param null_scenario_list Scenario parameters (null hypothesis)
//' @return List with power, type1, EN, conversion rate (including per-arm metrics)
//' @export
// [[Rcpp::export]]
List compute_hybrid_oc_cpp(int n_sim, List theta_list,
                            List base_args_list,
                            List scenario_params_list,
                            List null_scenario_list) {

  // Run under alternative (power)
  DataFrame alt_results = run_hybrid_simulations_cpp(n_sim, theta_list, base_args_list, scenario_params_list);
  LogicalVector alt_success = alt_results["is_success"];
  LogicalVector alt_success_exp = alt_results["is_success_exp"];
  LogicalVector alt_success_ref = alt_results["is_success_ref"];
  IntegerVector alt_n = alt_results["total_n"];
  IntegerVector alt_n_exp = alt_results["n_enrolled_exp"];
  IntegerVector alt_n_ref = alt_results["n_enrolled_ref"];
  LogicalVector alt_conv = alt_results["converted"];

  double power = 0.0, power_exp = 0.0, power_ref = 0.0;
  double en_alt = 0.0, en_alt_exp = 0.0, en_alt_ref = 0.0;
  double conv_rate = 0.0;
  for (int i = 0; i < n_sim; i++) {
    if (alt_success[i]) power += 1.0;
    if (alt_success_exp[i]) power_exp += 1.0;
    if (alt_success_ref[i]) power_ref += 1.0;
    en_alt += alt_n[i];
    en_alt_exp += alt_n_exp[i];
    en_alt_ref += alt_n_ref[i];
    if (alt_conv[i]) conv_rate += 1.0;
  }
  power /= n_sim;
  power_exp /= n_sim;
  power_ref /= n_sim;
  en_alt /= n_sim;
  en_alt_exp /= n_sim;
  en_alt_ref /= n_sim;
  conv_rate /= n_sim;

  // Run under null (type I error)
  DataFrame null_results = run_hybrid_simulations_cpp(n_sim, theta_list, base_args_list, null_scenario_list);
  LogicalVector null_success = null_results["is_success"];
  LogicalVector null_success_exp = null_results["is_success_exp"];
  LogicalVector null_success_ref = null_results["is_success_ref"];
  IntegerVector null_n = null_results["total_n"];
  IntegerVector null_n_exp = null_results["n_enrolled_exp"];
  IntegerVector null_n_ref = null_results["n_enrolled_ref"];
  LogicalVector null_conv = null_results["converted"];

  double type1 = 0.0, type1_exp = 0.0, type1_ref = 0.0;
  double en_null = 0.0, en_null_exp = 0.0, en_null_ref = 0.0;
  double conv_rate_null = 0.0;
  for (int i = 0; i < n_sim; i++) {
    if (null_success[i]) type1 += 1.0;
    if (null_success_exp[i]) type1_exp += 1.0;
    if (null_success_ref[i]) type1_ref += 1.0;
    en_null += null_n[i];
    en_null_exp += null_n_exp[i];
    en_null_ref += null_n_ref[i];
    if (null_conv[i]) conv_rate_null += 1.0;
  }
  type1 /= n_sim;
  type1_exp /= n_sim;
  type1_ref /= n_sim;
  en_null /= n_sim;
  en_null_exp /= n_sim;
  en_null_ref /= n_sim;
  conv_rate_null /= n_sim;

  return List::create(
    Named("power") = power,
    Named("power_exp") = power_exp,
    Named("power_ref") = power_ref,
    Named("type1") = type1,
    Named("type1_exp") = type1_exp,
    Named("type1_ref") = type1_ref,
    Named("EN_alt") = en_alt,
    Named("EN_alt_exp") = en_alt_exp,
    Named("EN_alt_ref") = en_alt_ref,
    Named("EN_null") = en_null,
    Named("EN_null_exp") = en_null_exp,
    Named("EN_null_ref") = en_null_ref,
    Named("EN_total") = (en_alt + en_null) / 2.0,
    Named("conversion_rate") = conv_rate,
    Named("conversion_rate_null") = conv_rate_null,
    Named("n_sim") = n_sim
  );
}
