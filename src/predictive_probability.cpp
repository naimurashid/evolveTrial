// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// =============================================================================
// HELPER: Simulate PWE survival times (vectorized batch)
// =============================================================================

// Internal: simulate n survival times from PWE model
arma::vec simulate_pwe_survival_batch_impl(
    int n,
    const arma::vec& lambda,
    const arma::vec& interval_cutpoints
) {
  int n_intervals = lambda.n_elem;

  // Handle edge case: all zero hazards
  bool all_zero = true;
  for (int j = 0; j < n_intervals; j++) {
    if (lambda(j) > 0) {
      all_zero = false;
      break;
    }
  }
  if (all_zero) {
    arma::vec result(n);
    result.fill(R_PosInf);
    return result;
  }

  // Pre-compute interval widths and cumulative hazards
  arma::vec int_widths(n_intervals);
  for (int j = 0; j < n_intervals; j++) {
    int_widths(j) = interval_cutpoints(j + 1) - interval_cutpoints(j);
  }

  // Safe lambda (clamp negatives to zero)
  arma::vec lambda_safe = arma::clamp(lambda, 0.0, R_PosInf);

  // Interval hazards and cumulative hazards at boundaries
  arma::vec interval_hazards = lambda_safe % int_widths;
  arma::vec cum_haz_boundaries(n_intervals + 1);
  cum_haz_boundaries(0) = 0.0;
  for (int j = 0; j < n_intervals; j++) {
    cum_haz_boundaries(j + 1) = cum_haz_boundaries(j) + interval_hazards(j);
  }

  // Generate survival times
  arma::vec surv_times(n);

  for (int i = 0; i < n; i++) {
    // Generate uniform and convert to target cumulative hazard
    double u = R::runif(0.0, 1.0);
    double target_cum_haz = -std::log(u);

    // Find which interval this falls in
    int interval_idx = 0;
    for (int j = 0; j < n_intervals; j++) {
      if (target_cum_haz >= cum_haz_boundaries(j)) {
        interval_idx = j;
      }
    }

    // Clamp to valid range
    if (interval_idx >= n_intervals) interval_idx = n_intervals - 1;

    // Check if beyond last interval
    if (target_cum_haz > cum_haz_boundaries(n_intervals)) {
      // Event in tail (after last cutpoint)
      if (lambda_safe(n_intervals - 1) <= 0) {
        surv_times(i) = R_PosInf;
      } else {
        double remaining = target_cum_haz - cum_haz_boundaries(n_intervals);
        surv_times(i) = interval_cutpoints(n_intervals) +
                        remaining / lambda_safe(n_intervals - 1);
      }
    } else if (lambda_safe(interval_idx) <= 0) {
      // Zero hazard interval - event at start of interval
      surv_times(i) = interval_cutpoints(interval_idx);
    } else {
      // Normal case: compute time within interval
      double remaining_haz = target_cum_haz - cum_haz_boundaries(interval_idx);
      surv_times(i) = interval_cutpoints(interval_idx) +
                      remaining_haz / lambda_safe(interval_idx);
    }
  }

  return surv_times;
}

//' Simulate PWE survival times (C++ vectorized)
//'
//' @param n Number of survival times to simulate
//' @param lambda Hazard rates per interval
//' @param interval_cutpoints Interval boundaries (length = n_intervals + 1)
//' @return Vector of n survival times
//' @export
// [[Rcpp::export]]
arma::vec simulate_pwe_survival_batch_cpp(
    int n,
    const arma::vec& lambda,
    const arma::vec& interval_cutpoints
) {
  return simulate_pwe_survival_batch_impl(n, lambda, interval_cutpoints);
}

// =============================================================================
// HELPER: Simulate future arm data
// =============================================================================

// Internal: simulate future events and exposure for one arm
// Returns list with events and exposure vectors
List simulate_future_arm_pwe_impl(
    int n_patients,
    const arma::vec& lambda,
    const arma::vec& interval_cutpoints,
    double accrual_rate,
    double followup
) {
  int n_intervals = lambda.n_elem;

  // Guard for zero patients
  if (n_patients <= 0) {
    arma::vec events(n_intervals, fill::zeros);
    arma::vec exposure(n_intervals, fill::zeros);
    return List::create(
      Named("events") = events,
      Named("exposure") = exposure
    );
  }

  // Enrollment times (uniform over enrollment period)
  double enrollment_duration = (double)n_patients / accrual_rate;
  arma::vec enrollment_times(n_patients);
  for (int i = 0; i < n_patients; i++) {
    enrollment_times(i) = R::runif(0.0, enrollment_duration);
  }
  enrollment_times = sort(enrollment_times);

  // Analysis time
  double analysis_time = enrollment_duration + followup;

  // Simulate survival times
  arma::vec surv_times = simulate_pwe_survival_batch_impl(n_patients, lambda, interval_cutpoints);

  // Compute observed times and events
  arma::vec time_available = analysis_time - enrollment_times;
  arma::vec observed_times = arma::min(surv_times, time_available);
  arma::uvec had_event = surv_times <= time_available;

  // Compute interval metrics
  arma::vec events(n_intervals, fill::zeros);
  arma::vec exposure(n_intervals, fill::zeros);

  for (int j = 0; j < n_intervals; j++) {
    double int_start = interval_cutpoints(j);
    double int_end = interval_cutpoints(j + 1);

    for (int i = 0; i < n_patients; i++) {
      // Exposure in this interval
      double exp_start = std::max(0.0, int_start);
      double exp_end = std::min(observed_times(i), int_end);
      double interval_exp = std::max(0.0, exp_end - exp_start);
      exposure(j) += interval_exp;

      // Event in this interval
      if (had_event(i) && surv_times(i) >= int_start && surv_times(i) < int_end) {
        events(j) += 1.0;
      }
    }
  }

  return List::create(
    Named("events") = events,
    Named("exposure") = exposure
  );
}

//' Simulate future arm data under PWE model (C++)
//'
//' @param n_patients Number of future patients
//' @param lambda Hazard rates per interval
//' @param interval_cutpoints Interval boundaries
//' @param accrual_rate Patients per time unit
//' @param followup Follow-up time after enrollment complete
//' @return List with events and exposure vectors
//' @export
// [[Rcpp::export]]
List simulate_future_arm_pwe_cpp(
    int n_patients,
    const arma::vec& lambda,
    const arma::vec& interval_cutpoints,
    double accrual_rate,
    double followup
) {
  return simulate_future_arm_pwe_impl(n_patients, lambda, interval_cutpoints,
                                       accrual_rate, followup);
}

// =============================================================================
// HELPER: Compute BA posterior probability P(HR < 1)
// =============================================================================

// Uses F-distribution on aggregated parameters
double compute_ba_posterior_impl(
    const arma::vec& a_exp,
    const arma::vec& b_exp,
    const arma::vec& a_ref,
    const arma::vec& b_ref
) {
  // Aggregate over intervals
  double a_exp_sum = arma::sum(a_exp);
  double b_exp_sum = arma::sum(b_exp);
  double a_ref_sum = arma::sum(a_ref);
  double b_ref_sum = arma::sum(b_ref);

  // Guard against zero values
  if (a_exp_sum <= 0 || a_ref_sum <= 0 || b_exp_sum <= 0 || b_ref_sum <= 0) {
    return NA_REAL;
  }

  // F-distribution for ratio of Gamma means
  double q = (b_exp_sum / b_ref_sum) * (a_ref_sum / a_exp_sum);
  double df1 = 2.0 * a_exp_sum;
  double df2 = 2.0 * a_ref_sum;

  return R::pf(q, df1, df2, 1, 0);  // lower.tail=TRUE, log.p=FALSE
}

//' Compute between-arm posterior probability (C++)
//'
//' @param a_exp Posterior shape for experimental arm
//' @param b_exp Posterior rate for experimental arm
//' @param a_ref Posterior shape for reference arm
//' @param b_ref Posterior rate for reference arm
//' @return P(HR < 1 | data)
//' @export
// [[Rcpp::export]]
double compute_ba_posterior_cpp(
    const arma::vec& a_exp,
    const arma::vec& b_exp,
    const arma::vec& a_ref,
    const arma::vec& b_ref
) {
  return compute_ba_posterior_impl(a_exp, b_exp, a_ref, b_ref);
}

// =============================================================================
// MAIN: Compute predictive probability (full MC with antithetic variates)
// =============================================================================

//' Compute predictive probability for hybrid trial conversion (C++)
//'
//' Full Monte Carlo computation with optional antithetic variates and early stopping.
//'
//' @param a_exp Current posterior shape for experimental arm (per interval)
//' @param b_exp Current posterior rate for experimental arm (per interval)
//' @param a_ref Current posterior shape for reference arm (per interval)
//' @param b_ref Current posterior rate for reference arm (per interval)
//' @param n_add Number of additional patients to simulate per arm
//' @param interval_cutpoints PWE interval boundaries
//' @param accrual_rate Patients per time unit
//' @param followup Follow-up time after enrollment
//' @param eff_ba Efficacy threshold for BA decision
//' @param pp_go PP threshold for GO decision
//' @param pp_nogo PP threshold for NO-GO decision
//' @param n_outer Number of outer MC samples
//' @param use_antithetic Use antithetic variates for variance reduction
//' @param use_early_stop Enable early stopping based on CI
//' @return Predictive probability of success
//' @export
// [[Rcpp::export]]
double compute_pp_predictive_cpp(
    const arma::vec& a_exp,
    const arma::vec& b_exp,
    const arma::vec& a_ref,
    const arma::vec& b_ref,
    int n_add,
    const arma::vec& interval_cutpoints,
    double accrual_rate,
    double followup,
    double eff_ba,
    double pp_go,
    double pp_nogo,
    int n_outer = 500,
    bool use_antithetic = true,
    bool use_early_stop = true
) {
  int n_intervals = a_exp.n_elem;

  // Early stopping parameters
  int min_samples_early_stop = 100;
  int early_stop_check_interval = 50;

  int success_count = 0;
  int checked = 0;

  if (use_antithetic) {
    // ANTITHETIC VARIATES: Process pairs of (U, 1-U)
    // Fix: For odd n_outer, only process antithetic sample on last pair if needed
    int n_pairs = n_outer / 2;
    bool process_extra = (n_outer % 2 == 1);  // Need one extra sample for odd n_outer

    for (int i = 0; i < n_pairs + (process_extra ? 1 : 0); i++) {
      // Generate uniform random numbers for inverse CDF transform
      arma::vec u_exp(n_intervals);
      arma::vec u_ref(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        u_exp(j) = R::runif(0.0, 1.0);
        u_ref(j) = R::runif(0.0, 1.0);
      }

      // Process both original and antithetic sample
      // For the extra iteration (odd n_outer), only process the original sample
      int max_anti = (process_extra && i == n_pairs) ? 1 : 2;
      for (int anti = 0; anti < max_anti; anti++) {
        // Draw "true" hazards from posterior using inverse CDF
        arma::vec lambda_exp_true(n_intervals);
        arma::vec lambda_ref_true(n_intervals);

        for (int j = 0; j < n_intervals; j++) {
          double u_e = (anti == 0) ? u_exp(j) : (1.0 - u_exp(j));
          double u_r = (anti == 0) ? u_ref(j) : (1.0 - u_ref(j));

          // qgamma: inverse CDF of Gamma distribution
          // R::qgamma(p, shape, scale, lower.tail, log.p)
          lambda_exp_true(j) = R::qgamma(u_e, a_exp(j), 1.0/b_exp(j), 1, 0);
          lambda_ref_true(j) = R::qgamma(u_r, a_ref(j), 1.0/b_ref(j), 1, 0);

          // Ensure valid hazards
          if (lambda_exp_true(j) < 1e-6) lambda_exp_true(j) = 1e-6;
          if (lambda_ref_true(j) < 1e-6) lambda_ref_true(j) = 1e-6;
        }

        // Simulate future events/exposure
        List future_exp = simulate_future_arm_pwe_impl(
          n_add, lambda_exp_true, interval_cutpoints, accrual_rate, followup
        );
        List future_ref = simulate_future_arm_pwe_impl(
          n_add, lambda_ref_true, interval_cutpoints, accrual_rate, followup
        );

        arma::vec events_exp = as<arma::vec>(future_exp["events"]);
        arma::vec exposure_exp = as<arma::vec>(future_exp["exposure"]);
        arma::vec events_ref = as<arma::vec>(future_ref["events"]);
        arma::vec exposure_ref = as<arma::vec>(future_ref["exposure"]);

        // Update posteriors with future data
        arma::vec a_exp_final = a_exp + events_exp;
        arma::vec b_exp_final = b_exp + exposure_exp;
        arma::vec a_ref_final = a_ref + events_ref;
        arma::vec b_ref_final = b_ref + exposure_ref;

        // Compute P(HR < 1 | final data)
        double p_between = compute_ba_posterior_impl(
          a_exp_final, b_exp_final, a_ref_final, b_ref_final
        );

        // Check success criterion
        if (!ISNA(p_between) && p_between >= eff_ba) {
          success_count++;
        }

        checked++;
      }

      // Early stopping check (only after complete pairs)
      if (use_early_stop && checked >= min_samples_early_stop &&
          (checked % early_stop_check_interval) < 2) {
        double current_pp = (double)success_count / checked;
        double se = std::sqrt(current_pp * (1.0 - current_pp) / checked);

        // Stop if clearly above pp_go or clearly below pp_nogo
        if (current_pp - 2.0 * se > pp_go || current_pp + 2.0 * se < pp_nogo) {
          break;
        }
      }
    }
  } else {
    // STANDARD MC: Independent samples
    for (int i = 0; i < n_outer; i++) {
      // Draw "true" hazards from posterior
      arma::vec lambda_exp_true(n_intervals);
      arma::vec lambda_ref_true(n_intervals);

      for (int j = 0; j < n_intervals; j++) {
        // rgamma: R::rgamma(shape, scale)
        lambda_exp_true(j) = R::rgamma(a_exp(j), 1.0/b_exp(j));
        lambda_ref_true(j) = R::rgamma(a_ref(j), 1.0/b_ref(j));

        // Ensure valid hazards
        if (lambda_exp_true(j) < 1e-6) lambda_exp_true(j) = 1e-6;
        if (lambda_ref_true(j) < 1e-6) lambda_ref_true(j) = 1e-6;
      }

      // Simulate future events/exposure
      List future_exp = simulate_future_arm_pwe_impl(
        n_add, lambda_exp_true, interval_cutpoints, accrual_rate, followup
      );
      List future_ref = simulate_future_arm_pwe_impl(
        n_add, lambda_ref_true, interval_cutpoints, accrual_rate, followup
      );

      arma::vec events_exp = as<arma::vec>(future_exp["events"]);
      arma::vec exposure_exp = as<arma::vec>(future_exp["exposure"]);
      arma::vec events_ref = as<arma::vec>(future_ref["events"]);
      arma::vec exposure_ref = as<arma::vec>(future_ref["exposure"]);

      // Update posteriors with future data
      arma::vec a_exp_final = a_exp + events_exp;
      arma::vec b_exp_final = b_exp + exposure_exp;
      arma::vec a_ref_final = a_ref + events_ref;
      arma::vec b_ref_final = b_ref + exposure_ref;

      // Compute P(HR < 1 | final data)
      double p_between = compute_ba_posterior_impl(
        a_exp_final, b_exp_final, a_ref_final, b_ref_final
      );

      // Check success criterion
      if (!ISNA(p_between) && p_between >= eff_ba) {
        success_count++;
      }

      checked = i + 1;

      // Early stopping check
      if (use_early_stop && checked >= min_samples_early_stop &&
          (checked % early_stop_check_interval) == 0) {
        double current_pp = (double)success_count / checked;
        double se = std::sqrt(current_pp * (1.0 - current_pp) / checked);

        // Stop if clearly above pp_go or clearly below pp_nogo
        if (current_pp - 2.0 * se > pp_go || current_pp + 2.0 * se < pp_nogo) {
          break;
        }
      }
    }
  }

  return (double)success_count / checked;
}

// =============================================================================
// SINGLE-ARM PP: Predictive probability for SA efficacy (vs historical)
// =============================================================================

//' Compute predictive probability for single-arm efficacy (C++)
//'
//' Monte Carlo computation of PP that single arm will reach efficacy threshold
//' vs historical control at the final analysis.
//'
//' @param a_arm Current posterior shape for arm (per interval)
//' @param b_arm Current posterior rate for arm (per interval)
//' @param hist_hazard Historical hazard rates (per interval)
//' @param hr_threshold Hazard ratio threshold for efficacy
//' @param n_add Number of additional patients to simulate
//' @param interval_cutpoints PWE interval boundaries
//' @param accrual_rate Patients per time unit
//' @param followup Follow-up time after enrollment
//' @param eff_threshold Efficacy posterior threshold (e.g., 0.90)
//' @param n_outer Number of outer MC samples
//' @param use_antithetic Use antithetic variates for variance reduction
//' @return Predictive probability of reaching efficacy
//' @export
// [[Rcpp::export]]
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
    int n_outer = 200,
    bool use_antithetic = true
) {
  int n_intervals = a_arm.n_elem;
  double hist_rate = arma::mean(hist_hazard);
  double threshold_rate = hr_threshold * hist_rate;

  int success_count = 0;
  int checked = 0;

  if (use_antithetic) {
    int n_pairs = n_outer / 2;
    bool process_extra = (n_outer % 2 == 1);

    for (int i = 0; i < n_pairs + (process_extra ? 1 : 0); i++) {
      // Generate uniform random numbers
      arma::vec u_arm(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        u_arm(j) = R::runif(0.0, 1.0);
      }

      int max_anti = (process_extra && i == n_pairs) ? 1 : 2;
      for (int anti = 0; anti < max_anti; anti++) {
        // Draw "true" hazards from posterior
        arma::vec lambda_true(n_intervals);
        for (int j = 0; j < n_intervals; j++) {
          double u = (anti == 0) ? u_arm(j) : (1.0 - u_arm(j));
          lambda_true(j) = R::qgamma(u, a_arm(j), 1.0/b_arm(j), 1, 0);
          if (lambda_true(j) < 1e-6) lambda_true(j) = 1e-6;
        }

        // Simulate future events/exposure
        List future = simulate_future_arm_pwe_impl(
          n_add, lambda_true, interval_cutpoints, accrual_rate, followup
        );

        arma::vec events = as<arma::vec>(future["events"]);
        arma::vec exposure = as<arma::vec>(future["exposure"]);

        // Update posteriors with future data
        arma::vec a_final = a_arm + events;
        arma::vec b_final = b_arm + exposure;

        // Compute P(lambda < threshold | final data)
        double a_sum = arma::sum(a_final);
        double b_sum = arma::sum(b_final);

        if (a_sum > 0 && b_sum > 0) {
          double p_eff = R::pgamma(threshold_rate, a_sum, 1.0 / b_sum, 1, 0);
          if (p_eff >= eff_threshold) {
            success_count++;
          }
        }
        checked++;
      }
    }
  } else {
    for (int i = 0; i < n_outer; i++) {
      // Draw "true" hazards from posterior
      arma::vec lambda_true(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        lambda_true(j) = R::rgamma(a_arm(j), 1.0/b_arm(j));
        if (lambda_true(j) < 1e-6) lambda_true(j) = 1e-6;
      }

      // Simulate future events/exposure
      List future = simulate_future_arm_pwe_impl(
        n_add, lambda_true, interval_cutpoints, accrual_rate, followup
      );

      arma::vec events = as<arma::vec>(future["events"]);
      arma::vec exposure = as<arma::vec>(future["exposure"]);

      // Update posteriors with future data
      arma::vec a_final = a_arm + events;
      arma::vec b_final = b_arm + exposure;

      // Compute P(lambda < threshold | final data)
      double a_sum = arma::sum(a_final);
      double b_sum = arma::sum(b_final);

      if (a_sum > 0 && b_sum > 0) {
        double p_eff = R::pgamma(threshold_rate, a_sum, 1.0 / b_sum, 1, 0);
        if (p_eff >= eff_threshold) {
          success_count++;
        }
      }
      checked++;
    }
  }

  return (checked > 0) ? (double)success_count / checked : 0.0;
}

// =============================================================================
// SINGLE-ARM PP: Predictive probability for SA futility (vs historical)
// =============================================================================

//' Compute predictive probability for single-arm futility (C++)
//'
//' Monte Carlo computation of PP that single arm will trigger futility
//' vs historical control at the final analysis.
//'
//' @param a_arm Current posterior shape for arm (per interval)
//' @param b_arm Current posterior rate for arm (per interval)
//' @param hist_hazard Historical hazard rates (per interval)
//' @param hr_threshold Hazard ratio threshold for efficacy
//' @param n_add Number of additional patients to simulate
//' @param interval_cutpoints PWE interval boundaries
//' @param accrual_rate Patients per time unit
//' @param followup Follow-up time after enrollment
//' @param fut_threshold Futility posterior threshold (e.g., 0.10)
//' @param n_outer Number of outer MC samples
//' @param use_antithetic Use antithetic variates for variance reduction
//' @return Predictive probability of triggering futility
//' @export
// [[Rcpp::export]]
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
    int n_outer = 200,
    bool use_antithetic = true
) {
  int n_intervals = a_arm.n_elem;
  double hist_rate = arma::mean(hist_hazard);
  double threshold_rate = hr_threshold * hist_rate;

  int futility_count = 0;
  int checked = 0;

  if (use_antithetic) {
    int n_pairs = n_outer / 2;
    bool process_extra = (n_outer % 2 == 1);

    for (int i = 0; i < n_pairs + (process_extra ? 1 : 0); i++) {
      // Generate uniform random numbers
      arma::vec u_arm(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        u_arm(j) = R::runif(0.0, 1.0);
      }

      int max_anti = (process_extra && i == n_pairs) ? 1 : 2;
      for (int anti = 0; anti < max_anti; anti++) {
        // Draw "true" hazards from posterior
        arma::vec lambda_true(n_intervals);
        for (int j = 0; j < n_intervals; j++) {
          double u = (anti == 0) ? u_arm(j) : (1.0 - u_arm(j));
          lambda_true(j) = R::qgamma(u, a_arm(j), 1.0/b_arm(j), 1, 0);
          if (lambda_true(j) < 1e-6) lambda_true(j) = 1e-6;
        }

        // Simulate future events/exposure
        List future = simulate_future_arm_pwe_impl(
          n_add, lambda_true, interval_cutpoints, accrual_rate, followup
        );

        arma::vec events = as<arma::vec>(future["events"]);
        arma::vec exposure = as<arma::vec>(future["exposure"]);

        // Update posteriors with future data
        arma::vec a_final = a_arm + events;
        arma::vec b_final = b_arm + exposure;

        // Compute P(lambda < threshold | final data)
        double a_sum = arma::sum(a_final);
        double b_sum = arma::sum(b_final);

        if (a_sum > 0 && b_sum > 0) {
          double p_eff = R::pgamma(threshold_rate, a_sum, 1.0 / b_sum, 1, 0);
          // Futility triggered when p_eff <= fut_threshold
          if (p_eff <= fut_threshold) {
            futility_count++;
          }
        }
        checked++;
      }
    }
  } else {
    for (int i = 0; i < n_outer; i++) {
      // Draw "true" hazards from posterior
      arma::vec lambda_true(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        lambda_true(j) = R::rgamma(a_arm(j), 1.0/b_arm(j));
        if (lambda_true(j) < 1e-6) lambda_true(j) = 1e-6;
      }

      // Simulate future events/exposure
      List future = simulate_future_arm_pwe_impl(
        n_add, lambda_true, interval_cutpoints, accrual_rate, followup
      );

      arma::vec events = as<arma::vec>(future["events"]);
      arma::vec exposure = as<arma::vec>(future["exposure"]);

      // Update posteriors with future data
      arma::vec a_final = a_arm + events;
      arma::vec b_final = b_arm + exposure;

      // Compute P(lambda < threshold | final data)
      double a_sum = arma::sum(a_final);
      double b_sum = arma::sum(b_final);

      if (a_sum > 0 && b_sum > 0) {
        double p_eff = R::pgamma(threshold_rate, a_sum, 1.0 / b_sum, 1, 0);
        // Futility triggered when p_eff <= fut_threshold
        if (p_eff <= fut_threshold) {
          futility_count++;
        }
      }
      checked++;
    }
  }

  return (checked > 0) ? (double)futility_count / checked : 0.0;
}
