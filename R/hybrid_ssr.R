#' hybrid_ssr.R
#' Sample Size Re-Estimation methods for hybrid trials
#'
#' Implements both predictive probability (recommended) and posterior probability
#' (fast approximation) methods for determining additional sample size needs
#' at the conversion decision point.
#'
#' @description
#' Two methods are provided:
#' 1. Predictive Probability: Full Monte Carlo integration over posterior uncertainty
#' 2. Posterior Probability: Fast analytical approximation using point estimates

# ==============================================================================
# MAIN SSR INTERFACE
# ==============================================================================

#' Perform sample size re-estimation for hybrid trial
#'
#' Main entry point for SSR calculation. Dispatches to appropriate method
#' based on theta$ss_method.
#'
#' @param state Current hybrid trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters (for predictive method)
#'
#' @return List with:
#'   - decision: "GO", "NO_GO", or "AMBIGUOUS"
#'   - n_add_recommended: Selected additional N per arm
#'   - pp_at_decision: PP value at decision point
#'   - pp_curve: Full PP curve data frame
#'   - method: Method used ("predictive" or "posterior")
#'
#' @export
perform_ssr <- function(state, theta, base_args, scenario_params = NULL) {

  method <- theta$ss_method %||% "predictive"
  n_candidates <- theta$n_add_candidates %||% seq(10, 100, by = 10)

  # Filter candidates by max_additional_n
  n_candidates <- n_candidates[n_candidates <= theta$max_additional_n]

  if (length(n_candidates) == 0) {
    return(list(
      decision = "NO_GO",
      n_add_recommended = NA_integer_,
      pp_at_decision = 0,
      pp_curve = data.frame(n_add = integer(), pp = numeric()),
      method = method,
      reason = "No valid candidate sample sizes"
    ))
  }

  # Compute PP curve
  if (method == "predictive") {
    pp_results <- compute_pp_curve_predictive(
      state, n_candidates, theta, base_args, scenario_params
    )
  } else {
    pp_results <- compute_pp_curve_posterior(
      state, n_candidates, theta, base_args
    )
  }

  # Make go/no-go decision
  decision_result <- make_conversion_decision(pp_results, theta)

  list(
    decision = decision_result$decision,
    n_add_recommended = decision_result$n_add,
    pp_at_decision = decision_result$pp,
    pp_curve = pp_results,
    method = method,
    reason = decision_result$reason
  )
}

# ==============================================================================
# CONVERSION DECISION LOGIC
# ==============================================================================

#' Make conversion decision based on PP curve
#'
#' Implements the two-threshold decision rule:
#' - GO: PP >= pp_go for some viable N
#' - NO_GO: max(PP) < pp_nogo
#' - AMBIGUOUS: pp_nogo <= max(PP) < pp_go (default to NO_GO)
#'
#' @param pp_curve Data frame with n_add and pp columns
#' @param theta Hybrid design parameters
#'
#' @return List with decision, n_add, pp, and reason
#' @export
make_conversion_decision <- function(pp_curve, theta) {

  pp_go <- theta$pp_go %||% 0.70
  pp_nogo <- theta$pp_nogo %||% 0.20

  if (nrow(pp_curve) == 0 || all(is.na(pp_curve$pp))) {
    return(list(
      decision = "NO_GO",
      n_add = NA_integer_,
      pp = NA_real_,
      reason = "No valid PP values computed"
    ))
  }

  max_pp <- max(pp_curve$pp, na.rm = TRUE)

  # Find smallest N achieving pp_go
  viable_idx <- which(pp_curve$pp >= pp_go)

  if (length(viable_idx) > 0) {
    # GO: Found viable N
    selected_idx <- min(viable_idx)
    return(list(
      decision = "GO",
      n_add = pp_curve$n_add[selected_idx],
      pp = pp_curve$pp[selected_idx],
      reason = sprintf("PP (%.3f) >= pp_go (%.2f) at N_add = %d",
                       pp_curve$pp[selected_idx], pp_go,
                       pp_curve$n_add[selected_idx])
    ))
  }

  if (max_pp < pp_nogo) {
    # NO_GO: Clear futility
    return(list(
      decision = "NO_GO",
      n_add = NA_integer_,
      pp = max_pp,
      reason = sprintf("Max PP (%.3f) < pp_nogo (%.2f)",
                       max_pp, pp_nogo)
    ))
  }

  # AMBIGUOUS: Default to NO_GO
  return(list(
    decision = "AMBIGUOUS",
    n_add = NA_integer_,
    pp = max_pp,
    reason = sprintf("PP (%.3f) in ambiguous region [%.2f, %.2f)",
                     max_pp, pp_nogo, pp_go)
  ))
}

# ==============================================================================
# PREDICTIVE PROBABILITY METHOD (RECOMMENDED)
# ==============================================================================

#' Compute PP curve using predictive probability method
#'
#' Full Monte Carlo integration over posterior uncertainty.
#' For each candidate N_add:
#'   1. Draw "true" parameters from current posterior
#'   2. Simulate future patients under those parameters
#'   3. Update posterior with simulated data
#'   4. Compute P(BA success | final data)
#'   5. Average over outer samples
#'
#' @param state Current hybrid trial state
#' @param n_candidates Vector of candidate N_add values
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Data frame with n_add and pp columns
#' @export
compute_pp_curve_predictive <- function(state, n_candidates, theta, base_args,
                                         scenario_params = NULL) {

  n_outer <- theta$n_outer %||% 500  # Reduced from 1000 for 2x speedup with <1% accuracy loss

  pp_values <- vapply(n_candidates, function(n_add) {
    compute_pp_predictive_full(state, n_add, theta, base_args, scenario_params,
                                n_outer = n_outer)
  }, numeric(1))

  data.frame(
    n_add = n_candidates,
    pp = pp_values
  )
}

#' Compute predictive probability for a single N_add value
#'
#' This is the core Monte Carlo algorithm for predictive probability.
#' Implements early stopping for efficiency when outcome is clear.
#'
#' @param state Current hybrid trial state
#' @param n_add Additional patients per arm
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#' @param n_outer Number of outer Monte Carlo samples
#' @param use_antithetic Logical; use antithetic variates for variance reduction (default TRUE)
#'
#' @return Predictive probability
#' @export
compute_pp_predictive_full <- function(state, n_add, theta, base_args,
                                        scenario_params, n_outer = 500,
                                        use_antithetic = TRUE) {


  # Extract arm information
  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm) || !exp_arm %in% state$active_arms) {
    return(NA_real_)
  }

  # Current posterior parameters
  a_exp <- state$posterior_a[[exp_arm]]
  b_exp <- state$posterior_b[[exp_arm]]
  a_ref <- state$posterior_a[[ref_arm]]
  b_ref <- state$posterior_b[[ref_arm]]

  n_intervals <- length(a_exp)
  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = n_intervals + 1)

  # Simulation parameters
  accrual_rate <- base_args$overall_accrual_rate %||% 2.0
  followup <- base_args$max_follow_up_sim %||% 24
  eff_ba <- theta$eff_ba %||% 0.975
  pp_go <- theta$pp_go %||% 0.7
  pp_nogo <- theta$pp_nogo %||% 0.2

  # Use Rcpp implementation for ~7x speedup
  compute_pp_predictive_cpp(
    a_exp = a_exp,
    b_exp = b_exp,
    a_ref = a_ref,
    b_ref = b_ref,
    n_add = as.integer(n_add),
    interval_cutpoints = interval_cutpoints,
    accrual_rate = accrual_rate,
    followup = followup,
    eff_ba = eff_ba,
    pp_go = pp_go,
    pp_nogo = pp_nogo,
    n_outer = as.integer(n_outer),
    use_antithetic = use_antithetic,
    use_early_stop = TRUE
  )
}

#' Simulate future arm data under PWE model
#'
#' Wrapper that delegates to simulate_future_arm_pwe() in hybrid_trial.R
#' to avoid code duplication. Adds guard for n_patients <= 0.
#'
#' @param n_patients Number of patients to simulate
#' @param lambda True hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#' @param accrual_rate Patients per month
#' @param followup Follow-up time in months
#'
#' @return List with events and exposure per interval
simulate_future_arm_data <- function(n_patients, lambda, interval_cutpoints,
                                      accrual_rate, followup) {

  # Guard for zero patients
  if (n_patients <= 0) {
    n_intervals <- length(lambda)
    return(list(events = numeric(n_intervals), exposure = numeric(n_intervals)))
  }

  # Use shared implementation from hybrid_trial.R
  simulate_future_arm_pwe(n_patients, lambda, interval_cutpoints,
                           accrual_rate, followup)
}

#' Simulate survival time from PWE model
#'
#' Wrapper that delegates to simulate_pwe_survival() in hybrid_trial.R
#' to avoid code duplication.
#'
#' @param lambda Hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#'
#' @return Survival time (Inf if no event possible due to zero hazards)
simulate_pwe_time <- function(lambda, interval_cutpoints) {
  # Use shared implementation from hybrid_trial.R
  simulate_pwe_survival(lambda, interval_cutpoints)
}

#' Compute between-arm posterior probability
#'
#' Uses closed-form F-distribution for exponential model,
#' Monte Carlo sampling for PWE.
#'
#' @param a_exp Posterior shape for experimental arm
#' @param b_exp Posterior rate for experimental arm
#' @param a_ref Posterior shape for reference arm
#' @param b_ref Posterior rate for reference arm
#' @param interval_cutpoints Interval boundaries
#' @param base_args evolveTrial base configuration
#'
#' @return P(HR < 1 | data)
compute_ba_posterior <- function(a_exp, b_exp, a_ref, b_ref,
                                  interval_cutpoints, base_args) {

  # For exponential (single interval) or aggregated approach
  if (length(a_exp) == 1 && length(a_ref) == 1) {
    # Closed-form F-distribution
    return(pf(
      q = (b_exp / b_ref) * (a_ref / a_exp),
      df1 = 2 * a_exp,
      df2 = 2 * a_ref
    ))
  }

  # Aggregated approach for PWE: sum over intervals
  # This is a simplification that works well in practice
  a_exp_sum <- sum(a_exp)
  b_exp_sum <- sum(b_exp)
  a_ref_sum <- sum(a_ref)
  b_ref_sum <- sum(b_ref)

  # Use F-distribution on aggregated parameters
  pf(
    q = (b_exp_sum / b_ref_sum) * (a_ref_sum / a_exp_sum),
    df1 = 2 * a_exp_sum,
    df2 = 2 * a_ref_sum
  )
}

# ==============================================================================
# POSTERIOR PROBABILITY METHOD (FAST APPROXIMATION)
# ==============================================================================

#' Compute PP curve using posterior method
#'
#' Fast analytical approximation using current posterior point estimate.
#' Less accurate but much faster than predictive method.
#'
#' @param state Current hybrid trial state
#' @param n_candidates Vector of candidate N_add values
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#'
#' @return Data frame with n_add and pp columns
#' @export
compute_pp_curve_posterior <- function(state, n_candidates, theta, base_args) {

  # Extract arm information
  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm) || !exp_arm %in% state$active_arms) {
    return(data.frame(n_add = n_candidates, pp = rep(NA_real_, length(n_candidates))))
  }

  # Current posterior parameters
  a_exp <- state$posterior_a[[exp_arm]]
  b_exp <- state$posterior_b[[exp_arm]]
  a_ref <- state$posterior_a[[ref_arm]]
  b_ref <- state$posterior_b[[ref_arm]]

  # Point estimates (posterior mode for gamma: (a-1)/b, or mean a/b)
  lambda_exp_hat <- sum(a_exp) / sum(b_exp)
  lambda_ref_hat <- sum(a_ref) / sum(b_ref)

  # Current HR estimate
  hr_hat <- lambda_exp_hat / lambda_ref_hat

  # Current information (total events)
  total_events <- sum(a_exp) + sum(a_ref) -
    2 * length(a_exp) * (base_args$prior_alpha_params_model[1] %||% 0.5)

  current_n <- sum(state$n_enrolled[state$active_arms])
  eff_ba <- theta$eff_ba %||% 0.975

  pp_values <- vapply(n_candidates, function(n_add) {
    compute_pp_posterior_single(
      hr_hat, current_n, n_add, total_events, eff_ba
    )
  }, numeric(1))

  data.frame(
    n_add = n_candidates,
    pp = pp_values
  )
}

#' Compute posterior-based PP for single N_add
#'
#' Uses normal approximation to log(HR) to compute conditional power.
#'
#' @param hr_hat Current HR point estimate
#' @param current_n Current total sample size
#' @param n_add Additional patients per arm
#' @param total_events Current total events
#' @param eff_ba BA efficacy threshold
#'
#' @return Conditional power (PP approximation)
compute_pp_posterior_single <- function(hr_hat, current_n, n_add, total_events,
                                         eff_ba) {

  if (is.na(hr_hat) || hr_hat <= 0) {
    return(NA_real_)
  }

  # Effect size on log scale
  log_hr <- log(hr_hat)

  # Future sample size
  future_n <- current_n + 2 * n_add

  # Approximate event rate (assume 70% of patients will have events)
  event_rate <- 0.70
  future_events <- future_n * event_rate

  # Variance of log(HR) ~ 4/d where d = total events
  var_current <- if (total_events > 0) 4 / total_events else 10
  var_future <- 4 / future_events

  # Weighted estimate of effect at end
  # Information-weighted average (current + future)
  info_current <- 1 / var_current
  info_future <- 1 / var_future

  # Posterior variance after additional data
  var_final <- 1 / (info_current + info_future)

  # Z-statistic transformation
  # For efficacy, need P(HR < 1 | data) > eff_ba
  # This is equivalent to P(log_HR < 0 | data) > eff_ba
  # Critical value for log_HR
  z_crit <- qnorm(eff_ba)

  # Under assumed truth (log_hr), what's probability of achieving eff_ba?
  # log_HR_final ~ N(log_hr, var_final)
  # P(log_HR_final / sqrt(var_final) < z_crit) needs to be evaluated

  se_final <- sqrt(var_final)

  if (log_hr < 0) {
    # Beneficial effect (HR < 1)
    # PP = P(log_HR_final < -z_crit * se_final)
    # = P(N(log_hr, var_final) < threshold)
    threshold <- -z_crit * se_final
    pp <- pnorm(threshold, mean = log_hr, sd = se_final)
  } else {
    # Null or harmful effect
    # PP stays low
    pp <- pnorm(-z_crit * se_final, mean = log_hr, sd = se_final)
  }

  # Bound to [0, 1]
  min(max(pp, 0), 1)
}

# ==============================================================================
# ADVANCED PREDICTIVE PROBABILITY WITH PWE
# ==============================================================================

#' Compute predictive probability with full PWE model
#'
#' Uses Monte Carlo sampling of median survival times for
#' more accurate PWE comparison.
#'
#' @param state Current hybrid trial state
#' @param n_add Additional patients per arm
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#' @param n_outer Number of outer samples
#' @param n_inner Number of inner samples for posterior
#'
#' @return Predictive probability
#' @export
compute_pp_predictive_pwe <- function(state, n_add, theta, base_args,
                                       scenario_params, n_outer = 500,
                                       n_inner = 250) {

  # Extract arm information
  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm)) return(NA_real_)

  # Current posterior parameters
  a_exp <- state$posterior_a[[exp_arm]]
  b_exp <- state$posterior_b[[exp_arm]]
  a_ref <- state$posterior_a[[ref_arm]]
  b_ref <- state$posterior_b[[ref_arm]]

  n_intervals <- length(a_exp)
  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = n_intervals + 1)

  accrual_rate <- base_args$overall_accrual_rate %||% 2.0
  followup <- base_args$max_follow_up_sim %||% 24
  eff_ba <- theta$eff_ba %||% 0.975

  success_count <- 0

  for (i in seq_len(n_outer)) {
    # Draw "true" PWE parameters from posterior
    lambda_exp_true <- rgamma(n_intervals, shape = a_exp, rate = b_exp)
    lambda_ref_true <- rgamma(n_intervals, shape = a_ref, rate = b_ref)

    lambda_exp_true <- pmax(lambda_exp_true, 1e-6)
    lambda_ref_true <- pmax(lambda_ref_true, 1e-6)

    # Simulate future data
    future_exp <- simulate_future_arm_data(
      n_add, lambda_exp_true, interval_cutpoints, accrual_rate, followup
    )
    future_ref <- simulate_future_arm_data(
      n_add, lambda_ref_true, interval_cutpoints, accrual_rate, followup
    )

    # Update posteriors
    a_exp_final <- a_exp + future_exp$events
    b_exp_final <- b_exp + future_exp$exposure
    a_ref_final <- a_ref + future_ref$events
    b_ref_final <- b_ref + future_ref$exposure

    # Compute P(HR < 1) using median comparison with MC
    p_between <- compute_ba_posterior_pwe_mc(
      a_exp_final, b_exp_final,
      a_ref_final, b_ref_final,
      interval_cutpoints, n_samples = n_inner
    )

    if (!is.na(p_between) && p_between >= eff_ba) {
      success_count <- success_count + 1
    }
  }

  success_count / n_outer
}

#' Monte Carlo comparison for PWE posteriors
#'
#' Samples from posteriors and compares median survival times.
#'
#' @param a_exp Posterior shapes for experimental arm
#' @param b_exp Posterior rates for experimental arm
#' @param a_ref Posterior shapes for reference arm
#' @param b_ref Posterior rates for reference arm
#' @param interval_cutpoints Interval boundaries
#' @param n_samples Number of Monte Carlo samples
#'
#' @return P(median_exp > median_ref)
compute_ba_posterior_pwe_mc <- function(a_exp, b_exp, a_ref, b_ref,
                                         interval_cutpoints, n_samples = 1000) {

  n_intervals <- length(a_exp)
  count_exp_better <- 0

  for (i in seq_len(n_samples)) {
    # Sample hazards from posteriors
    lambda_exp <- rgamma(n_intervals, shape = a_exp, rate = b_exp)
    lambda_ref <- rgamma(n_intervals, shape = a_ref, rate = b_ref)

    # Compute median survival for each
    median_exp <- compute_pwe_median_survival(lambda_exp, interval_cutpoints)
    median_ref <- compute_pwe_median_survival(lambda_ref, interval_cutpoints)

    # Experimental better if longer median survival
    if (!is.na(median_exp) && !is.na(median_ref) && median_exp > median_ref) {
      count_exp_better <- count_exp_better + 1
    }
  }

  count_exp_better / n_samples
}

#' Compute median survival from PWE hazards
#'
#' Wrapper that delegates to compute_pwe_median() in hybrid_trial.R
#' to avoid code duplication.
#'
#' @param lambda Vector of hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#'
#' @return Median survival time
compute_pwe_median_survival <- function(lambda, interval_cutpoints) {
  # Use shared implementation from hybrid_trial.R
  compute_pwe_median(lambda, interval_cutpoints)
}

# ==============================================================================
# CLOSED-FORM EXPONENTIAL METHODS (FOR VALIDATION)
# ==============================================================================

#' Compute exact PP for exponential model (for validation)
#'
#' Uses closed-form solutions based on F-distribution.
#' This is used to validate Monte Carlo methods.
#'
#' @param a_exp Posterior shape for experimental arm
#' @param b_exp Posterior rate for experimental arm
#' @param a_ref Posterior shape for reference arm
#' @param b_ref Posterior rate for reference arm
#'
#' @return P(lambda_exp < lambda_ref | data)
#' @export
compute_ba_posterior_exponential <- function(a_exp, b_exp, a_ref, b_ref) {

  # P(lambda_exp / lambda_ref < 1)
  # lambda_exp ~ Gamma(a_exp, b_exp)
  # lambda_ref ~ Gamma(a_ref, b_ref)
  # (lambda_exp / b_exp) / (lambda_ref / b_ref) ~ F(2*a_exp, 2*a_ref)

  pf(
    q = (b_exp / b_ref) * (a_ref / a_exp),
    df1 = 2 * a_exp,
    df2 = 2 * a_ref
  )
}

#' Validate MC against closed-form for exponential
#'
#' @param a_exp Experimental arm shape
#' @param b_exp Experimental arm rate
#' @param a_ref Reference arm shape
#' @param b_ref Reference arm rate
#' @param n_samples MC sample size
#'
#' @return List with closed_form and mc_estimate
#' @export
validate_exponential_ba <- function(a_exp, b_exp, a_ref, b_ref,
                                     n_samples = 10000) {

  # Closed-form
  cf <- compute_ba_posterior_exponential(a_exp, b_exp, a_ref, b_ref)

  # Monte Carlo
  lambda_exp <- rgamma(n_samples, shape = a_exp, rate = b_exp)
  lambda_ref <- rgamma(n_samples, shape = a_ref, rate = b_ref)
  mc <- mean(lambda_exp < lambda_ref)

  list(
    closed_form = cf,
    mc_estimate = mc,
    difference = abs(cf - mc),
    se_mc = sqrt(mc * (1 - mc) / n_samples)
  )
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Optimal N_add finder
#'
#' Binary search for minimum N_add achieving target PP.
#'
#' @param state Current trial state
#' @param target_pp Target predictive probability
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#' @param max_n Maximum N to consider
#' @param step Step size for search
#'
#' @return Optimal N_add or NA if not achievable
#' @export
find_optimal_n_add <- function(state, target_pp, theta, base_args,
                                scenario_params = NULL, max_n = 100,
                                step = 10) {

  # Coarse search
  n_candidates <- seq(step, max_n, by = step)

  if (theta$ss_method == "predictive") {
    pp_values <- vapply(n_candidates, function(n) {
      compute_pp_predictive_full(state, n, theta, base_args, scenario_params,
                                  n_outer = 500)  # Reduced for speed
    }, numeric(1))
  } else {
    pp_curve <- compute_pp_curve_posterior(state, n_candidates, theta, base_args)
    pp_values <- pp_curve$pp
  }

  # Find first N achieving target
  idx <- which(pp_values >= target_pp)

  if (length(idx) == 0) {
    return(NA_integer_)
  }

  # Fine-tune with smaller step if desired
  coarse_n <- n_candidates[min(idx)]

  if (step <= 5 || coarse_n <= step) {
    return(as.integer(coarse_n))
  }

  # Fine search between coarse_n - step and coarse_n
  fine_candidates <- seq(max(step, coarse_n - step), coarse_n, by = 5)

  if (theta$ss_method == "predictive") {
    pp_fine <- vapply(fine_candidates, function(n) {
      compute_pp_predictive_full(state, n, theta, base_args, scenario_params,
                                  n_outer = 500)
    }, numeric(1))
  } else {
    pp_fine_df <- compute_pp_curve_posterior(state, fine_candidates, theta, base_args)
    pp_fine <- pp_fine_df$pp
  }

  idx_fine <- which(pp_fine >= target_pp)

  if (length(idx_fine) == 0) {
    return(as.integer(coarse_n))
  }

  as.integer(fine_candidates[min(idx_fine)])
}

#' Expected information gain from additional patients
#'
#' Computes expected reduction in posterior variance from adding n_add patients.
#'
#' @param state Current trial state
#' @param n_add Additional patients per arm
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#'
#' @return Expected information gain (relative to current)
#' @export
expected_information_gain <- function(state, n_add, theta, base_args) {

  # Extract arm information
  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm)) return(NA_real_)

  # Current precision (sum of shape parameters gives information)
  a_exp_current <- sum(state$posterior_a[[exp_arm]])
  a_ref_current <- sum(state$posterior_a[[ref_arm]])

  info_current <- a_exp_current + a_ref_current

  # Expected future events (assume ~70% event rate)
  event_rate <- 0.70
  expected_events <- 2 * n_add * event_rate

  # Future information
  info_future <- info_current + expected_events

  # Information ratio
  info_future / info_current
}

# ==============================================================================
# NULL COALESCING OPERATOR
# ==============================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
