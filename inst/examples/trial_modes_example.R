#!/usr/bin/env Rscript
# Example: Using Different Trial Modes in evolveTrial
#
# This example demonstrates the four trial modes available in the C++ simulator:
# 1. single_arm - Single-arm vs historical control
# 2. between_arm - Randomized comparison (BA only)
# 3. hybrid - SA→conversion→BA seamless design
# 4. dual_single_arm - Two independent single-arm evaluations

library(evolveTrial)

# ============================================================================
# SETUP: Common configuration for all examples
# ============================================================================

# Base trial parameters
base_args <- list(
  interval_cutpoints_sim = c(0, 6, 12, 24),  # 3 intervals: 0-6, 6-12, 12-24 months
  overall_accrual_rate = 2.0,                 # 2 patients/month
  max_follow_up_sim = 12,                     # 12 months follow-up
  look_interval = 3.0,                        # Interim looks every 3 months
  max_trial_time = 72,                        # Max 72 months total
  prior_alpha_params_model = rep(0.5, 3),     # Weakly informative Gamma prior
  prior_beta_params_model = rep(0.5, 3)
)

# Design thresholds
theta <- list(
  # Single-arm thresholds
  eff_sa = 0.90,           # SA efficacy: P(HR < threshold | data) > 0.90

  fut_sa = 0.10,           # SA futility: P(HR < threshold | data) < 0.10
  hr_threshold_sa = 0.70,  # Target hazard ratio vs historical
  ev_sa = 10,              # Minimum 10 events for SA decision
  nmax_sa = 40,            # Max 40 patients in SA phase

  # Between-arm thresholds
  eff_ba = 0.90,           # BA efficacy: P(exp better | data) > 0.90
  fut_ba = 0.10,           # BA futility: P(exp better | data) < 0.10
  ev_ba = 20,              # Minimum 20 events for BA decision
  nmax_ba = 80,            # Max 80 patients in BA phase

  # Conversion parameters (for hybrid mode)
  pp_go = 0.50,            # PP > 0.50 → convert to BA
  pp_nogo = 0.30,          # PP < 0.30 → stop trial

  # Futility action
  futility_action = "drop_arm"  # Options: "drop_arm", "stop_trial", "continue"
)

# Hazard rates (piecewise exponential)
lambda_exp <- c(0.04, 0.05, 0.06)   # Experimental: better survival
lambda_ref <- c(0.08, 0.09, 0.10)   # Reference: standard of care
lambda_hist <- c(0.05, 0.06, 0.07)  # Historical control

cat("=======================================================\n")
cat("       evolveTrial: Trial Mode Examples\n")
cat("=======================================================\n\n")

# ============================================================================
# EXAMPLE 1: Single-Arm Trial (vs Historical Control)
# ============================================================================

cat("--- Example 1: Single-Arm Trial ---\n")
cat("Design: Single experimental arm compared to historical control\n")
cat("Use case: Phase II oncology trials with historical response rate\n\n")

set.seed(123)
result_sa <- compute_oc_lambda(
  theta = theta,
  base_args = base_args,
  lambda_exp = lambda_exp,
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "single_arm"
)

cat("Results:\n")
cat("  Power:", round(result_sa$power, 3), "\n")
cat("  Type I Error:", round(result_sa$type1, 3), "\n")
cat("  Expected N (alt):", round(result_sa$EN_alt, 1), "\n")
cat("  Expected N (null):", round(result_sa$EN_null, 1), "\n\n")

# ============================================================================
# EXAMPLE 2: Between-Arm (Randomized) Trial
# ============================================================================

cat("--- Example 2: Between-Arm (Randomized) Trial ---\n")
cat("Design: Direct randomized comparison (experimental vs control)\n")
cat("Use case: Phase III confirmatory trials\n\n")

set.seed(456)
result_ba <- compute_oc_lambda(
  theta = theta,
  base_args = base_args,
  lambda_exp = lambda_exp,
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "between_arm"
)

cat("Results:\n")
cat("  Power:", round(result_ba$power, 3), "\n")
cat("  Type I Error:", round(result_ba$type1, 3), "\n")
cat("  Expected N (alt):", round(result_ba$EN_alt, 1), "\n")
cat("  Expected N (null):", round(result_ba$EN_null, 1), "\n\n")

# ============================================================================
# EXAMPLE 3: Hybrid (Seamless SA→BA) Trial
# ============================================================================

cat("--- Example 3: Hybrid (Seamless SA→BA) Trial ---\n")
cat("Design: Start as SA vs historical, convert to BA if promising\n")
cat("Use case: Adaptive Phase II/III designs\n\n")

set.seed(789)
result_hybrid <- compute_oc_lambda(
  theta = theta,
  base_args = base_args,
  lambda_exp = lambda_exp,
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "hybrid"
)

cat("Results:\n")
cat("  Power:", round(result_hybrid$power, 3), "\n")
cat("  Type I Error:", round(result_hybrid$type1, 3), "\n")
cat("  Expected N (alt):", round(result_hybrid$EN_alt, 1), "\n")
cat("  Conversion Rate (alt):", round(result_hybrid$conversion_rate, 3), "\n\n")

# ============================================================================
# EXAMPLE 4: Dual Single-Arm Trial
# ============================================================================

cat("--- Example 4: Dual Single-Arm Trial ---\n")
cat("Design: Both arms evaluated independently vs historical\n")
cat("Use case: Two-arm phase II with external control\n\n")

set.seed(111)
result_dual <- compute_oc_lambda(
  theta = theta,
  base_args = base_args,
  lambda_exp = lambda_exp,
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "dual_single_arm"
)

cat("Results (per-arm metrics):\n")
cat("  Power (experimental):", round(result_dual$power_exp, 3), "\n")
cat("  Power (reference):", round(result_dual$power_ref, 3), "\n")
cat("  Type I (experimental):", round(result_dual$type1_exp, 3), "\n")
cat("  Type I (reference):", round(result_dual$type1_ref, 3), "\n\n")

# ============================================================================
# EXAMPLE 5: Predictive Probability for Decisions
# ============================================================================

cat("--- Example 5: Using Predictive Probability ---\n")
cat("Design: PP-based efficacy decisions, posterior-based futility\n\n")

set.seed(222)
result_pp <- compute_oc_lambda(
  theta = theta,
  base_args = base_args,
  lambda_exp = lambda_exp,
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "single_arm",
  efficacy_method = "predictive",  # PP for efficacy
  futility_method = "posterior"    # Posterior for futility
)

cat("Results (PP efficacy):\n")
cat("  Power:", round(result_pp$power, 3), "\n")
cat("  Type I Error:", round(result_pp$type1, 3), "\n")
cat("  Expected N (alt):", round(result_pp$EN_alt, 1), "\n\n")

# ============================================================================
# EXAMPLE 6: Futility Action Comparison
# ============================================================================

cat("--- Example 6: Futility Action Impact ---\n")
cat("Comparing futility_action = 'drop_arm' vs 'continue'\n\n")

# Scenario with high futility probability
theta_futility <- theta
theta_futility$fut_sa <- 0.40  # Higher futility threshold

set.seed(333)
result_drop <- compute_oc_lambda(
  theta = c(theta_futility, list(futility_action = "drop_arm")),
  base_args = base_args,
  lambda_exp = c(0.07, 0.08, 0.09),  # Similar to historical (bad)
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "single_arm"
)

set.seed(333)
result_continue <- compute_oc_lambda(
  theta = c(theta_futility, list(futility_action = "continue")),
  base_args = base_args,
  lambda_exp = c(0.07, 0.08, 0.09),
  lambda_ref = lambda_ref,
  lambda_hist = lambda_hist,
  num_simulations = 500,
  trial_mode = "single_arm"
)

cat("Results:\n")
cat("  With 'drop_arm':  EN =", round(result_drop$EN_alt, 1), "\n")
cat("  With 'continue':  EN =", round(result_continue$EN_alt, 1), "\n")
cat("  Difference:", round(result_continue$EN_alt - result_drop$EN_alt, 1), "patients\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("=======================================================\n")
cat("                     SUMMARY\n")
cat("=======================================================\n")
cat("Trial Mode        | Power  | Type I | EN(alt)\n")
cat("------------------|--------|--------|--------\n")
cat(sprintf("Single-arm        | %0.3f  | %0.3f  | %5.1f\n",
            result_sa$power, result_sa$type1, result_sa$EN_alt))
cat(sprintf("Between-arm       | %0.3f  | %0.3f  | %5.1f\n",
            result_ba$power, result_ba$type1, result_ba$EN_alt))
cat(sprintf("Hybrid            | %0.3f  | %0.3f  | %5.1f\n",
            result_hybrid$power, result_hybrid$type1, result_hybrid$EN_alt))
cat(sprintf("Dual single-arm   | %0.3f* | %0.3f* | %5.1f\n",
            result_dual$power_exp, result_dual$type1_exp, result_dual$EN_alt))
cat("*Per experimental arm\n\n")
cat("Performance: ~200 simulations/second on typical hardware\n")
