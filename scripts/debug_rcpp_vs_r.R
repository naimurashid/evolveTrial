#!/usr/bin/env Rscript
#' Debug script to compare Rcpp vs R hybrid trial simulation
#' Identifies specific differences in behavior

library(evolveTrial)

# Set seed for reproducibility
set.seed(42)

# Create matched parameters
hybrid_theta <- list(
  eff_sa = 0.90,
  fut_sa = 0.10,
  hr_threshold_sa = 0.80,
  ev_sa = 15L,
  nmax_sa = 50L,
  conversion_trigger = "any_single_success",
  pp_go = 0.70,
  pp_nogo = 0.20,
  ss_method = "posterior",
  max_additional_n = 60L,
  n_add_candidates = seq(10, 100, by = 10),
  eff_ba = 0.975,
  fut_ba = 0.05,
  ev_ba = 15L,
  nmax_ba = 80L,
  futility_action = "drop_arm",
  prior_strength = 0.5,
  n_outer = 200L,
  n_inner = 50L
)

base_args <- list(
  n_intervals = 4L,
  interval_cutpoints_sim = c(0, 6, 12, 18, 24),
  overall_accrual_rate = 2.0,
  max_follow_up_sim = 24,
  max_trial_time = 72,
  prior_alpha_params_model = rep(0.5, 4),
  prior_beta_params_model = rep(0.5, 4)
)

# Scenario with treatment effect
scenario_params <- list(
  historical_median = 6,
  ref_median = 7.5,
  exp_median = 9  # Treatment benefit
)

# Convert to lambda for C++
n_intervals <- 4
lambda_hist <- rep(log(2) / scenario_params$historical_median, n_intervals)
lambda_ref <- rep(log(2) / scenario_params$ref_median, n_intervals)
lambda_exp <- rep(log(2) / scenario_params$exp_median, n_intervals)
lambda_null <- lambda_ref  # Null: same as reference

cat("=== COMPARISON: Rcpp vs R Hybrid Trial Simulation ===\n\n")

cat("1. PARAMETER SUMMARY:\n")
cat("   Efficacy SA:", hybrid_theta$eff_sa, "\n")
cat("   Futility SA:", hybrid_theta$fut_sa, "\n")
cat("   PP Go:", hybrid_theta$pp_go, "\n")
cat("   PP No-Go:", hybrid_theta$pp_nogo, "\n")
cat("   nmax_sa:", hybrid_theta$nmax_sa, "\n")
cat("   nmax_ba:", hybrid_theta$nmax_ba, "\n")
cat("\n")

# Run R simulations
cat("2. RUNNING R SIMULATIONS (200 reps)...\n")

# Need to source the R implementation
# For now, let's use the pure R simulation if available
run_pure_r_simulation <- function(hybrid_theta, base_args, scenario_params, n_sim = 200) {

  arm_names <- c("Reference", "Experimental")
  reference_arm <- "Reference"

  results <- data.frame(
    outcome = character(n_sim),
    total_n = integer(n_sim),
    is_success = logical(n_sim),
    converted = logical(n_sim),
    conversion_decision = character(n_sim),
    sa_efficacy = logical(n_sim),
    pp_at_conversion = numeric(n_sim),
    stringsAsFactors = FALSE
  )

  for (i in 1:n_sim) {
    # Create initial state
    state <- create_hybrid_state(arm_names, reference_arm, hybrid_theta, base_args)

    # Simulate trial
    look_interval <- 3
    current_time <- 0
    accrual_rate <- base_args$overall_accrual_rate
    max_trial_time <- base_args$max_trial_time

    # Lambda vectors for simulation
    lambda_exp <- rep(log(2) / scenario_params$exp_median, base_args$n_intervals)
    lambda_ref <- rep(log(2) / scenario_params$ref_median, base_args$n_intervals)

    while (state$current_state != "STATE_STOP" && current_time < max_trial_time) {
      current_time <- current_time + look_interval
      state$current_time <- current_time
      state$interim_count <- state$interim_count + 1

      # Enroll patients
      n_per_arm <- round(accrual_rate * look_interval / 2)

      for (arm in state$active_arms) {
        max_n <- if (state$current_state == "STATE_SINGLE") {
          hybrid_theta$nmax_sa
        } else {
          hybrid_theta$nmax_ba
        }

        can_enroll <- min(n_per_arm, max_n - state$n_enrolled[arm])
        if (can_enroll <= 0) next

        # Simulate patients
        lambda <- if (arm == "Reference") lambda_ref else lambda_exp
        surv_times <- simulate_pwe_survival_batch(can_enroll, lambda, base_args$interval_cutpoints_sim)

        for (j in 1:can_enroll) {
          enroll_time <- current_time - runif(1, 0, look_interval)
          state$registries[[arm]] <- rbind(state$registries[[arm]], data.frame(
            patient_id = nrow(state$registries[[arm]]) + 1,
            enrollment_time = enroll_time,
            event_time = enroll_time + surv_times[j],
            observed_time = min(surv_times[j], max_trial_time - enroll_time),
            event = surv_times[j] <= (max_trial_time - enroll_time)
          ))
          state$n_enrolled[arm] <- state$n_enrolled[arm] + 1
          if (state$current_state == "STATE_SINGLE") {
            state$n_enrolled_phase1[arm] <- state$n_enrolled_phase1[arm] + 1
          }
        }
      }

      # Update state
      state <- update_hybrid_state(state, hybrid_theta, base_args, scenario_params)
    }

    # Record results
    results$outcome[i] <- state$trial_outcome
    results$total_n[i] <- sum(state$n_enrolled)
    results$is_success[i] <- is_trial_success(state, hybrid_theta)
    results$converted[i] <- state$conversion_decision == "GO"
    results$conversion_decision[i] <- state$conversion_decision %||% "NA"
    results$sa_efficacy[i] <- any(state$sa_efficacy_reached)
    results$pp_at_conversion[i] <- state$pp_at_conversion %||% NA_real_
  }

  results
}

# Run the R simulation
set.seed(42)
r_results <- tryCatch({
  run_pure_r_simulation(hybrid_theta, base_args, scenario_params, n_sim = 200)
}, error = function(e) {
  cat("   R simulation error:", e$message, "\n")
  NULL
})

if (!is.null(r_results)) {
  cat("   R Power:", mean(r_results$is_success), "\n")
  cat("   R EN_alt:", mean(r_results$total_n), "\n")
  cat("   R P(conversion):", mean(r_results$converted), "\n")
  cat("   R SA efficacy rate:", mean(r_results$sa_efficacy), "\n")
  cat("\n")
  cat("   Outcome distribution:\n")
  print(table(r_results$outcome))
  cat("\n")
  cat("   Conversion decision distribution:\n")
  print(table(r_results$conversion_decision))
}

# Run Rcpp simulations
cat("\n3. RUNNING RCPP SIMULATIONS (200 reps)...\n")

# Build C++ compatible parameters
theta_cpp <- list(
  eff_sa = hybrid_theta$eff_sa,
  fut_sa = hybrid_theta$fut_sa,
  hr_threshold_sa = hybrid_theta$hr_threshold_sa,
  ev_sa = as.integer(hybrid_theta$ev_sa),
  nmax_sa = as.integer(hybrid_theta$nmax_sa),
  conversion_trigger = hybrid_theta$conversion_trigger,
  pp_go = hybrid_theta$pp_go,
  pp_nogo = hybrid_theta$pp_nogo,
  ss_method = hybrid_theta$ss_method,
  max_additional_n = as.integer(hybrid_theta$max_additional_n),
  eff_ba = hybrid_theta$eff_ba,
  fut_ba = hybrid_theta$fut_ba,
  ev_ba = as.integer(hybrid_theta$ev_ba),
  nmax_ba = as.integer(hybrid_theta$nmax_ba),
  futility_action = hybrid_theta$futility_action,
  prior_strength = hybrid_theta$prior_strength,
  n_outer = as.integer(hybrid_theta$n_outer)
)

scenario_cpp <- list(
  lambda_hist = lambda_hist,
  lambda_ref = lambda_ref,
  lambda_exp = lambda_exp
)

set.seed(42)
cpp_results <- tryCatch({
  run_hybrid_simulations_cpp(
    n_sim = 200L,
    theta_list = theta_cpp,
    base_args_list = base_args,
    scenario_params_list = scenario_cpp
  )
}, error = function(e) {
  cat("   C++ simulation error:", e$message, "\n")
  NULL
})

if (!is.null(cpp_results)) {
  cat("   C++ Power:", mean(cpp_results$is_success), "\n")
  cat("   C++ EN_alt:", mean(cpp_results$total_n), "\n")
  cat("   C++ P(conversion):", mean(cpp_results$converted), "\n")
  cat("   C++ SA efficacy rate:", mean(cpp_results$sa_efficacy), "\n")
  cat("\n")
  cat("   Outcome distribution:\n")
  print(table(cpp_results$outcome))
}

# Summary comparison
cat("\n4. COMPARISON SUMMARY:\n")
if (!is.null(r_results) && !is.null(cpp_results)) {
  cat("   Power diff (Rcpp - R):", mean(cpp_results$is_success) - mean(r_results$is_success), "\n")
  cat("   EN_alt diff (Rcpp - R):", mean(cpp_results$total_n) - mean(r_results$total_n), "\n")
  cat("   Conversion diff (Rcpp - R):", mean(cpp_results$converted) - mean(r_results$converted), "\n")
}

cat("\n5. CHECKING KEY LOGIC DIFFERENCES:\n")
cat("   Issue 1: R computes PP for multiple n_add candidates, C++ uses single value\n")
cat("   Issue 2: R immediately stops on ambiguous PP, C++ may continue\n")
cat("   Issue 3: R uses 'conversion_evaluated' flag, C++ re-evaluates each look\n")
