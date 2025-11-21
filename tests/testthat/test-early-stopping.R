
library(testthat)
library(evolveTrial)

test_that("Early stopping is triggered when lenient thresholds are set", {
  # Scenario with high probability of early stopping for futility
  base_args <- list(
    num_simulations = 100,
    arm_names = c("Arm A"),
    reference_arm_name = "Arm A",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c("Arm A" = 1),
    weibull_median_true_arms = c("Arm A" = 1), # Very short survival
    null_median_arms = c("Arm A" = 12),
    futility_median_arms = c("Arm A" = 10),    # Futility threshold
    interval_cutpoints_sim = seq(0, 24, by = 3),
    max_follow_up_sim = 24,
    censor_max_time_sim = 24,
    prior_alpha_params_model = rep(0.5, 8),
    prior_beta_params_model = rep(0.5, 8),
    num_posterior_draws = 500,
    cohort_size_per_arm = 1,
    max_total_patients_per_arm = c("Arm A" = 20),
    futility_stopping_rule_hc = TRUE,
    futility_threshold_hc_prob = 0.99, # High probability for futility
    overall_accrual_rate = 10,
    randomization_probs = c("Arm A" = 1),
    min_follow_up_at_final = 0,
    min_events_hc = 0,
    min_median_followup_hc = 0,
    interim_calendar_beat = 1
  )

  # With min_patients_for_analysis > max_total_patients_per_arm, PET should be 0
  args_bug <- base_args
  args_bug$min_patients_for_analysis <- 30
  
  results_bug <- do.call(run_simulation_pure, args_bug)
  
  expect_equal(results_bug$PET_Futility, 0)

  # With min_patients_for_analysis = 0, PET should be > 0
  args_fix <- base_args
  args_fix$min_patients_for_analysis <- 0
  
  results_fix <- do.call(run_simulation_pure, args_fix)

  expect_true(results_fix$PET_Futility > 0)
})
