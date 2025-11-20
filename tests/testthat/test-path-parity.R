test_that("single-arm and vs-reference produce same output structure", {
  set.seed(9999)

  base_args <- list(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100
  )

  # Single-arm configuration
  result_hc <- do.call(run_simulation_pure, c(base_args, list(
    compare_arms_option = FALSE,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE
  )))

  # Vs-reference configuration
  result_vsref <- do.call(run_simulation_pure, c(base_args, list(
    compare_arms_option = TRUE,
    min_events_per_arm = 5,
    min_median_followup_per_arm = 2,
    min_person_time_frac_per_arm = 0.15,
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    efficacy_stopping_rule_hc = FALSE,
    futility_stopping_rule_hc = FALSE,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )))

  # Both should return data frames with same structure
  expect_s3_class(result_hc, "data.frame")
  expect_s3_class(result_vsref, "data.frame")
  expect_equal(nrow(result_hc), 2)
  expect_equal(nrow(result_vsref), 2)
  expect_equal(names(result_hc), names(result_vsref))

  # Check that all required columns exist
  expected_cols <- c("Arm_Name", "True_Median", "Type_I_Error_or_Power",
                     "PET_Efficacy", "PET_Futility", "Pr_Reach_Max_N",
                     "Pr_Final_Efficacy", "Pr_Final_Futility",
                     "Pr_Final_Inconclusive", "Exp_N")
  expect_true(all(expected_cols %in% names(result_hc)))
  expect_true(all(expected_cols %in% names(result_vsref)))
})

test_that("both paths respect max_follow_up_sim and max_total_patients_per_arm", {
  set.seed(8888)

  shared_params <- list(
    num_simulations = 25,
    arm_names = c("A", "B"),
    reference_arm_name = "A",
    weibull_shape_true_arms = c(A = 1.5, B = 1.5),
    weibull_median_true_arms = c(A = 10, B = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(A = 0.5, B = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    max_total_patients_per_arm = c(A = 40, B = 40),
    cohort_size_per_arm = c(A = 5, B = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85
  )

  result_hc <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = FALSE,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(A = 8, B = 8),
    futility_median_arms = c(A = 9, B = 9),
    median_pfs_success_threshold_arms = c(A = 11, B = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE
  )))

  result_vsref <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = TRUE,
    min_events_per_arm = 5,
    min_median_followup_per_arm = 2,
    min_person_time_frac_per_arm = 0.15,
    median_pfs_success_threshold_arms = c(A = 11, B = 11),
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    efficacy_stopping_rule_hc = FALSE,
    futility_stopping_rule_hc = FALSE,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )))

  # Both paths should respect max sample size
  expect_true(all(result_hc$Exp_N <= 40))
  expect_true(all(result_vsref$Exp_N <= 40))

  # Expected N should be positive
  expect_true(all(result_hc$Exp_N > 0))
  expect_true(all(result_vsref$Exp_N > 0))
})

test_that("both paths handle parallel execution", {
  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  set.seed(7778)

  shared_params <- list(
    num_simulations = 20,
    arm_names = c("X", "Y"),
    reference_arm_name = "X",
    weibull_shape_true_arms = c(X = 1.5, Y = 1.5),
    weibull_median_true_arms = c(X = 10, Y = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(X = 0.5, Y = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    max_total_patients_per_arm = c(X = 50, Y = 50),
    cohort_size_per_arm = c(X = 5, Y = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100,
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK",
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85
  )

  # Single-arm with parallel
  result_hc_par <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = FALSE,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(X = 8, Y = 8),
    futility_median_arms = c(X = 9, Y = 9),
    median_pfs_success_threshold_arms = c(X = 11, Y = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE
  )))

  # Vs-reference with parallel
  result_vsref_par <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = TRUE,
    min_events_per_arm = 5,
    min_median_followup_per_arm = 2,
    min_person_time_frac_per_arm = 0.15,
    median_pfs_success_threshold_arms = c(X = 11, Y = 11),
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    efficacy_stopping_rule_hc = FALSE,
    futility_stopping_rule_hc = FALSE,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )))

  # Both should complete without errors
  expect_s3_class(result_hc_par, "data.frame")
  expect_s3_class(result_vsref_par, "data.frame")
  expect_equal(nrow(result_hc_par), 2)
  expect_equal(nrow(result_vsref_par), 2)

  # All probabilities should be valid
  expect_true(all(result_hc_par$Type_I_Error_or_Power >= 0 & result_hc_par$Type_I_Error_or_Power <= 1))
  expect_true(all(result_vsref_par$Type_I_Error_or_Power >= 0 & result_vsref_par$Type_I_Error_or_Power <= 1))
})

test_that("both paths work with rebalance_after_events", {
  set.seed(6665)

  shared_params <- list(
    num_simulations = 20,
    arm_names = c("C", "D"),
    reference_arm_name = "C",
    weibull_shape_true_arms = c(C = 1.5, D = 1.5),
    weibull_median_true_arms = c(C = 12, D = 14),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(C = 0.5, D = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    max_total_patients_per_arm = c(C = 40, D = 40),
    cohort_size_per_arm = c(C = 5, D = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    rebalance_after_events = 15,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85
  )

  result_hc_rebal <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = FALSE,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(C = 8, D = 8),
    futility_median_arms = c(C = 9, D = 9),
    median_pfs_success_threshold_arms = c(C = 11, D = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE
  )))

  result_vsref_rebal <- do.call(run_simulation_pure, c(shared_params, list(
    compare_arms_option = TRUE,
    min_events_per_arm = 5,
    min_median_followup_per_arm = 2,
    min_person_time_frac_per_arm = 0.15,
    median_pfs_success_threshold_arms = c(C = 11, D = 11),
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    efficacy_stopping_rule_hc = FALSE,
    futility_stopping_rule_hc = FALSE,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )))

  # Both should complete successfully with rebalancing
  expect_s3_class(result_hc_rebal, "data.frame")
  expect_s3_class(result_vsref_rebal, "data.frame")
  expect_true(all(!is.na(result_hc_rebal$Exp_N)))
  expect_true(all(!is.na(result_vsref_rebal$Exp_N)))
})
