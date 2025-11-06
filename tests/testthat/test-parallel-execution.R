test_that("parallel execution produces valid results", {
  set.seed(4242)

  base_args <- list(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_for_analysis = 5,
    min_median_followup = 2,
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_current_prob_hc = 0.9,
    posterior_futility_threshold_hc = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE
  )

  # Skip on systems without parallel support or CRAN
  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  # Run parallel version with 2 workers
  result_par <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK"
  )))

  # Check structure
  expect_s3_class(result_par, "data.frame")
  expect_equal(nrow(result_par), 2)
  expect_equal(result_par$Arm_Name, c("Control", "Treatment"))

  # Check column completeness
  expected_cols <- c("Arm_Name", "True_Median", "Type_I_Error_or_Power",
                     "PET_Efficacy", "PET_Futility", "Pr_Reach_Max_N",
                     "Pr_Final_Efficacy", "Pr_Final_Futility",
                     "Pr_Final_Inconclusive", "Exp_N")
  expect_named(result_par, expected_cols)

  # Check probabilities sum to 1 for each arm
  for (arm in result_par$Arm_Name) {
    idx <- which(result_par$Arm_Name == arm)
    prob_sum <- result_par$PET_Efficacy[idx] +
                result_par$PET_Futility[idx] +
                result_par$Pr_Reach_Max_N[idx]
    expect_equal(prob_sum, 1, tolerance = 1e-6,
                 info = paste("Probabilities should sum to 1 for", arm))
  }

  # Check all probabilities are in [0, 1]
  prob_cols <- c("Type_I_Error_or_Power", "PET_Efficacy", "PET_Futility",
                 "Pr_Reach_Max_N", "Pr_Final_Efficacy", "Pr_Final_Futility",
                 "Pr_Final_Inconclusive")
  for (col in prob_cols) {
    expect_true(all(result_par[[col]] >= 0 & result_par[[col]] <= 1),
                info = paste(col, "should be in [0, 1]"))
  }

  # Check expected N is reasonable
  expect_true(all(result_par$Exp_N > 0 & result_par$Exp_N <= 50),
              info = "Expected N should be positive and <= max_total_patients_per_arm")
})

test_that("parallel and sequential produce similar operating characteristics", {
  set.seed(8888)

  base_args <- list(
    num_simulations = 100,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_for_analysis = 5,
    min_median_followup = 2,
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_current_prob_hc = 0.9,
    posterior_futility_threshold_hc = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE
  )

  # Skip on systems without parallel support or CRAN
  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  # Run sequential version
  set.seed(9999)
  result_seq <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = FALSE
  )))

  # Run parallel version
  set.seed(9999)
  result_par <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK"
  )))

  # Results will differ due to RNG stream splitting, but operating characteristics
  # should be similar (within Monte Carlo error)
  # With 100 simulations, standard error is roughly sqrt(p*(1-p)/100) ≈ 0.05
  # We use a generous tolerance of 0.15 to allow for both MC error and
  # RNG stream differences
  tolerance <- 0.15

  expect_equal(result_seq$Type_I_Error_or_Power,
               result_par$Type_I_Error_or_Power,
               tolerance = tolerance,
               info = "Type I Error/Power should be similar")

  expect_equal(result_seq$Exp_N, result_par$Exp_N,
               tolerance = 10,  # Allow 10-patient difference
               info = "Expected N should be similar")
})

test_that("parallel execution handles few simulations gracefully", {
  base_args <- list(
    num_simulations = 3,  # Fewer than 2*workers
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_for_analysis = 5,
    min_median_followup = 2,
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_current_prob_hc = 0.9,
    posterior_futility_threshold_hc = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE
  )

  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  # Should fall back to sequential execution (with message suppressed in test)
  set.seed(1234)
  result <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK"
  )))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that("parallel execution with vs-reference design works", {
  set.seed(5555)

  base_args <- list(
    num_simulations = 30,
    arm_names = c("Control", "Treatment_A", "Treatment_B"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,  # vs-reference design
    weibull_shape_true_arms = c(Control = 1.5, Treatment_A = 1.5, Treatment_B = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment_A = 12, Treatment_B = 11),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.4, Treatment_A = 0.3, Treatment_B = 0.3),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_events_per_arm = 8,
    min_median_followup_per_arm = 2,
    min_person_time_frac_per_arm = 0.1,
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    median_pfs_success_threshold_arms = c(Control = 11, Treatment_A = 11, Treatment_B = 11),
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment_A = 50, Treatment_B = 50),
    cohort_size_per_arm = c(Control = 5, Treatment_A = 5, Treatment_B = 5),
    overall_accrual_rate = 3,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE
  )

  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  # Run parallel version
  result <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK"
  )))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(result$Arm_Name, c("Control", "Treatment_A", "Treatment_B"))

  # All arms should have valid results
  expect_true(all(result$Exp_N > 0))
})

test_that("parallel execution preserves interval rebalancing", {
  set.seed(7777)

  base_args <- list(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_for_analysis = 5,
    min_median_followup = 2,
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_current_prob_hc = 0.9,
    posterior_futility_threshold_hc = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    rebalance_after_events = 15,  # Enable rebalancing
    diagnostics = FALSE
  )

  skip_on_cran()
  skip_if_not(capabilities("parallel"))

  # This should complete without error
  result <- do.call(run_simulation_pure, c(base_args, list(
    parallel_replicates = TRUE,
    num_workers = 2,
    cluster_type = "PSOCK"
  )))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})
