# Tests for multiple trial modes in hybrid trial simulation
# Tests single_arm, between_arm, dual_single_arm, and hybrid modes

test_that("single_arm mode runs correctly", {
  skip_on_cran()

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
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15L,
    nmax_ba = 80L,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6,
    ref_median = 7.5,
    exp_median = 9,
    null_between_improvement = 1.25
  )

  set.seed(123)
  result <- compute_hybrid_oc_rcpp(
    hybrid_theta = hybrid_theta,
    base_args = base_args,
    scenario_params = scenario_params,
    num_simulations = 300,
    seed = 123,
    trial_mode = "single_arm"
  )

  # Check that key metrics are valid
  expect_true(result$power >= 0 && result$power <= 1,
              info = "Power should be between 0 and 1")
  expect_true(result$type1 >= 0 && result$type1 <= 1,
              info = "Type I should be between 0 and 1")
  expect_true(result$EN_alt > 0,
              info = "EN_alt should be positive")

  # single_arm should not have conversion (P_conversion should be 0 or NA)
  expect_true(result$P_conversion == 0 || is.na(result$P_conversion),
              info = "single_arm mode should not have conversion")

  # Verify trial_mode is recorded

expect_equal(result$trial_mode, "single_arm")
})

test_that("between_arm mode runs correctly", {
  skip_on_cran()

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
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15L,
    nmax_ba = 80L,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6,
    ref_median = 7.5,
    exp_median = 9,
    null_between_improvement = 1.25
  )

  set.seed(456)
  result <- compute_hybrid_oc_rcpp(
    hybrid_theta = hybrid_theta,
    base_args = base_args,
    scenario_params = scenario_params,
    num_simulations = 300,
    seed = 456,
    trial_mode = "between_arm"
  )

  # Check that key metrics are valid
  expect_true(result$power >= 0 && result$power <= 1,
              info = "Power should be between 0 and 1")
  expect_true(result$type1 >= 0 && result$type1 <= 1,
              info = "Type I should be between 0 and 1")
  expect_true(result$EN_alt > 0,
              info = "EN_alt should be positive")

  # between_arm should skip SA phase and go directly to BA
  # P_conversion should be 1 (immediate GO) or handled differently
  expect_equal(result$trial_mode, "between_arm")
})

test_that("dual_single_arm mode reports both arms", {
  skip_on_cran()

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
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15L,
    nmax_ba = 80L,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6,
    ref_median = 7.5,
    exp_median = 9,
    null_between_improvement = 1.25
  )

  set.seed(789)
  result <- compute_hybrid_oc_rcpp(
    hybrid_theta = hybrid_theta,
    base_args = base_args,
    scenario_params = scenario_params,
    num_simulations = 300,
    seed = 789,
    trial_mode = "dual_single_arm"
  )

  # Check that per-arm metrics are present
  expect_true(!is.null(result$power_exp),
              info = "power_exp should be present for dual_single_arm")
  expect_true(!is.null(result$power_ref),
              info = "power_ref should be present for dual_single_arm")
  expect_true(!is.null(result$type1_exp),
              info = "type1_exp should be present for dual_single_arm")
  expect_true(!is.null(result$type1_ref),
              info = "type1_ref should be present for dual_single_arm")

  # Check that per-arm metrics are valid
  expect_true(result$power_exp >= 0 && result$power_exp <= 1,
              info = "power_exp should be between 0 and 1")
  expect_true(result$power_ref >= 0 && result$power_ref <= 1,
              info = "power_ref should be between 0 and 1")
  expect_true(result$type1_exp >= 0 && result$type1_exp <= 1,
              info = "type1_exp should be between 0 and 1")
  expect_true(result$type1_ref >= 0 && result$type1_ref <= 1,
              info = "type1_ref should be between 0 and 1")

  # Verify trial_mode is recorded
  expect_equal(result$trial_mode, "dual_single_arm")

  # Check that note is present for dual_single_arm
  expect_true(!is.null(result$note),
              info = "Note should be present for dual_single_arm mode")
})

test_that("hybrid mode still works correctly (regression test)", {
  skip_on_cran()

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
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15L,
    nmax_ba = 80L,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6,
    ref_median = 7.5,
    exp_median = 9,
    null_between_improvement = 1.25
  )

  set.seed(111)
  result <- compute_hybrid_oc_rcpp(
    hybrid_theta = hybrid_theta,
    base_args = base_args,
    scenario_params = scenario_params,
    num_simulations = 300,
    seed = 111,
    trial_mode = "hybrid"
  )

  # Check that key metrics are valid
  expect_true(result$power >= 0 && result$power <= 1,
              info = "Power should be between 0 and 1")
  expect_true(result$type1 >= 0 && result$type1 <= 1,
              info = "Type I should be between 0 and 1")
  expect_true(result$EN_alt > 0,
              info = "EN_alt should be positive")

  # Hybrid mode should have conversion
  expect_true(result$P_conversion >= 0 && result$P_conversion <= 1,
              info = "P_conversion should be between 0 and 1 for hybrid mode")

  # Verify trial_mode is recorded
  expect_equal(result$trial_mode, "hybrid")
})

test_that("efficacy_method and futility_method parameters are validated", {
  skip_on_cran()

  hybrid_theta <- list(
    eff_sa = 0.90, fut_sa = 0.10, hr_threshold_sa = 0.80,
    ev_sa = 15L, nmax_sa = 50L,
    conversion_trigger = "any_single_success",
    pp_go = 0.70, pp_nogo = 0.20, ss_method = "posterior",
    max_additional_n = 60L,
    eff_ba = 0.975, fut_ba = 0.05, ev_ba = 15L, nmax_ba = 80L,
    futility_action = "drop_arm", prior_strength = 0.5, n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6, ref_median = 7.5, exp_median = 9
  )

  # Test invalid efficacy_method
  expect_error(
    compute_hybrid_oc_rcpp(
      hybrid_theta = hybrid_theta,
      base_args = base_args,
      scenario_params = scenario_params,
      num_simulations = 10,
      efficacy_method = "invalid"
    ),
    "posterior.*predictive"
  )

  # Test invalid futility_method
  expect_error(
    compute_hybrid_oc_rcpp(
      hybrid_theta = hybrid_theta,
      base_args = base_args,
      scenario_params = scenario_params,
      num_simulations = 10,
      futility_method = "invalid"
    ),
    "posterior.*predictive"
  )

  # Test invalid trial_mode
  expect_error(
    compute_hybrid_oc_rcpp(
      hybrid_theta = hybrid_theta,
      base_args = base_args,
      scenario_params = scenario_params,
      num_simulations = 10,
      trial_mode = "invalid_mode"
    ),
    "trial_mode"
  )
})

test_that("single_arm mode skips BA phase entirely", {
  skip_on_cran()

  hybrid_theta <- list(
    eff_sa = 0.90, fut_sa = 0.10, hr_threshold_sa = 0.80,
    ev_sa = 15L, nmax_sa = 50L,
    conversion_trigger = "any_single_success",
    pp_go = 0.70, pp_nogo = 0.20, ss_method = "posterior",
    max_additional_n = 60L,
    eff_ba = 0.975, fut_ba = 0.05, ev_ba = 15L, nmax_ba = 80L,
    futility_action = "drop_arm", prior_strength = 0.5, n_outer = 200L
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

  scenario_params <- list(
    historical_median = 6, ref_median = 7.5, exp_median = 9
  )

  # Run raw simulations to check outcomes
  n_intervals <- 4
  lambda_hist <- rep(log(2) / 6, n_intervals)
  lambda_ref <- rep(log(2) / 7.5, n_intervals)
  lambda_exp <- rep(log(2) / 9, n_intervals)

  theta_cpp <- list(
    eff_sa = 0.90, fut_sa = 0.10, hr_threshold_sa = 0.80,
    ev_sa = 15L, nmax_sa = 50L,
    conversion_trigger = "any_single_success",
    pp_go = 0.70, pp_nogo = 0.20, ss_method = "posterior",
    max_additional_n = 60L,
    eff_ba = 0.975, fut_ba = 0.05, ev_ba = 15L, nmax_ba = 80L,
    futility_action = "drop_arm", prior_strength = 0.5, n_outer = 200L,
    trial_mode = "single_arm",
    efficacy_method = "posterior",
    futility_method = "posterior"
  )

  base_args_cpp <- list(
    n_intervals = n_intervals,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, n_intervals),
    prior_beta_params_model = rep(0.5, n_intervals)
  )

  scenario_cpp <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  set.seed(222)
  results <- run_hybrid_simulations_cpp(
    n_sim = 200L,
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = scenario_cpp
  )

  # In single_arm mode, converted should always be FALSE
  expect_true(all(results$converted == FALSE),
              info = "single_arm mode should never convert to BA phase")

  # BA efficacy should always be FALSE
  expect_true(all(results$ba_efficacy == FALSE),
              info = "single_arm mode should never have BA efficacy")
})

test_that("between_arm mode skips SA phase entirely", {
  skip_on_cran()

  n_intervals <- 4
  lambda_hist <- rep(log(2) / 6, n_intervals)
  lambda_ref <- rep(log(2) / 7.5, n_intervals)
  lambda_exp <- rep(log(2) / 9, n_intervals)

  theta_cpp <- list(
    eff_sa = 0.90, fut_sa = 0.10, hr_threshold_sa = 0.80,
    ev_sa = 15L, nmax_sa = 50L,
    conversion_trigger = "any_single_success",
    pp_go = 0.70, pp_nogo = 0.20, ss_method = "posterior",
    max_additional_n = 60L,
    eff_ba = 0.975, fut_ba = 0.05, ev_ba = 15L, nmax_ba = 80L,
    futility_action = "drop_arm", prior_strength = 0.5, n_outer = 200L,
    trial_mode = "between_arm",
    efficacy_method = "posterior",
    futility_method = "posterior"
  )

  base_args_cpp <- list(
    n_intervals = n_intervals,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, n_intervals),
    prior_beta_params_model = rep(0.5, n_intervals)
  )

  scenario_cpp <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  set.seed(333)
  results <- run_hybrid_simulations_cpp(
    n_sim = 200L,
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = scenario_cpp
  )

  # In between_arm mode, SA efficacy should always be FALSE (skipped)
  expect_true(all(results$sa_efficacy == FALSE),
              info = "between_arm mode should skip SA phase")

  # Converted should be TRUE (immediate GO) or handled as direct BA
  expect_true(all(results$converted == TRUE),
              info = "between_arm mode should always 'convert' to BA phase")
})

test_that("compute_oc_lambda works with direct hazard input", {
  skip_on_cran()

  base_args <- list(
    interval_cutpoints_sim = c(0, 6, 12, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 12,
    look_interval = 3.0,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 3),
    prior_beta_params_model = rep(0.5, 3)
  )

  theta <- list(
    eff_sa = 0.90, fut_sa = 0.10,
    eff_ba = 0.90, fut_ba = 0.10,
    ev_sa = 10, ev_ba = 20,
    nmax_sa = 40, nmax_ba = 80,
    hr_threshold_sa = 0.7,
    pp_go = 0.5, pp_nogo = 0.3
  )

  set.seed(444)
  result <- compute_oc_lambda(
    theta = theta,
    base_args = base_args,
    lambda_exp = c(0.04, 0.05, 0.06),
    lambda_ref = c(0.08, 0.09, 0.10),
    lambda_hist = c(0.05, 0.06, 0.07),
    num_simulations = 200,
    trial_mode = "hybrid"
  )

  expect_true(result$power >= 0 && result$power <= 1,
              info = "Power should be between 0 and 1")
  expect_true(result$type1 >= 0 && result$type1 <= 1,
              info = "Type I should be between 0 and 1")
  expect_equal(result$trial_mode, "hybrid")
})

test_that("compute_oc_lambda validates lambda lengths", {
  skip_on_cran()

  base_args <- list(
    interval_cutpoints_sim = c(0, 6, 12, 24),  # 3 intervals
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 12,
    prior_alpha_params_model = rep(0.5, 3),
    prior_beta_params_model = rep(0.5, 3)
  )

  theta <- list(
    eff_sa = 0.90, fut_sa = 0.10,
    eff_ba = 0.90, fut_ba = 0.10
  )

  # Wrong length lambda_exp (2 instead of 3)
  expect_error(
    compute_oc_lambda(
      theta = theta,
      base_args = base_args,
      lambda_exp = c(0.04, 0.05),  # Wrong length
      lambda_ref = c(0.08, 0.09, 0.10),
      lambda_hist = c(0.05, 0.06, 0.07),
      num_simulations = 10
    ),
    "lambda_exp must have length"
  )
})

test_that("predictive probability methods work for efficacy", {
  skip_on_cran()

  base_args <- list(
    interval_cutpoints_sim = c(0, 6, 12, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 12,
    look_interval = 3.0,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 3),
    prior_beta_params_model = rep(0.5, 3)
  )

  theta <- list(
    eff_sa = 0.90, fut_sa = 0.10,
    eff_ba = 0.90, fut_ba = 0.10,
    ev_sa = 10, ev_ba = 20,
    nmax_sa = 40, nmax_ba = 80,
    hr_threshold_sa = 0.7,
    pp_go = 0.5, pp_nogo = 0.3
  )

  # Test with predictive for efficacy, posterior for futility
  set.seed(555)
  result_pp <- compute_oc_lambda(
    theta = theta,
    base_args = base_args,
    lambda_exp = c(0.04, 0.05, 0.06),
    lambda_ref = c(0.08, 0.09, 0.10),
    lambda_hist = c(0.05, 0.06, 0.07),
    num_simulations = 100,
    trial_mode = "single_arm",
    efficacy_method = "predictive",
    futility_method = "posterior"
  )

  expect_equal(result_pp$efficacy_method, "predictive")
  expect_equal(result_pp$futility_method, "posterior")
  expect_true(result_pp$power >= 0 && result_pp$power <= 1)
})

test_that("futility_action parameter affects trial behavior", {
  skip_on_cran()

  base_args <- list(
    interval_cutpoints_sim = c(0, 6, 12, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 12,
    look_interval = 3.0,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 3),
    prior_beta_params_model = rep(0.5, 3)
  )

  # Scenario with easy futility trigger (exp similar to hist)
  theta_drop <- list(
    eff_sa = 0.99, fut_sa = 0.50,  # Easy futility
    eff_ba = 0.90, fut_ba = 0.10,
    ev_sa = 5, ev_ba = 10,
    nmax_sa = 30, nmax_ba = 60,
    hr_threshold_sa = 0.7,
    pp_go = 0.5, pp_nogo = 0.3,
    futility_action = "drop_arm"
  )

  theta_stop <- theta_drop
  theta_stop$futility_action <- "stop_trial"

  theta_continue <- theta_drop
  theta_continue$futility_action <- "continue"

  # Scenario where experimental is worse than historical
  set.seed(666)
  result_drop <- compute_oc_lambda(
    theta = theta_drop, base_args = base_args,
    lambda_exp = c(0.10, 0.10, 0.10),  # Worse than historical
    lambda_ref = c(0.05, 0.05, 0.05),
    lambda_hist = c(0.05, 0.05, 0.05),
    num_simulations = 100,
    trial_mode = "single_arm"
  )

  set.seed(666)
  result_continue <- compute_oc_lambda(
    theta = theta_continue, base_args = base_args,
    lambda_exp = c(0.10, 0.10, 0.10),
    lambda_ref = c(0.05, 0.05, 0.05),
    lambda_hist = c(0.05, 0.05, 0.05),
    num_simulations = 100,
    trial_mode = "single_arm"
  )

  # "continue" should allow trial to run longer (higher EN)
  expect_true(result_continue$EN_alt >= result_drop$EN_alt,
              info = "continue futility_action should lead to higher EN than drop_arm")
})
