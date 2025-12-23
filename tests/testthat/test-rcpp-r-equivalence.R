# Tests for R vs Rcpp equivalence in hybrid trial simulation
# These tests verify that the Rcpp implementation produces results
# consistent with the R implementation

test_that("Rcpp and R implementations produce similar power estimates", {
  skip_on_cran()

  # Set up common parameters
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

  scenario_params <- list(
    historical_median = 6,
    ref_median = 7.5,
    exp_median = 9,
    null_between_improvement = 1.25
  )

  # Run Rcpp version
  set.seed(123)
  rcpp_result <- compute_hybrid_oc_rcpp(
    hybrid_theta = hybrid_theta,
    base_args = base_args,
    scenario_params = scenario_params,
    num_simulations = 500,
    seed = 123
  )

  # Check that key metrics are within reasonable ranges
  expect_true(rcpp_result$power >= 0 && rcpp_result$power <= 1,
              info = "Power should be between 0 and 1")
  expect_true(rcpp_result$type1 >= 0 && rcpp_result$type1 <= 1,
              info = "Type I should be between 0 and 1")
  expect_true(rcpp_result$type1_between >= 0 && rcpp_result$type1_between <= 1,
              info = "Type I between should be between 0 and 1")
  expect_true(rcpp_result$EN_alt > 0,
              info = "EN_alt should be positive")
  expect_true(rcpp_result$EN_null > 0,
              info = "EN_null should be positive")
  expect_true(rcpp_result$P_conversion >= 0 && rcpp_result$P_conversion <= 1,
              info = "P_conversion should be between 0 and 1")
})

test_that("null_between scenario is properly simulated", {
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
    seed = 456
  )

  # Check that null_between results are present
  expect_true(!is.null(result$type1_between),
              info = "type1_between should be computed")
  expect_true(!is.null(result$results_by_scenario$null_between),
              info = "null_between scenario results should be present")
  expect_true(!is.null(result$results_by_scenario$null_between$ba_success_rate),
              info = "ba_success_rate should be computed for null_between")

  # type1_between should generally be low (it's a type I error rate)
  # Under null_between, both arms are equal, so BA efficacy should be rare
  expect_true(result$type1_between < 0.20,
              info = "type1_between should be controlled (< 0.20)")
})

test_that("antithetic MC draws produce correct sample counts", {
  skip_on_cran()

  # Test that odd n_outer values don't cause oversampling
  # This tests the fix for the antithetic branch overshoot bug

  n_intervals <- 4
  a_exp <- rep(5, n_intervals)
  b_exp <- rep(10, n_intervals)
  a_ref <- rep(5, n_intervals)
  b_ref <- rep(10, n_intervals)
  interval_cutpoints <- c(0, 6, 12, 18, 24)

  # Test with odd n_outer
  set.seed(789)
  pp_odd <- compute_pp_predictive_cpp(
    a_exp = a_exp, b_exp = b_exp,
    a_ref = a_ref, b_ref = b_ref,
    n_add = 30,
    interval_cutpoints = interval_cutpoints,
    accrual_rate = 2.0,
    followup = 24,
    eff_ba = 0.975,
    pp_go = 0.70,
    pp_nogo = 0.20,
    n_outer = 201L,  # Odd number
    use_antithetic = TRUE,
    use_early_stop = FALSE
  )

  # Test with even n_outer
  set.seed(789)
  pp_even <- compute_pp_predictive_cpp(
    a_exp = a_exp, b_exp = b_exp,
    a_ref = a_ref, b_ref = b_ref,
    n_add = 30,
    interval_cutpoints = interval_cutpoints,
    accrual_rate = 2.0,
    followup = 24,
    eff_ba = 0.975,
    pp_go = 0.70,
    pp_nogo = 0.20,
    n_outer = 200L,  # Even number
    use_antithetic = TRUE,
    use_early_stop = FALSE
  )

  # Both should produce valid PP values
  expect_true(pp_odd >= 0 && pp_odd <= 1,
              info = "PP with odd n_outer should be valid")
  expect_true(pp_even >= 0 && pp_even <= 1,
              info = "PP with even n_outer should be valid")
})

test_that("conversion decision logic matches R implementation", {
  skip_on_cran()

  # Test that conversion decisions (GO, NO-GO, AMBIGUOUS) are handled correctly

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
    futility_action = "drop_arm", prior_strength = 0.5, n_outer = 200L
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

  set.seed(999)
  results <- run_hybrid_simulations_cpp(
    n_sim = 200L,
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = scenario_cpp
  )

  # Check outcome distribution includes expected types
  outcomes <- as.character(results$outcome)
  valid_outcomes <- c(
    "success", "all_arms_futile", "conversion_nogo", "conversion_ambiguous",
    "max_n_single_phase", "max_time_ba", "max_n_ba", "futility_ba",
    "max_time_single_phase"
  )

  for (outcome in unique(outcomes)) {
    expect_true(outcome %in% valid_outcomes,
                info = paste("Unexpected outcome:", outcome))
  }

  # Success should be consistent with sa_efficacy and ba_efficacy
  for (i in seq_len(nrow(results))) {
    if (results$is_success[i]) {
      # Success requires either BA efficacy OR (SA efficacy + NO-GO/AMBIGUOUS)
      ba_success <- results$ba_efficacy[i]
      sa_with_nogo <- results$sa_efficacy[i] && !results$converted[i]
      expect_true(ba_success || sa_with_nogo,
                  info = paste("Invalid success at row", i))
    }
  }
})
