# Tests for Rcpp vs R implementation equivalence
# These tests verify that the C++ implementations produce
# statistically equivalent results to the R implementations

test_that("simulate_pwe_survival_batch_cpp matches R version", {
  skip_if_not_installed("evolveTrial")

  # Test parameters
  n <- 1000
  lambda <- c(0.05, 0.08, 0.03, 0.06)
  interval_cutpoints <- c(0, 6, 12, 18, 24)

  # Run both with same seed
  set.seed(42)
  r_times <- simulate_pwe_survival_batch(n, lambda, interval_cutpoints)

  set.seed(42)
  cpp_times <- simulate_pwe_survival_batch_cpp(n, lambda, interval_cutpoints)

  # Should be identical with same seed (drop dimensions from arma::vec)
  expect_equal(as.vector(cpp_times), as.vector(r_times), tolerance = 1e-10)

  # Test with different seeds - distributions should match
  set.seed(123)
  r_times2 <- simulate_pwe_survival_batch(10000, lambda, interval_cutpoints)

  set.seed(456)
  cpp_times2 <- simulate_pwe_survival_batch_cpp(10000, lambda, interval_cutpoints)

  # Means should be close (within 5%)
  expect_equal(mean(r_times2), mean(cpp_times2), tolerance = 0.05 * mean(r_times2))

  # Medians should be close
  expect_equal(median(r_times2), median(cpp_times2), tolerance = 0.05 * median(r_times2))
})

test_that("simulate_pwe_survival_batch_cpp handles edge cases", {
  skip_if_not_installed("evolveTrial")

  # Zero hazards should return Inf
  lambda_zero <- c(0, 0, 0)
  cutpoints <- c(0, 6, 12, 18)

  result <- simulate_pwe_survival_batch_cpp(100, lambda_zero, cutpoints)
  expect_true(all(is.infinite(result)))

  # Mixed hazards with zero middle interval
  lambda_mixed <- c(0.1, 0, 0.1)
  set.seed(42)
  result_mixed <- simulate_pwe_survival_batch_cpp(1000, lambda_mixed, cutpoints)

  # Should have some finite values
  expect_true(any(is.finite(result_mixed)))

  # Events should only occur in intervals 1 and 3
  finite_times <- result_mixed[is.finite(result_mixed)]
  expect_true(all(finite_times < 6 | finite_times >= 12))
})

test_that("simulate_future_arm_pwe_cpp matches R version", {
  skip_if_not_installed("evolveTrial")

  n_patients <- 50
  lambda <- c(0.05, 0.08, 0.03, 0.06)
  interval_cutpoints <- c(0, 6, 12, 18, 24)
  accrual_rate <- 2.0
  followup <- 12.0

  # Run with same seed
  set.seed(42)
  r_result <- simulate_future_arm_pwe_vectorized(n_patients, lambda, interval_cutpoints,
                                                  accrual_rate, followup)

  set.seed(42)
  cpp_result <- simulate_future_arm_pwe_cpp(n_patients, lambda, interval_cutpoints,
                                             accrual_rate, followup)

  # Should be identical with same seed (drop dimensions from arma::vec)
  expect_equal(as.vector(cpp_result$events), as.vector(r_result$events), tolerance = 1e-10)
  expect_equal(as.vector(cpp_result$exposure), as.vector(r_result$exposure), tolerance = 1e-10)

  # Test with many runs - average should converge
  n_runs <- 100
  r_events_total <- numeric(4)
  cpp_events_total <- numeric(4)

  for (i in 1:n_runs) {
    r_result <- simulate_future_arm_pwe_vectorized(n_patients, lambda, interval_cutpoints,
                                                    accrual_rate, followup)
    r_events_total <- r_events_total + r_result$events

    cpp_result <- simulate_future_arm_pwe_cpp(n_patients, lambda, interval_cutpoints,
                                               accrual_rate, followup)
    cpp_events_total <- cpp_events_total + cpp_result$events
  }

  # Totals should be within 10%
  expect_equal(sum(cpp_events_total), sum(r_events_total),
               tolerance = 0.1 * sum(r_events_total))
})

test_that("compute_ba_posterior_cpp matches R version", {
  skip_if_not_installed("evolveTrial")

  # Reference R implementation (aggregated F-distribution approach)
  compute_ba_posterior_r <- function(a_exp, b_exp, a_ref, b_ref) {
    a_exp_sum <- sum(a_exp)
    b_exp_sum <- sum(b_exp)
    a_ref_sum <- sum(a_ref)
    b_ref_sum <- sum(b_ref)

    pf(
      q = (b_exp_sum / b_ref_sum) * (a_ref_sum / a_exp_sum),
      df1 = 2 * a_exp_sum,
      df2 = 2 * a_ref_sum
    )
  }

  # Single interval (exponential model)
  a_exp <- 10
  b_exp <- 50
  a_ref <- 12
  b_ref <- 60

  r_result <- compute_ba_posterior_r(a_exp, b_exp, a_ref, b_ref)
  cpp_result <- compute_ba_posterior_cpp(a_exp, b_exp, a_ref, b_ref)

  expect_equal(cpp_result, r_result, tolerance = 1e-10)

  # Multi-interval (PWE model)
  a_exp_pwe <- c(5, 8, 3, 6)
  b_exp_pwe <- c(25, 40, 15, 30)
  a_ref_pwe <- c(6, 7, 4, 5)
  b_ref_pwe <- c(30, 35, 20, 25)

  r_result_pwe <- compute_ba_posterior_r(a_exp_pwe, b_exp_pwe, a_ref_pwe, b_ref_pwe)
  cpp_result_pwe <- compute_ba_posterior_cpp(a_exp_pwe, b_exp_pwe, a_ref_pwe, b_ref_pwe)

  expect_equal(cpp_result_pwe, r_result_pwe, tolerance = 1e-10)
})

test_that("compute_pp_predictive_cpp is statistically equivalent to R version", {
  skip_if_not_installed("evolveTrial")

  # This test compares statistical properties, not exact values
  # since both use Monte Carlo sampling

  # Test parameters
  a_exp <- c(5, 8, 3, 6)
  b_exp <- c(25, 40, 15, 30)
  a_ref <- c(6, 7, 4, 5)
  b_ref <- c(30, 35, 20, 25)
  n_add <- 30
  interval_cutpoints <- c(0, 6, 12, 18, 24)
  accrual_rate <- 2.0
  followup <- 12.0
  eff_ba <- 0.95
  pp_go <- 0.7
  pp_nogo <- 0.2

  # Run multiple times to get distribution
  n_reps <- 20
  cpp_results <- numeric(n_reps)
  r_results <- numeric(n_reps)

  # Create mock state for R version
  state <- list(
    active_arms = c("Exp", "Ref"),
    reference_arm = "Ref",
    posterior_a = list(Exp = a_exp, Ref = a_ref),
    posterior_b = list(Exp = b_exp, Ref = b_ref)
  )

  theta <- list(eff_ba = eff_ba, pp_go = pp_go, pp_nogo = pp_nogo)
  base_args <- list(
    interval_cutpoints_sim = interval_cutpoints,
    overall_accrual_rate = accrual_rate,
    max_follow_up_sim = followup
  )
  scenario_params <- list()

  for (i in 1:n_reps) {
    cpp_results[i] <- compute_pp_predictive_cpp(
      a_exp, b_exp, a_ref, b_ref, n_add,
      interval_cutpoints, accrual_rate, followup,
      eff_ba, pp_go, pp_nogo,
      n_outer = 200, use_antithetic = FALSE, use_early_stop = FALSE
    )

    r_results[i] <- compute_pp_predictive_full(
      state, n_add, theta, base_args, scenario_params,
      n_outer = 200, use_antithetic = FALSE
    )
  }

  # Means should be close (within 0.1 absolute)
  expect_equal(mean(cpp_results), mean(r_results), tolerance = 0.1)

  # Standard deviations should be similar
  expect_equal(sd(cpp_results), sd(r_results), tolerance = 0.05)
})

test_that("compute_pp_predictive_cpp with antithetic variates reduces variance", {
  skip_if_not_installed("evolveTrial")

  # Test parameters
  a_exp <- c(5, 8, 3, 6)
  b_exp <- c(25, 40, 15, 30)
  a_ref <- c(6, 7, 4, 5)
  b_ref <- c(30, 35, 20, 25)
  n_add <- 30
  interval_cutpoints <- c(0, 6, 12, 18, 24)
  accrual_rate <- 2.0
  followup <- 12.0
  eff_ba <- 0.95
  pp_go <- 0.7
  pp_nogo <- 0.2

  n_reps <- 30

  # Without antithetic
  results_standard <- numeric(n_reps)
  for (i in 1:n_reps) {
    results_standard[i] <- compute_pp_predictive_cpp(
      a_exp, b_exp, a_ref, b_ref, n_add,
      interval_cutpoints, accrual_rate, followup,
      eff_ba, pp_go, pp_nogo,
      n_outer = 200, use_antithetic = FALSE, use_early_stop = FALSE
    )
  }

  # With antithetic
  results_antithetic <- numeric(n_reps)
  for (i in 1:n_reps) {
    results_antithetic[i] <- compute_pp_predictive_cpp(
      a_exp, b_exp, a_ref, b_ref, n_add,
      interval_cutpoints, accrual_rate, followup,
      eff_ba, pp_go, pp_nogo,
      n_outer = 200, use_antithetic = TRUE, use_early_stop = FALSE
    )
  }

  # Antithetic should have lower or equal variance
  # (may not always be strictly lower due to sampling variation)
  var_standard <- var(results_standard)
  var_antithetic <- var(results_antithetic)

  # Antithetic variance should be at most 120% of standard (allowing for noise)
  expect_true(var_antithetic <= 1.2 * var_standard)

  # Means should be similar
  expect_equal(mean(results_standard), mean(results_antithetic), tolerance = 0.1)
})

test_that("compute_pp_predictive_cpp is fast", {
  skip_if_not_installed("evolveTrial")
  skip_on_cran()  # Skip on CRAN to avoid flaky timing-based tests

  # Test parameters
  a_exp <- c(5, 8, 3, 6)
  b_exp <- c(25, 40, 15, 30)
  a_ref <- c(6, 7, 4, 5)
  b_ref <- c(30, 35, 20, 25)
  n_add <- 30
  interval_cutpoints <- c(0, 6, 12, 18, 24)
  accrual_rate <- 2.0
  followup <- 12.0
  eff_ba <- 0.95
  pp_go <- 0.7
  pp_nogo <- 0.2

  # Time C++ version
  cpp_time <- system.time({
    for (i in 1:10) {
      compute_pp_predictive_cpp(
        a_exp, b_exp, a_ref, b_ref, n_add,
        interval_cutpoints, accrual_rate, followup,
        eff_ba, pp_go, pp_nogo,
        n_outer = 500, use_antithetic = TRUE, use_early_stop = FALSE
      )
    }
  })["elapsed"]

  cat(sprintf("\nC++ time for 10 runs: %.3f s (%.1f ms per call)\n",
              cpp_time, cpp_time * 100))

  # Use generous threshold (2 seconds) to avoid flaky failures on slow CI
  expect_true(cpp_time < 2.0)
})

# ==============================================================================
# FULL TRIAL SIMULATION RCPP TESTS
# ==============================================================================

test_that("simulate_hybrid_trial_cpp runs without error", {
  skip_if_not_installed("evolveTrial")

  theta <- list(
    eff_sa = 0.90,
    fut_sa = 0.10,
    hr_threshold_sa = 0.80,
    ev_sa = 10,
    nmax_sa = 40,
    conversion_trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    ss_method = "posterior",
    max_additional_n = 60,
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15,
    nmax_ba = 80,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 100
  )

  base_args <- list(
    n_intervals = 4,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 4),
    prior_beta_params_model = rep(0.5, 4)
  )

  # Compute lambda from median: lambda = log(2)/median (exponential approx)
  lambda_hist <- rep(log(2) / 6, 4)      # historical median = 6
  lambda_ref <- rep(log(2) / 7.5, 4)     # ref median = 7.5
  lambda_exp <- rep(log(2) / 9, 4)       # exp median = 9

  scenario_params <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  set.seed(42)
  result <- simulate_hybrid_trial_cpp(theta, base_args, scenario_params)

  expect_true(is.list(result))
  expect_true("trial_outcome" %in% names(result))
  expect_true("total_n" %in% names(result))
  expect_true("converted" %in% names(result))
  expect_true(result$total_n > 0)
})

test_that("run_hybrid_simulations_cpp returns correct structure", {
  skip_if_not_installed("evolveTrial")

  theta <- list(
    eff_sa = 0.90,
    fut_sa = 0.10,
    hr_threshold_sa = 0.80,
    ev_sa = 10,
    nmax_sa = 40,
    conversion_trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    ss_method = "posterior",
    max_additional_n = 60,
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15,
    nmax_ba = 80,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 100
  )

  base_args <- list(
    n_intervals = 4,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 4),
    prior_beta_params_model = rep(0.5, 4)
  )

  # Compute lambda from median
  lambda_hist <- rep(log(2) / 6, 4)
  lambda_ref <- rep(log(2) / 7.5, 4)
  lambda_exp <- rep(log(2) / 9, 4)

  scenario_params <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  set.seed(42)
  results <- run_hybrid_simulations_cpp(100, theta, base_args, scenario_params)

  expect_true(is.data.frame(results))
  expect_equal(nrow(results), 100)
  expect_true("outcome" %in% names(results))
  expect_true("total_n" %in% names(results))
  expect_true("converted" %in% names(results))
})

test_that("compute_hybrid_oc_cpp returns correct operating characteristics", {
  skip_if_not_installed("evolveTrial")

  theta <- list(
    eff_sa = 0.90,
    fut_sa = 0.10,
    hr_threshold_sa = 0.80,
    ev_sa = 10,
    nmax_sa = 40,
    conversion_trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    ss_method = "posterior",
    max_additional_n = 60,
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15,
    nmax_ba = 80,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 100
  )

  base_args <- list(
    n_intervals = 4,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 4),
    prior_beta_params_model = rep(0.5, 4)
  )

  # Alternative scenario (strong effect)
  lambda_hist <- rep(log(2) / 6, 4)
  lambda_ref <- rep(log(2) / 7.5, 4)
  alt_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = rep(log(2) / 12, 4)  # Large effect (median 12)
  )

  # Null scenario (no effect)
  null_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = rep(log(2) / 7.5, 4)  # No effect (same as ref)
  )

  set.seed(42)
  oc <- compute_hybrid_oc_cpp(500, theta, base_args, alt_scenario, null_scenario)

  expect_true(is.list(oc))
  expect_true("power" %in% names(oc))
  expect_true("type1" %in% names(oc))
  expect_true("EN_alt" %in% names(oc))
  expect_true("EN_null" %in% names(oc))
  expect_true("conversion_rate" %in% names(oc))

  # Power should be > 0.5 with strong effect
  expect_gt(oc$power, 0.5)

  # Type I should be < 0.3 (not too high)
  expect_lt(oc$type1, 0.3)
})

test_that("compute_hybrid_oc_cpp is fast", {
  skip_if_not_installed("evolveTrial")
  skip_on_cran()  # Skip on CRAN to avoid flaky timing-based tests

  theta <- list(
    eff_sa = 0.90,
    fut_sa = 0.10,
    hr_threshold_sa = 0.80,
    ev_sa = 10,
    nmax_sa = 40,
    conversion_trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    ss_method = "posterior",
    max_additional_n = 60,
    eff_ba = 0.975,
    fut_ba = 0.05,
    ev_ba = 15,
    nmax_ba = 80,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 100
  )

  base_args <- list(
    n_intervals = 4,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 2.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 4),
    prior_beta_params_model = rep(0.5, 4)
  )

  lambda_hist <- rep(log(2) / 6, 4)
  lambda_ref <- rep(log(2) / 7.5, 4)

  alt_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = rep(log(2) / 12, 4)
  )

  null_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = rep(log(2) / 7.5, 4)
  )

  # Time 1000 simulation OC calculation
  cpp_time <- system.time({
    set.seed(42)
    oc <- compute_hybrid_oc_cpp(1000, theta, base_args, alt_scenario, null_scenario)
  })["elapsed"]

  cat(sprintf("\nRcpp OC time for 1000 sims: %.2f s\n", cpp_time))

  # Should complete in under 15 seconds
  expect_true(cpp_time < 15)
})
