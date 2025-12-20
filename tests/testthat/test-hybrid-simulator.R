# test-hybrid-simulator.R
# Unit tests for hybrid single-arm to between-arm trial simulator
#
# These tests rely on the package being loaded via testthat.R
# which calls library(evolveTrial) or devtools::load_all()

# ==============================================================================
# TEST SETUP
# ==============================================================================

# Helper to create test configuration
create_test_config <- function() {
  list(
    theta = create_hybrid_theta(
      eff_sa = 0.90,
      fut_sa = 0.10,
      hr_threshold_sa = 0.80,
      ev_sa = 10,
      nmax_sa = 40,
      conversion_trigger = "any_single_success",
      pp_go = 0.70,
      pp_nogo = 0.20,
      ss_method = "posterior",  # Use faster method for tests
      max_additional_n = 60,
      eff_ba = 0.975,
      fut_ba = 0.05,
      ev_ba = 15,
      nmax_ba = 80,
      futility_action = "drop_arm",
      prior_strength = 0.5,
      n_outer = 100,  # Reduced for speed
      n_inner = 50
    ),
    base_args = list(
      n_intervals = 8,
      interval_cutpoints_sim = seq(0, 24, length.out = 9),
      overall_accrual_rate = 2.0,
      max_follow_up_sim = 24,
      max_trial_time = 72,
      prior_alpha_params_model = rep(0.5, 8),
      prior_beta_params_model = rep(0.5, 8)
    ),
    scenario_params = list(
      historical_median = 6,
      ref_median = 7.5,
      exp_median = 9,
      weibull_shape = 1.0
    )
  )
}

# ==============================================================================
# PARAMETER VALIDATION TESTS
# ==============================================================================

test_that("create_hybrid_theta validates conversion trigger", {
  expect_error(
    create_hybrid_theta(conversion_trigger = "invalid"),
    "should be one of"
  )
})

test_that("create_hybrid_theta validates SSR method", {
  expect_error(
    create_hybrid_theta(ss_method = "invalid"),
    "should be one of"
  )
})

test_that("create_hybrid_theta validates futility action", {
  expect_error(
    create_hybrid_theta(futility_action = "invalid"),
    "should be one of"
  )
})

test_that("create_hybrid_theta has correct class", {
  theta <- create_hybrid_theta()
  expect_s3_class(theta, "hybrid_theta")
})

test_that("create_hybrid_theta uses defaults correctly", {
  theta <- create_hybrid_theta()
  expect_equal(theta$eff_sa, 0.90)
  expect_equal(theta$fut_sa, 0.10)
  expect_equal(theta$pp_go, 0.70)
  expect_equal(theta$pp_nogo, 0.20)
  expect_equal(theta$eff_ba, 0.975)
})

# ==============================================================================
# STATE INITIALIZATION TESTS
# ==============================================================================

test_that("create_hybrid_state initializes correctly", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  expect_s3_class(state, "hybrid_trial_state")
  expect_equal(state$current_state, "STATE_SINGLE")
  expect_equal(length(state$arm_names), 2)
  expect_true("Reference" %in% state$active_arms)
  expect_true("Experimental" %in% state$active_arms)
  expect_equal(length(state$dropped_arms), 0)
})

test_that("state has correct posterior structure", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  # Check posterior parameters are initialized
  expect_length(state$posterior_a$Reference, config$base_args$n_intervals)
  expect_length(state$posterior_b$Reference, config$base_args$n_intervals)
  expect_length(state$posterior_a$Experimental, config$base_args$n_intervals)
})

test_that("state has correct registry structure", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  # Registries should be empty data frames
  expect_equal(nrow(state$registries$Reference), 0)
  expect_equal(nrow(state$registries$Experimental), 0)
  expect_true("patient_id" %in% names(state$registries$Reference))
  expect_true("event_time" %in% names(state$registries$Reference))
})

# ==============================================================================
# CONVERSION TRIGGER TESTS
# ==============================================================================

test_that("any_single_success trigger works", {
  config <- create_test_config()
  config$theta$conversion_trigger <- "any_single_success"

  # One arm has efficacy
  sa_efficacy <- c(Reference = FALSE, Experimental = TRUE)
  active_arms <- c("Reference", "Experimental")

  result <- evaluate_conversion_trigger(
    sa_efficacy_reached = sa_efficacy,
    active_arms = active_arms,
    trigger = "any_single_success"
  )

  expect_true(result$triggered)
  expect_equal(result$efficacy_count, 1)
})

test_that("all_single_success trigger requires all arms", {
  # Only one arm has efficacy
  sa_efficacy <- c(Reference = FALSE, Experimental = TRUE)
  active_arms <- c("Reference", "Experimental")

  result <- evaluate_conversion_trigger(
    sa_efficacy_reached = sa_efficacy,
    active_arms = active_arms,
    trigger = "all_single_success"
  )

  expect_false(result$triggered)

  # Both arms have efficacy
  sa_efficacy <- c(Reference = TRUE, Experimental = TRUE)

  result <- evaluate_conversion_trigger(
    sa_efficacy_reached = sa_efficacy,
    active_arms = active_arms,
    trigger = "all_single_success"
  )

  expect_true(result$triggered)
})

test_that("k_of_K trigger works with k=1", {
  sa_efficacy <- c(Reference = TRUE, Experimental = FALSE)
  active_arms <- c("Reference", "Experimental")

  result <- evaluate_conversion_trigger(
    sa_efficacy_reached = sa_efficacy,
    active_arms = active_arms,
    trigger = "k_of_K",
    k_required = 1
  )

  expect_true(result$triggered)
})

# ==============================================================================
# DECISION RULE TESTS
# ==============================================================================

test_that("SA efficacy decision works", {
  result <- evaluate_sa_efficacy(
    p_single = 0.92,
    eff_threshold = 0.90,
    min_events = 15,
    current_events = 20
  )

  expect_true(result$efficacy)
})

test_that("SA efficacy requires minimum events", {
  result <- evaluate_sa_efficacy(
    p_single = 0.92,
    eff_threshold = 0.90,
    min_events = 15,
    current_events = 10  # Below minimum
  )

  expect_false(result$efficacy)
})

test_that("SA futility decision works", {
  result <- evaluate_sa_futility(
    p_single = 0.08,
    fut_threshold = 0.10,
    min_events = 15,
    current_events = 20
  )

  expect_true(result$futility)
})

test_that("BA efficacy decision works", {
  result <- evaluate_ba_efficacy(
    p_between = 0.98,
    eff_threshold = 0.975,
    min_events = 30,
    current_events = 40
  )

  expect_true(result$efficacy)
})

# ==============================================================================
# PP CURVE COMPUTATION TESTS
# ==============================================================================

test_that("PP curve is monotonically increasing", {
  # PP should generally increase with sample size
  # (assuming positive effect)

  config <- create_test_config()

  # Create a state with favorable posterior
  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  # Simulate favorable posterior (lower hazard for experimental)
  state$posterior_a$Experimental <- rep(10, 8)
  state$posterior_b$Experimental <- rep(100, 8)  # Low hazard
  state$posterior_a$Reference <- rep(10, 8)
  state$posterior_b$Reference <- rep(80, 8)  # Higher hazard

  # Compute PP for different N values using posterior method
  pp_curve <- compute_pp_curve_posterior(
    state,
    n_candidates = c(20, 40, 60),
    theta = config$theta,
    base_args = config$base_args
  )

  expect_true(all(diff(pp_curve$pp) >= -0.1))  # Allow small noise
})

# ==============================================================================
# CONVERSION DECISION TESTS
# ==============================================================================

test_that("conversion GO decision when PP >= pp_go", {
  pp_curve <- data.frame(
    n_add = c(20, 40, 60),
    pp = c(0.50, 0.72, 0.85)
  )

  theta <- create_hybrid_theta(pp_go = 0.70, pp_nogo = 0.20)

  result <- make_conversion_decision(pp_curve, theta)

  expect_equal(result$decision, "GO")
  expect_equal(result$n_add, 40)  # First N achieving pp_go
})
test_that("conversion NO_GO decision when max PP < pp_nogo", {
  pp_curve <- data.frame(
    n_add = c(20, 40, 60),
    pp = c(0.05, 0.10, 0.15)
  )

  theta <- create_hybrid_theta(pp_go = 0.70, pp_nogo = 0.20)

  result <- make_conversion_decision(pp_curve, theta)

  expect_equal(result$decision, "NO_GO")
  expect_true(is.na(result$n_add))
})

test_that("conversion AMBIGUOUS decision in middle range", {
  pp_curve <- data.frame(
    n_add = c(20, 40, 60),
    pp = c(0.30, 0.45, 0.55)
  )

  theta <- create_hybrid_theta(pp_go = 0.70, pp_nogo = 0.20)

  result <- make_conversion_decision(pp_curve, theta)

  expect_equal(result$decision, "AMBIGUOUS")
})

# ==============================================================================
# CLOSED-FORM VALIDATION TESTS
# ==============================================================================

test_that("closed-form F-distribution matches Monte Carlo", {
  # Test the between-arm comparison for exponential model

  a_exp <- 20
  b_exp <- 100
  a_ref <- 20
  b_ref <- 80

  validation <- validate_exponential_ba(a_exp, b_exp, a_ref, b_ref, n_samples = 10000)

  # Closed form should be within 2 SEs of MC
  expect_lt(validation$difference, 2 * validation$se_mc)
})

test_that("gamma CDF gives correct SA probability", {

  # P(lambda < c * lambda_hist)

  # Test case: posterior with low hazard should have high P(HR < 0.8)
  # Posterior mean = a/b = 15/200 = 0.075
  # Historical hazard = 0.12
  # Threshold = 0.8 * 0.12 = 0.096
  # With posterior mean 0.075 < threshold 0.096, probability should be high
  a_post <- 15
  b_post <- 200  # Low hazard (rate = 0.075)
  hist_hazard <- 0.12  # log(2)/6 months
  hr_threshold <- 0.80

  threshold <- hr_threshold * hist_hazard
  p_single <- pgamma(threshold, shape = a_post, rate = b_post)

  # Should be reasonably high probability
  expect_gt(p_single, 0.5)
})

# ==============================================================================
# RESULT COMPILATION TESTS
# ==============================================================================

test_that("compile_hybrid_results returns expected fields", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  # Set terminal state
  state$current_state <- "STATE_STOP"
  state$trial_outcome <- "ba_efficacy"
  state$stop_reason <- "Test completion"
  state$sa_efficacy_reached <- c(Reference = TRUE, Experimental = TRUE)
  state$ba_efficacy_reached <- TRUE
  state$conversion_decision <- "GO"
  state$n_enrolled <- c(Reference = 50, Experimental = 50)
  state$n_enrolled_phase1 <- c(Reference = 30, Experimental = 30)

  result <- compile_hybrid_results(state, config$theta)

  expect_true("decision" %in% names(result))
  expect_true("total_n" %in% names(result))
  expect_true("converted" %in% names(result))
  expect_true("sa_efficacy" %in% names(result))
  expect_true("ba_efficacy" %in% names(result))
  expect_equal(result$total_n, 100)
  expect_true(result$converted)
})

# ==============================================================================
# TRIAL OUTCOME CLASSIFICATION TESTS
# ==============================================================================

test_that("is_trial_success correctly identifies BA efficacy", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  state$ba_efficacy_reached <- TRUE

  expect_true(is_trial_success(state, config$theta))
})

test_that("is_trial_success identifies SA efficacy with NO_GO", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  state$sa_efficacy_reached <- c(Reference = TRUE, Experimental = TRUE)
  state$ba_efficacy_reached <- FALSE
  state$conversion_decision <- "NO_GO"

  expect_true(is_trial_success(state, config$theta))
})

test_that("is_trial_success returns FALSE for futility", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  state$sa_efficacy_reached <- c(Reference = FALSE, Experimental = FALSE)
  state$sa_futility_reached <- c(Reference = TRUE, Experimental = TRUE)
  state$ba_efficacy_reached <- FALSE
  state$conversion_decision <- NA

  expect_false(is_trial_success(state, config$theta))
})

# ==============================================================================
# PWE MEDIAN COMPUTATION TESTS
# ==============================================================================

test_that("PWE median computation is correct", {
  # Single interval (exponential) case
  lambda <- 0.1  # hazard rate
  interval_cutpoints <- c(0, 100)  # Wide interval

  median <- compute_pwe_median_survival(lambda, interval_cutpoints)

  # For exponential, median = log(2)/lambda
  expected_median <- log(2) / lambda
  expect_equal(median, expected_median, tolerance = 0.01)
})

test_that("PWE median handles multiple intervals", {
  # Two intervals with different hazards
  lambda <- c(0.1, 0.2)
  interval_cutpoints <- c(0, 5, 100)

  median <- compute_pwe_median_survival(lambda, interval_cutpoints)

  # Median should be positive and finite
  expect_true(median > 0)
  expect_true(is.finite(median))
})

# ==============================================================================
# FUTILITY ACTION TESTS
# ==============================================================================

test_that("futility action stop_trial works", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  result <- apply_futility_action(state, "Experimental", "stop_trial")

  expect_true(result$trial_stopped)
  expect_equal(result$state$current_state, "STATE_STOP")
})

test_that("futility action drop_arm works", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  result <- apply_futility_action(state, "Experimental", "drop_arm")

  expect_true(result$arm_dropped)
  expect_false("Experimental" %in% result$state$active_arms)
  expect_true("Experimental" %in% result$state$dropped_arms)
})

test_that("futility action continue does not modify state", {
  config <- create_test_config()

  state <- create_hybrid_state(
    arm_names = c("Reference", "Experimental"),
    reference_arm = "Reference",
    theta = config$theta,
    base_args = config$base_args
  )

  result <- apply_futility_action(state, "Experimental", "continue")

  expect_false(result$trial_stopped)
  expect_false(result$arm_dropped)
  expect_true("Experimental" %in% result$state$active_arms)
})
