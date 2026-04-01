# tests/testthat/test-return-variance.R
#
# Tests for return_variance = TRUE in run_simulation_pure() and run_scenarios().
# Covers:
#   - backward compatibility (return_variance = FALSE, the default)
#   - presence and non-negativity of all 10 Var_ columns
#   - analytical correctness of proportion variances
#   - continuous variance decreases with more replicates
#   - compatibility with return_percentiles = TRUE
#   - run_scenarios() forward pass

# ---------------------------------------------------------------------------
# Shared minimal base_args (single-arm, fast)
# ---------------------------------------------------------------------------

base_args <- list(
  num_simulations                        = 100,
  arm_names                              = "Experimental",
  reference_arm_name                     = NULL,
  compare_arms_option                    = FALSE,
  weibull_shape_true_arms                = c(Experimental = 1.3),
  weibull_median_true_arms               = c(Experimental = 6.0),
  null_median_arms                       = c(Experimental = 6.0),
  futility_median_arms                   = c(Experimental = 6.0),
  interval_cutpoints_sim                 = c(0, 8, 16, 24),
  max_follow_up_sim                      = 24,
  censor_max_time_sim                    = 24,
  prior_alpha_params_model               = rep(0.3, 3),
  prior_beta_params_model                = rep(0.3, 3),
  num_posterior_draws                    = 500,
  cohort_size_per_arm                    = 1L,
  max_total_patients_per_arm             = c(Experimental = 60L),
  min_patients_for_analysis              = 0L,
  overall_accrual_rate                   = 2.0,
  randomization_probs                    = c(Experimental = 1.0),
  efficacy_stopping_rule_hc              = TRUE,
  futility_stopping_rule_hc              = TRUE,
  efficacy_threshold_hc_prob             = c(Experimental = 0.95),
  futility_threshold_hc_prob             = c(Experimental = 0.10),
  efficacy_stopping_rule_vs_ref          = FALSE,
  futility_stopping_rule_vs_ref          = FALSE,
  efficacy_threshold_vs_ref_prob         = 0.98,
  futility_threshold_vs_ref_prob         = 0.60,
  use_ph_model_vs_ref                    = FALSE,
  ph_loghr_prior_mean                    = 0,
  ph_loghr_prior_sd                      = 1,
  min_events_hc                          = 10L,
  min_median_followup_hc                 = 0,
  interim_calendar_beat                  = 3.0,
  min_events_per_arm                     = 10L,
  min_median_followup_per_arm            = 3,
  min_person_time_frac_per_arm           = 0,
  min_follow_up_at_final                 = 0,
  rebalance_after_events                 = 50L,
  pred_success_pp_threshold_hc           = 0.90,
  pred_futility_pp_threshold_hc          = 0.20,
  num_posterior_draws_pred               = 500L,
  predictive_fast                        = FALSE,
  median_pfs_success_threshold_arms      = c(Experimental = 6.0),
  median_pfs_futility_threshold_arms     = c(Experimental = 6.0),
  final_success_posterior_prob_threshold = 0.95,
  final_futility_posterior_prob_threshold = 0.10,
  diagnostics                            = FALSE,
  parallel_replicates                    = FALSE
)

# ---------------------------------------------------------------------------
# Backward compatibility: return_variance = FALSE (default)
# ---------------------------------------------------------------------------

test_that("return_variance=FALSE (default) returns plain data.frame with no Var_ columns", {
  set.seed(42L)
  result <- do.call(run_simulation_pure, base_args)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  var_cols <- grep("^Var_", names(result), value = TRUE)
  expect_length(var_cols, 0L)
})

test_that("return_variance=FALSE returns same structure as before the change", {
  set.seed(1L)
  args <- c(base_args, list(return_variance = FALSE))
  result <- do.call(run_simulation_pure, args)

  expect_s3_class(result, "data.frame")
  expected_cols <- c("Arm_Name", "True_Median", "Type_I_Error_or_Power",
                     "PET_Efficacy", "PET_Futility", "Pr_Reach_Max_N",
                     "Pr_Final_Efficacy", "Pr_Final_Futility",
                     "Pr_Final_Inconclusive", "Exp_N", "Exp_Events", "Exp_Time")
  expect_true(all(expected_cols %in% names(result)))
  var_cols <- grep("^Var_", names(result), value = TRUE)
  expect_length(var_cols, 0L)
})

# ---------------------------------------------------------------------------
# Variance columns present and well-formed
# ---------------------------------------------------------------------------

test_that("return_variance=TRUE adds exactly 10 Var_ columns", {
  set.seed(10L)
  args <- c(base_args, list(return_variance = TRUE))
  result <- do.call(run_simulation_pure, args)

  expect_s3_class(result, "data.frame")
  var_cols <- grep("^Var_", names(result), value = TRUE)
  expect_setequal(var_cols, c(
    "Var_Exp_N", "Var_Exp_Events", "Var_Exp_Time",
    "Var_Type_I_Error_or_Power", "Var_PET_Efficacy", "Var_PET_Futility",
    "Var_Pr_Reach_Max_N", "Var_Pr_Final_Efficacy", "Var_Pr_Final_Futility",
    "Var_Pr_Final_Inconclusive"
  ))
})

test_that("all Var_ columns are non-negative", {
  set.seed(11L)
  args <- c(base_args, list(return_variance = TRUE))
  result <- do.call(run_simulation_pure, args)

  var_cols <- grep("^Var_", names(result), value = TRUE)
  for (col in var_cols) {
    expect_true(all(result[[col]] >= 0, na.rm = TRUE),
                label = paste("Var_ column non-negative:", col))
  }
})

test_that("original columns are unchanged when return_variance=TRUE", {
  set.seed(20L); base   <- do.call(run_simulation_pure, base_args)
  set.seed(20L); with_v <- do.call(run_simulation_pure,
                                   utils::modifyList(base_args, list(return_variance = TRUE)))

  core_cols <- c("Type_I_Error_or_Power", "PET_Efficacy", "PET_Futility",
                 "Pr_Reach_Max_N", "Exp_N", "Exp_Events", "Exp_Time")
  for (col in core_cols) {
    expect_equal(base[[col]], with_v[[col]],
                 label = paste("core column unchanged:", col))
  }
})

# ---------------------------------------------------------------------------
# Analytical correctness of proportion variances
# ---------------------------------------------------------------------------

test_that("proportion variances equal p*(1-p)/n analytically", {
  set.seed(30L)
  n <- base_args$num_simulations
  args <- c(base_args, list(return_variance = TRUE))
  result <- do.call(run_simulation_pure, args)

  prop_pairs <- list(
    c("Type_I_Error_or_Power", "Var_Type_I_Error_or_Power"),
    c("PET_Efficacy",          "Var_PET_Efficacy"),
    c("PET_Futility",          "Var_PET_Futility"),
    c("Pr_Reach_Max_N",        "Var_Pr_Reach_Max_N"),
    c("Pr_Final_Efficacy",     "Var_Pr_Final_Efficacy"),
    c("Pr_Final_Futility",     "Var_Pr_Final_Futility"),
    c("Pr_Final_Inconclusive", "Var_Pr_Final_Inconclusive")
  )

  for (pair in prop_pairs) {
    p        <- result[[pair[1]]]
    expected <- p * (1 - p) / n
    actual   <- result[[pair[2]]]
    expect_equal(actual, expected, tolerance = 1e-10,
                 label = paste("analytical variance:", pair[2]))
  }
})

# ---------------------------------------------------------------------------
# Continuous variance decreases with more replicates
# ---------------------------------------------------------------------------

test_that("Var_Exp_N decreases as num_simulations increases", {
  args_small <- utils::modifyList(base_args, list(num_simulations = 100L, return_variance = TRUE))
  args_large <- utils::modifyList(base_args, list(num_simulations = 500L, return_variance = TRUE))

  set.seed(50L); r_small <- do.call(run_simulation_pure, args_small)
  set.seed(50L); r_large <- do.call(run_simulation_pure, args_large)

  # Larger n → smaller variance of the mean
  expect_lt(r_large$Var_Exp_N, r_small$Var_Exp_N)
})

# ---------------------------------------------------------------------------
# Compatibility with return_percentiles = TRUE
# ---------------------------------------------------------------------------

test_that("return_variance and return_percentiles can both be TRUE", {
  set.seed(60L)
  args <- c(base_args, list(return_variance = TRUE, return_percentiles = TRUE))
  result <- do.call(run_simulation_pure, args)

  # return_percentiles wraps output in a list
  expect_true(is.list(result))
  expect_named(result, c("summary", "percentiles"))

  # Variance columns present on summary
  var_cols <- grep("^Var_", names(result$summary), value = TRUE)
  expect_length(var_cols, 10L)
  expect_true(all(vapply(var_cols,
    function(col) all(result$summary[[col]] >= 0, na.rm = TRUE), logical(1))))
})

test_that("return_variance=TRUE with return_percentiles=FALSE returns a data.frame (not list)", {
  set.seed(61L)
  args <- c(base_args, list(return_variance = TRUE, return_percentiles = FALSE))
  result <- do.call(run_simulation_pure, args)

  expect_s3_class(result, "data.frame")
})

# ---------------------------------------------------------------------------
# run_scenarios() forward pass
# ---------------------------------------------------------------------------

test_that("run_scenarios with return_variance=FALSE returns plain data.table (no Var_ cols)", {
  scens <- list(
    list(weibull_median_true_arms = c(Experimental = 6.0)),
    list(weibull_median_true_arms = c(Experimental = 9.0))
  )
  set.seed(70L)
  result <- run_scenarios(base_args, scens, parallel = FALSE, seed = 70L)

  expect_true(data.table::is.data.table(result) || is.data.frame(result))
  expect_equal(nrow(result), 2L)
  var_cols <- grep("^Var_", names(result), value = TRUE)
  expect_length(var_cols, 0L)
})

test_that("run_scenarios with return_variance=TRUE adds Var_ columns for all scenarios", {
  scens <- list(
    list(weibull_median_true_arms = c(Experimental = 6.0)),
    list(weibull_median_true_arms = c(Experimental = 9.0))
  )
  set.seed(71L)
  result <- run_scenarios(base_args, scens, parallel = FALSE,
                          seed = 71L, return_variance = TRUE)

  expect_equal(nrow(result), 2L)
  var_cols <- grep("^Var_", names(result), value = TRUE)
  expect_length(var_cols, 10L)
  expect_true(all(vapply(var_cols,
    function(col) all(result[[col]] >= 0, na.rm = TRUE), logical(1))))
})

test_that("run_scenarios: scenario column present and Var_ columns consistent across scenarios", {
  scens <- list(
    list(weibull_median_true_arms = c(Experimental = 6.0)),
    list(weibull_median_true_arms = c(Experimental = 9.0))
  )
  set.seed(72L)
  result <- run_scenarios(base_args, scens, parallel = FALSE,
                          seed = 72L, return_variance = TRUE)

  expect_true("scenario" %in% names(result))
  expect_setequal(result$scenario, c(1L, 2L))

  # Both scenarios should have non-negative Var_Exp_N
  expect_true(all(result$Var_Exp_N >= 0))
})
