# Regression guard for the F1 sign bug (commit e3b1f40).
#
# The between-arm (BA) futility rule must use the *benefit-side* margin
#   Pr(HR >= 1 - delta_HR) >= p_F      (Web Appendix A.3.4)
# i.e. the buffered comparator is the minimal-worthwhile-benefit hazard
# ratio HR = 1 - delta_HR, which lies BELOW the null HR = 1.
# The old bug used the harm-side comparator HR = 1 + delta_HR.
#
# calculate_current_probs_vs_ref() computes its posterior draws internally
# via sample_vs_ref_medians(), so we mock that binding to inject a fixed,
# known logHR vector and then assert the futility comparator direction.

test_that("BA futility margin uses the benefit-side comparator HR = 1 - delta_HR", {
  # Fixed logHR draws chosen to straddle both candidate comparators:
  #   log(0.90) = -0.10536   (benefit side, correct)
  #   log(1.10) =  0.09531   (harm side, the old bug)
  # Several draws fall strictly between these two thresholds, so the
  # benefit-side and harm-side futility probabilities genuinely differ.
  logHR <- c(-0.30, -0.20, -0.08, -0.02, 0.00, 0.03, 0.08, 0.20, 0.40, 0.60)
  delta <- 0.10

  args <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = delta,
    num_posterior_draws = length(logHR)
  )

  testthat::local_mocked_bindings(
    sample_vs_ref_medians = function(...) list(logHR = logHR),
    .package = "evolveTrial"
  )

  res <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args, ctrl_cache = NULL
  )

  benefit_side <- mean(logHR >= log(1 - delta))
  harm_side    <- mean(logHR >= log(1 + delta))

  # The two comparators must actually differ, otherwise the test pins nothing.
  expect_false(isTRUE(all.equal(benefit_side, harm_side)))

  # Assertion: futility probability equals the benefit-side comparator.
  expect_equal(res$pr_fut, benefit_side)

  # Negative control: it must NOT equal the harm-side comparator (old F1 bug).
  expect_false(isTRUE(all.equal(res$pr_fut, harm_side)))

  # Efficacy probability is unaffected by the margin: Pr(HR < 1).
  expect_equal(res$pr_eff, mean(logHR < 0))
})

test_that("BA futility margin defaults to the null comparator HR = 1 when delta is NULL or 0", {
  logHR <- c(-0.30, -0.20, -0.08, -0.02, 0.00, 0.03, 0.08, 0.20, 0.40, 0.60)
  null_comparator <- mean(logHR >= 0)  # Pr(HR >= 1)

  testthat::local_mocked_bindings(
    sample_vs_ref_medians = function(...) list(logHR = logHR),
    .package = "evolveTrial"
  )

  # delta_HR absent (NULL): comparator is the null HR = 1.
  args_null <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = NULL,
    num_posterior_draws = length(logHR)
  )
  res_null <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args_null, ctrl_cache = NULL
  )
  expect_equal(res_null$pr_fut, null_comparator)

  # delta_HR = 0: log1p(-0) = 0, so the comparator is also the null HR = 1.
  args_zero <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = 0,
    num_posterior_draws = length(logHR)
  )
  res_zero <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args_zero, ctrl_cache = NULL
  )
  expect_equal(res_zero$pr_fut, null_comparator)
})
