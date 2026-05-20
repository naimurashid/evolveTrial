# Regression guard for D2: the rule_variant system must actually buffer the
# between-arm (BA) efficacy comparator, not just futility.
#
# The BA decision compares HR (efficacy = HR < 1). delta_HR buffers the
# comparator toward benefit, HR = 1 - delta_HR (below the null HR = 1).
#
#   rule_variant   BA efficacy comparator        BA futility comparator
#   current        Pr(HR < 1)                    Pr(HR >= 1 - delta_HR)
#   symmetric      Pr(HR < 1 - delta_HR)         Pr(HR >= 1 - delta_HR)
#   no_margin      Pr(HR < 1)                    Pr(HR >= 1)
#
# Before this fix `symmetric` was a silent no-op for BA/multi-arm: the
# efficacy comparator was hardcoded as Pr(HR < 1). The efficacy-side margin
# `compare_arms_hr_eff_margin` (NULL except under `symmetric`) now buffers it.
#
# calculate_current_probs_vs_ref() draws its posterior internally via
# sample_vs_ref_medians(), so we mock that binding to inject a fixed logHR.

# Fixed logHR draws chosen so the buffered efficacy comparator log(0.90)
# genuinely differs from the unbuffered HR < 1: several draws fall strictly
# between log(0.90) = -0.10536 and 0, so Pr(HR < 1) != Pr(HR < 0.90).
logHR_fixture <- c(-0.30, -0.20, -0.08, -0.02, 0.00, 0.03, 0.08, 0.20, 0.40, 0.60)

test_that("current: efficacy unbuffered (HR < 1), futility benefit-side", {
  delta <- 0.10
  args <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = delta,        # futility margin
    compare_arms_hr_eff_margin = NULL,     # current: no efficacy buffering
    num_posterior_draws = length(logHR_fixture)
  )

  testthat::local_mocked_bindings(
    sample_vs_ref_medians = function(...) list(logHR = logHR_fixture),
    .package = "evolveTrial"
  )

  res <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args, ctrl_cache = NULL
  )

  expect_equal(res$pr_eff, mean(logHR_fixture < 0))
  expect_equal(res$pr_fut, mean(logHR_fixture >= log(1 - delta)))
})

test_that("symmetric: efficacy buffered (HR < 1 - delta), differs from current", {
  delta <- 0.10
  args <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = delta,        # futility margin
    compare_arms_hr_eff_margin = delta,    # symmetric: buffer efficacy too
    num_posterior_draws = length(logHR_fixture)
  )

  testthat::local_mocked_bindings(
    sample_vs_ref_medians = function(...) list(logHR = logHR_fixture),
    .package = "evolveTrial"
  )

  res <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args, ctrl_cache = NULL
  )

  buffered_eff   <- mean(logHR_fixture < log(1 - delta))
  unbuffered_eff <- mean(logHR_fixture < 0)

  # The buffered and unbuffered efficacy comparators must genuinely differ,
  # otherwise the test pins nothing.
  expect_false(isTRUE(all.equal(buffered_eff, unbuffered_eff)))

  # Assertion that proves symmetric is no longer a no-op for BA/multi-arm.
  expect_equal(res$pr_eff, buffered_eff)
  expect_false(isTRUE(all.equal(res$pr_eff, unbuffered_eff)))

  # Futility comparator is the benefit-side margin (unchanged from current).
  expect_equal(res$pr_fut, mean(logHR_fixture >= log(1 - delta)))
})

test_that("no_margin: both comparators at the null HR = 1", {
  args <- list(
    use_ph_model_vs_ref = TRUE,
    compare_arms_hr_margin = 0,            # no_margin: delta = 0
    compare_arms_hr_eff_margin = NULL,
    num_posterior_draws = length(logHR_fixture)
  )

  testthat::local_mocked_bindings(
    sample_vs_ref_medians = function(...) list(logHR = logHR_fixture),
    .package = "evolveTrial"
  )

  res <- evolveTrial:::calculate_current_probs_vs_ref(
    slCtrl = NULL, slTrt = NULL, args = args, ctrl_cache = NULL
  )

  expect_equal(res$pr_eff, mean(logHR_fixture < 0))
  expect_equal(res$pr_fut, mean(logHR_fixture >= 0))
})
