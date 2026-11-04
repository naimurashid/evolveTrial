# HANDOFF-2026-11-04d open item 3: an arm with no final success threshold falls back to Inf and can never succeed at
# the final analysis. run_simulation_pure() warns about it once per R session.

test_that("unset final success threshold warns once per session, and a set threshold does not warn", {
  state <- getFromNamespace(".evolveTrial_warn_state", "evolveTrial")
  old <- state$final_success_unset
  on.exit(state$final_success_unset <- old, add = TRUE)

  run <- function(thr) {
    set.seed(1)
    run_simulation_pure(
      num_simulations = 2, arm_names = "Treatment", reference_arm_name = "Treatment", compare_arms_option = FALSE,
      weibull_shape_true_arms = c(Treatment = 1.3), weibull_median_true_arms = c(Treatment = 9),
      interval_cutpoints_sim = c(0, 6, 12, 24), max_follow_up_sim = 24, censor_max_time_sim = 24,
      randomization_probs = c(Treatment = 1), prior_alpha_params_model = c(0.5, 0.5, 0.5),
      prior_beta_params_model = c(0.5, 0.5, 0.5), num_posterior_draws = 200,
      max_total_patients_per_arm = c(Treatment = 10), cohort_size_per_arm = c(Treatment = 1),
      overall_accrual_rate = 2, median_pfs_success_threshold_arms = thr, progress = FALSE
    )
  }

  state$final_success_unset <- NULL
  expect_silent(run(c(Treatment = 6)))                              # set: no warning, flag not consumed
  expect_warning(run(NULL), "unset for arm\\(s\\) Treatment.*falls back to Inf")
  expect_no_warning(run(NULL))                                      # second unset call: silent
})
