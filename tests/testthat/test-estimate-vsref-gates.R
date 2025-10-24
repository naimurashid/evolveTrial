test_that("estimate_vsref_gate_timing reproduces heuristics", {
  args <- list(
    arm_names = c("Doublet", "Triplet"),
    reference_arm_name = "Doublet",
    overall_accrual_rate = 3,
    randomization_probs = c(Doublet = 0.5, Triplet = 0.5),
    max_total_patients_per_arm = c(Doublet = 70, Triplet = 70),
    max_follow_up_sim = 24,
    min_events_per_arm = 8,
    min_median_followup_per_arm = 3,
    min_person_time_frac_per_arm = 0.15
  )

  diag <- estimate_vsref_gate_timing(args)
  expect_named(diag, c("per_arm", "joint_lower_bound"))

  per_arm <- diag$per_arm
  expect_equal(nrow(per_arm), 2L)
  expect_equal(per_arm$AccrualRate, rep(1.5, 2))
  expect_equal(per_arm$MinEvents, rep(8, 2))
  expect_equal(per_arm$TimeForEvents, rep(8 / 1.5, 2))
  expect_equal(per_arm$MinMedianFollowup, rep(3, 2))
  expect_equal(per_arm$TimeForMedianFollowup, rep(6, 2))
  expect_equal(per_arm$MinPersonTimeMonths, rep(70 * 24 * 0.15, 2))
  expect_equal(per_arm$TimeForPersonTime, rep(sqrt(2 * 70 * 24 * 0.15 / 1.5), 2))

  expect_equal(diag$joint_lower_bound, sqrt(2 * 70 * 24 * 0.15 / 1.5))
})
