test_that("vs-reference gates honor positional vectors over full arm list", {
  args <- list(
    arm_names = c("A", "B", "C"),
    reference_arm_name = "A",
    randomization_probs = c(A = 0.4, B = 0.3, C = 0.3),
    min_events_per_arm = c(5, 1, 2),
    min_median_followup_per_arm = 0,
    min_person_time_frac_per_arm = 0,
    max_total_patients_per_arm = c(A = 50, B = 50, C = 50),
    max_follow_up_sim = 10
  )

  sl_ctrl <- list(metrics = list(
    events_total = 5,
    median_followup = 0,
    person_time_total = 0
  ))
  sl_trt <- list(metrics = list(
    events_total = 2,
    median_followup = 0,
    person_time_total = 0
  ))

  expect_true(evolveTrial:::gates_pass_for_both_arms(
    sl_ctrl, sl_trt, args, trt_name = "C"
  ))
})
