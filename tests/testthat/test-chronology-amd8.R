# Amendment 8, item 1 (evolveTrial 0.0.0.9002): calendar cutoff and
# final-analysis chronology. PI ruling: "Month 72 closes accrual, interim
# analyses and the final analysis. The final analysis is at last enrollment +
# 6 months of follow-up, capped at month 72; interims run only before that
# time; no stop time exceeds 72."
#
# Chronology is captured per trial by amd8_chronology() (helper-chronology.R),
# which traces interim_check() and slice_arm_data_at_time(); the same
# instrumentation runs on 0.0.0.9001, where the first four tests fail.
# Pinned reference values in the last two tests were computed with the
# 0.0.0.9001 engine (commit b6797fe).

CUT <- 72
FU <- 6

check_chronology <- function(ch, cut = CUT, fu = FU) {
  live <- !ch$stopped
  # no enrolment at or after the cutoff
  expect_equal(sum(ch$n_enroll_ge_cut), 0)
  # every interim strictly before T_final (trials not stopped early) and
  # strictly before the cutoff (all trials)
  expect_true(all(ch$max_interim[live] < ch$final_time[live]))
  expect_true(all(ch$max_interim < cut))
  # T_final = min(last enrolment + fu, cut), with the last enrolment before it
  expect_equal(ch$final_time[live], pmin(ch$last_enroll[live] + fu, cut))
  expect_true(all(ch$last_enroll[live] < ch$final_time[live]))
  # N = number of patients enrolled before the final-analysis (or stop) time
  expect_equal(ch$N, ch$n_enroll_before_final)
  # no stop or final-analysis time exceeds the cutoff
  expect_true(all(ch$final_time <= cut))
}

test_that("amd8: slow accrual under H1 (cutoff binds) respects the 72-month chronology", {
  ch <- amd8_chronology(amd8_args(lambda = 0.92, true_median = 9, fu_final = FU), 1:40)
  # nondegenerate: the cutoff binds in most trials
  expect_gt(mean(ch$final_time == CUT & !ch$stopped), 0.5)
  check_chronology(ch)
})

test_that("amd8: slow accrual under H0 respects the 72-month chronology (stop times)", {
  ch <- amd8_chronology(amd8_args(lambda = 0.92, true_median = 6, fu_final = FU), 1:80)
  expect_gt(sum(ch$stopped), 0)
  check_chronology(ch)
})

test_that("amd8: the run_scenarios path enforces the cutoff on aggregate OCs", {
  a <- amd8_args(lambda = 0.92, true_median = 9, fu_final = FU, num_simulations = 100L)
  res <- as.data.frame(run_scenarios(a, list(list()), parallel = FALSE, seed = 7,
                                     return_variance = TRUE))
  expect_lte(res$Exp_Time, CUT)
  # with mean accrual 0.92/month at most ~66 arrivals before month 72; N cap 117
  expect_lt(res$Exp_N, 0.92 * CUT + 4 * sqrt(0.92 * CUT))
})

test_that("amd8: min_follow_up_at_final = 0 drops the look at the moment enrolment ends", {
  # fast accrual, no calendar cap: T_final = last enrolment, so no interim may
  # run at the last enrolment time (0.0.0.9001 ran one there)
  ch <- amd8_chronology(amd8_args(theta = amd8_theta_t50, lambda = 3, true_median = 9,
                                  fu_final = 0, max_trial_time = Inf), 1:30)
  live <- !ch$stopped
  expect_gt(sum(live), 0)
  expect_equal(ch$final_time[live], ch$last_enroll[live])
  expect_true(all(ch$max_interim[live] < ch$final_time[live]))
})

test_that("amd8: max_trial_time = Inf with min_follow_up_at_final > 0 is unchanged from 0.0.0.9001", {
  # Breakage inventory item 2: seeded multi-replicate run, identical to b6797fe.
  bc <- amd8_args(lambda = 1.59, max_trial_time = Inf, fu_final = FU, num_simulations = 200L)
  res <- as.data.frame(run_scenarios(
    bc, list(list(weibull_median_true_arms = c(Experimental = 6)),
             list(weibull_median_true_arms = c(Experimental = 9))),
    parallel = FALSE, seed = 20261004, return_variance = TRUE))
  cols <- c("Type_I_Error_or_Power","PET_Efficacy","PET_Futility","Pr_Final_Efficacy",
            "Pr_Final_Futility","Pr_Final_Inconclusive","Exp_N","Exp_Events","Exp_Time",
            "Var_Exp_N","Var_Exp_Time")
  old <- setNames(list(c(0.029999999999999999, 0.93000000000000005), c(0, 0), c(0.92000000000000004, 
0.065000000000000002), c(0.029999999999999999, 0.93000000000000005
), c(0, 0), c(0.050000000000000003, 0.0050000000000000001), c(51.195, 
112.355), c(40.130000000000003, 93.540000000000006), c(33.489639161583803, 
77.069586523199746), c(3.6480752512562815, 1.6003466080402009
), c(1.6692293179666797, 1.0462310149938689)), cols)
  for (cc in cols) expect_identical(res[[cc]], old[[cc]], label = cc)
})

test_that("amd8: fast accrual where 72 never binds matches 0.0.0.9001 up to T_final", {
  # Breakage inventory item 3. Per-trial seeds (one trial per seed) so that
  # extra 0.0.0.9001 looks at or after T_final, which consume posterior draws,
  # cannot shift the stream of later trials. 0.0.0.9001 ran post-enrolment
  # looks up to month 72 in every trial not stopped early; 0.0.0.9002 runs
  # none at or after T_final. Everything else is equal: N, final time, early
  # stops, final decisions, and the looks strictly before T_final.
  ch <- amd8_chronology(amd8_args(theta = amd8_theta_t50, lambda = 3, true_median = 9,
                                  fu_final = FU), 1:40)
  old <- list(N = c(65, 65, 65, 51, 65, 65, 65, 65, 65, 65, 65, 65, 65, 
65, 65, 65, 65, 65, 33, 65, 41, 37, 65, 65, 65, 40, 48, 65, 65, 
65, 65, 65, 65, 65, 46, 63, 65, 65, 65, 46), final_time = c(24.58053906534467, 
29.08728057694719, 26.596759917166569, 15.159268281221197, 26.585335273983258, 
26.26434908038857, 26.797751772603995, 31.861388389303322, 31.847083883724824, 
26.21912917932314, 24.798331052677156, 27.797639774899316, 27.682634284484742, 
25.401730153853507, 22.170668451809696, 23.862031523778636, 27.257967582947892, 
31.199895432439707, 12.721847025717633, 32.61707641427477, 10.616430515081037, 
10.578362542593011, 27.508429584090106, 30.289469587974665, 31.869070107222697, 
12.64231985604985, 19.313187833660891, 27.263639323384453, 28.566325893838307, 
27.364849039386975, 23.993073793367905, 25.227339095867364, 31.049161956060153, 
26.384456032877114, 12.917082669861657, 21.156149218667238, 24.964494611525652, 
27.624325919620034, 24.784396173441436, 13.042735830370811), 
    stopped = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, 
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, 
    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    TRUE, TRUE, FALSE, FALSE, FALSE, TRUE), final_eff = c(1, 
    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 
    0), final_fut = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0), n_interims_before_final = c(11L, 
    14L, 12L, 6L, 12L, 12L, 12L, 15L, 15L, 12L, 11L, 13L, 13L, 
    12L, 10L, 11L, 13L, 15L, 5L, 15L, 4L, 4L, 13L, 14L, 15L, 
    5L, 8L, 13L, 13L, 13L, 11L, 12L, 14L, 12L, 5L, 9L, 12L, 13L, 
    11L, 5L))
  expect_identical(ch$N, old$N)
  expect_identical(ch$final_time, old$final_time)
  expect_identical(ch$stopped, old$stopped)
  expect_identical(ch$final_eff, old$final_eff)
  expect_identical(ch$final_fut, old$final_fut)
  expect_identical(ch$n_interims_before_final, old$n_interims_before_final)
  expect_identical(ch$n_interims, ch$n_interims_before_final + as.integer(ch$stopped))
})
