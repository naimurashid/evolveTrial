# test-hybrid-ba-margin.R
# Tests for the D3 fix: PH-model benefit-side delta_HR futility margin in the
# between-arm (Stage 2) decision of the C++ hybrid simulator.
#
# The BA decision now uses the proportional-hazards logHR posterior
# (ph_beta_mode_var_cpp) instead of the independent-arms median-comparison
# model. Even hr_margin = 0 switches the posterior model, so there is no inert
# default to regression-test against; instead we test the *new* behavior:
#   (a) runs complete and return finite, in-range OCs;
#   (b) raising hr_margin moves the BA futility rate in the benefit-side
#       direction (more aggressive futility) for a near-null experimental arm.

# Build theta / base_args / scenario lists directly in the C++-expected format.
# We call run_hybrid_simulations_cpp / compute_hybrid_oc_cpp directly because
# the compute_hybrid_oc_rcpp wrapper rebuilds theta and would drop the new
# hr_margin field.

make_ba_theta <- function(hr_margin = 0.0, hr_eff_margin = 0.0) {
  list(
    eff_sa = 0.90,
    fut_sa = 0.10,
    hr_threshold_sa = 0.80,
    ev_sa = 10L,
    nmax_sa = 40L,
    conversion_trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    ss_method = "posterior",
    max_additional_n = 60L,
    eff_ba = 0.975,
    fut_ba = 0.50,        # p_F role: futility if Pr(HR >= 1 - delta_HR) >= fut_ba
    ev_ba = 15L,
    nmax_ba = 80L,
    futility_action = "drop_arm",
    prior_strength = 0.5,
    n_outer = 100L,
    trial_mode = "between_arm",
    efficacy_method = "posterior",
    futility_method = "posterior",
    hr_margin = hr_margin,
    hr_eff_margin = hr_eff_margin
  )
}

make_ba_base_args <- function() {
  list(
    n_intervals = 4L,
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),
    overall_accrual_rate = 3.0,
    max_follow_up_sim = 24,
    max_trial_time = 72,
    prior_alpha_params_model = rep(0.5, 4),
    prior_beta_params_model = rep(0.5, 4)
  )
}

# Near-null experimental arm: exp barely better than ref. With a benefit-side
# margin the futility comparator HR = 1 - delta_HR moves below 1, so
# Pr(HR >= 1 - delta_HR) grows -> futility becomes more aggressive.
make_near_null_scenario <- function() {
  list(
    lambda_hist = rep(log(2) / 7.5, 4),
    lambda_ref  = rep(log(2) / 7.5, 4),
    lambda_exp  = rep(log(2) / 7.7, 4)  # essentially equal to ref
  )
}

ba_futility_rate <- function(theta, base_args, scenario, n_sim, seed) {
  set.seed(seed)
  res <- run_hybrid_simulations_cpp(
    n_sim = as.integer(n_sim),
    theta_list = theta,
    base_args_list = base_args,
    scenario_params_list = scenario
  )
  mean(res$outcome == "futility_ba")
}

# ==============================================================================
# TEST 1: runs complete and return finite, in-range OCs for both margins
# ==============================================================================

test_that("hybrid BA with PH delta_HR margin returns finite in-range OCs", {
  skip_on_cran()

  base_args <- make_ba_base_args()
  scenario <- make_near_null_scenario()
  # Null scenario: experimental equal to reference.
  null_scenario <- list(
    lambda_hist = rep(log(2) / 7.5, 4),
    lambda_ref  = rep(log(2) / 7.5, 4),
    lambda_exp  = rep(log(2) / 7.5, 4)
  )

  for (m in c(0.0, 0.10)) {
    theta <- make_ba_theta(hr_margin = m)
    set.seed(2026)
    oc <- compute_hybrid_oc_cpp(
      n_sim = 300L,
      theta_list = theta,
      base_args_list = base_args,
      scenario_params_list = scenario,
      null_scenario_list = null_scenario
    )

    expect_true(is.finite(oc$power),
                info = sprintf("power finite (hr_margin = %.2f)", m))
    expect_true(is.finite(oc$type1),
                info = sprintf("type1 finite (hr_margin = %.2f)", m))
    expect_gte(oc$power, 0)
    expect_lte(oc$power, 1)
    expect_gte(oc$type1, 0)
    expect_lte(oc$type1, 1)
    expect_true(is.finite(oc$EN_alt) && oc$EN_alt > 0,
                info = sprintf("EN_alt positive (hr_margin = %.2f)", m))
  }
})

# ==============================================================================
# TEST 2: raising hr_margin makes BA futility more aggressive (benefit-side)
# ==============================================================================

test_that("raising hr_margin increases BA futility rate for a near-null arm", {
  skip_on_cran()

  base_args <- make_ba_base_args()
  scenario <- make_near_null_scenario()

  # With a near-null experimental arm, the PH logHR posterior is centered close
  # to 0. Pr(HR >= 1 - delta_HR) = pnorm(log(1 - delta_HR), mean, sd, lower=F).
  # log(1 - delta_HR) < 0 for delta_HR > 0, so the lower bound of the tail moves
  # left -> the tail probability is strictly larger -> futility (fut_ba gate on
  # pr_fut) fires more often. This is the benefit-side direction.
  fut_rate_0   <- ba_futility_rate(make_ba_theta(hr_margin = 0.00),
                                   base_args, scenario, n_sim = 400, seed = 4242)
  fut_rate_010 <- ba_futility_rate(make_ba_theta(hr_margin = 0.10),
                                   base_args, scenario, n_sim = 400, seed = 4242)

  # The shift should be material, not just noise, for a 0.10 margin.
  expect_gt(fut_rate_010, fut_rate_0)
})

# ==============================================================================
# TEST 3: direct unit check of the PH benefit-side direction
# ==============================================================================
# ph_beta_mode_var_cpp is exported, so we can verify the benefit-side sign
# convention directly without a custom test shim: for a near-null logHR
# posterior, P(HR >= 1 - delta) must increase as delta increases.

test_that("PH posterior tail Pr(HR >= 1 - delta) increases with delta", {
  skip_on_cran()

  # Reference (control) and experimental (treatment) per-interval counts:
  # roughly equal event rates -> logHR posterior centered near 0.
  E_C <- c(8L, 6L, 4L, 2L)
  PT_C <- c(120, 90, 60, 30)
  E_T <- c(7L, 6L, 4L, 2L)
  PT_T <- c(122, 92, 61, 31)
  alpha_prior <- rep(0.5, 4)
  beta_prior  <- rep(0.5, 4)

  bp <- ph_beta_mode_var_cpp(E_C, PT_C, E_T, PT_T,
                             alpha_prior, beta_prior,
                             mu = 0, sigma = 1)
  mean <- bp$mean
  sd <- bp$sd
  expect_true(is.finite(mean) && is.finite(sd) && sd > 0)

  pr_fut <- function(delta) {
    log_fut <- log1p(-min(max(0, delta), 0.95))
    pnorm(log_fut, mean = mean, sd = sd, lower.tail = FALSE)
  }

  # Benefit-side: larger delta_HR -> lower futility comparator HR = 1 - delta
  # -> log_fut more negative -> larger upper-tail probability.
  expect_gt(pr_fut(0.10), pr_fut(0.00))
  expect_gt(pr_fut(0.20), pr_fut(0.10))
})

# ==============================================================================
# TEST 4 (Task 1.4b): conversion-PP forecast uses the PH + delta_HR rule
# ==============================================================================
# The SA->BA conversion-PP forecast now scores each forecast replicate's BA
# outcome with the same PH logHR model as the real Stage-2 BA decision: a
# replicate "succeeds" when pr_eff >= eff_ba. pr_eff depends on the benefit-side
# efficacy comparator HR = 1 - delta_HR_eff. Raising hr_eff_margin makes that
# comparator strictly harder (P(HR < 1 - delta_HR_eff) <= P(HR < 1)), so every
# forecast replicate is less likely to clear eff_ba -> the predictive
# probability at conversion falls -> the GO conversion rate falls. This verifies
# the forecast and the real BA decision share the same rule.

# Moderate-effect experimental arm: clearly better than reference, so the PP
# forecast sits in a sensitive range where a stricter efficacy comparator
# visibly lowers it.
make_moderate_effect_scenario <- function() {
  list(
    lambda_hist = rep(log(2) / 7.5, 4),
    lambda_ref  = rep(log(2) / 7.5, 4),
    lambda_exp  = rep(log(2) / 11.0, 4)  # ~30% lower hazard than reference
  )
}

conversion_pp_summary <- function(theta, base_args, scenario, n_sim, seed) {
  # The conversion-PP forecast is only invoked in the full SA -> conversion ->
  # BA "hybrid" mode. make_ba_theta() uses "between_arm" (which skips the SA
  # phase and the PP forecast entirely), so override the mode here.
  theta$trial_mode <- "hybrid"
  set.seed(seed)
  res <- run_hybrid_simulations_cpp(
    n_sim = as.integer(n_sim),
    theta_list = theta,
    base_args_list = base_args,
    scenario_params_list = scenario
  )
  # pp_conversion is NA when conversion was never PP-evaluated; restrict the
  # mean PP to trials that actually reached the PP-based conversion decision.
  evaluated <- res$pp_conversion[is.finite(res$pp_conversion)]
  list(
    conv_rate = mean(res$converted),
    mean_pp   = if (length(evaluated) > 0) mean(evaluated) else NA_real_,
    n_eval    = length(evaluated)
  )
}

test_that("conversion-PP forecast runs and returns in-range PP for both margins", {
  skip_on_cran()

  base_args <- make_ba_base_args()
  scenario <- make_moderate_effect_scenario()

  for (m in c(0.0, 0.10)) {
    theta <- make_ba_theta(hr_eff_margin = m)
    s <- conversion_pp_summary(theta, base_args, scenario,
                               n_sim = 300, seed = 909 + as.integer(m * 100))
    expect_gte(s$conv_rate, 0)
    expect_lte(s$conv_rate, 1)
    if (!is.na(s$mean_pp)) {
      expect_gte(s$mean_pp, 0)
      expect_lte(s$mean_pp, 1)
    }
  }
})

test_that("raising hr_eff_margin lowers the conversion-PP forecast", {
  skip_on_cran()

  base_args <- make_ba_base_args()
  scenario <- make_moderate_effect_scenario()

  # Same seed -> same SA-phase data and same forecast random draws; only the
  # efficacy comparator changes, isolating the rule effect.
  s0   <- conversion_pp_summary(make_ba_theta(hr_eff_margin = 0.00),
                                base_args, scenario, n_sim = 500, seed = 7171)
  s015 <- conversion_pp_summary(make_ba_theta(hr_eff_margin = 0.15),
                                base_args, scenario, n_sim = 500, seed = 7171)

  # Need PP-evaluated trials in both runs for the comparison to be meaningful.
  skip_if(s0$n_eval < 20 || s015$n_eval < 20,
          "too few PP-evaluated conversions to compare")

  # A stricter benefit-side efficacy comparator can only lower (never raise)
  # each replicate's success probability, so the mean conversion PP must not
  # increase; with delta_HR_eff = 0.15 the drop should be material.
  expect_lt(s015$mean_pp, s0$mean_pp)
})
