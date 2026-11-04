# Helpers for test-chronology-amd8.R (Amendment 8, item 1; evolveTrial 0.0.0.9002).
#
# amd8_args(): wrapper-style single-arm configuration, mirroring
# context-bo-trial-paper/R/evolve_sim.R (EVOLVE cohort B, futility-only
# interims, 24-month PWE horizon, 120-month follow-up, 72-month cutoff).
#
# amd8_chronology(): runs one trial per seed (num_simulations = 1,
# set.seed(seed) before each call) and records each trial's chronology by
# tracing two internal functions that exist unchanged in 0.0.0.9001 and
# 0.0.0.9002, so the same instrumentation runs on both engines:
#   * interim_check(): the calendar time of every interim look entered with
#     at least one arm recruiting (looks with all arms stopped are no-ops);
#   * slice_arm_data_at_time(): the last call in a trial is the final
#     analysis slice; its registry gives every enrolment time.
# N (Exp_N) and the final-analysis or stop time (Exp_Time) are read from the
# returned results data frame, which for one replicate are that trial's values.

amd8_theta_t95 <- list(total_n = 117, eff = 0.895954670800146, fut = 0.399967861412818,
                       ev = 20, beat = 7.09081887219977, margin = 0.114915937521501,
                       prior_strength = 0.265756152558728, n_intervals = 6)

amd8_args <- function(theta = amd8_theta_t95, lambda = 0.92, kappa = 1.3, m0 = 6,
                      true_median = 9, fu_final = 6, max_trial_time = 72,
                      num_posterior_draws = 1000L, num_simulations = 1L) {
  n_int <- as.integer(theta$n_intervals)
  list(
    num_simulations = num_simulations, diagnostics = FALSE,
    parallel_replicates = FALSE, progress = FALSE,
    arm_names = "Experimental", reference_arm_name = NULL, compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Experimental = kappa),
    weibull_median_true_arms = c(Experimental = true_median),
    null_median_arms = c(Experimental = m0),
    futility_median_arms = c(Experimental = m0 + theta$margin),
    interval_cutpoints_sim = seq(0, 24, length.out = n_int + 1L),
    max_follow_up_sim = 120, censor_max_time_sim = 120,
    max_trial_time = max_trial_time,
    prior_alpha_params_model = rep(theta$prior_strength, n_int),
    prior_beta_params_model = rep(theta$prior_strength, n_int),
    num_posterior_draws = num_posterior_draws,
    cohort_size_per_arm = 1L,
    max_total_patients_per_arm = c(Experimental = as.integer(theta$total_n)),
    min_patients_per_arm = 10L, min_patients_for_analysis = 0L,
    overall_accrual_rate = lambda, randomization_probs = c(Experimental = 1.0),
    efficacy_stopping_rule_hc = FALSE,
    efficacy_threshold_hc_prob = c(Experimental = theta$eff),
    futility_threshold_hc_prob = c(Experimental = theta$fut),
    futility_stopping_rule_hc = TRUE,
    enable_efficacy_stopping = TRUE, enable_futility_stopping = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE, futility_stopping_rule_vs_ref = FALSE,
    efficacy_threshold = 0.98, futility_threshold = 0.60,
    efficacy_threshold_vs_ref_prob = 0.98, futility_threshold_vs_ref_prob = 0.60,
    use_ph_model_vs_ref = FALSE, ph_loghr_prior_mean = 0, ph_loghr_prior_sd = 1,
    min_events_hc = as.integer(theta$ev), min_median_followup_hc = 0,
    interim_calendar_beat = theta$beat,
    min_events_per_arm = as.integer(theta$ev), min_median_followup_per_arm = 3,
    min_person_time_frac_per_arm = 0,
    min_follow_up_at_final = fu_final,
    rebalance_after_events = 50L,
    pred_success_pp_threshold_hc = 0.90, pred_futility_pp_threshold_hc = 0.20,
    num_posterior_draws_pred = 1000L, predictive_fast = FALSE,
    median_pfs_success_threshold_arms = c(Experimental = m0),
    median_pfs_futility_threshold_arms = c(Experimental = m0 + theta$margin),
    final_success_posterior_prob_threshold = theta$eff,
    final_futility_posterior_prob_threshold = theta$fut
  )
}

amd8_filter <- function(args) {
  # as run_scenarios() does: keep only formals of run_simulation_pure
  args[intersect(names(args), names(formals(evolveTrial::run_simulation_pure)))]
}

amd8_chronology <- function(args, seeds) {
  # P3-1: the traces below see only the R body; since P3-1 both C++ engines admit these frames, so pin the R body here
  old_opt <- options(evolveTrial.use_cpp_engine = FALSE)
  on.exit(options(old_opt), add = TRUE)
  ns <- asNamespace("evolveTrial")
  args <- amd8_filter(args)
  rec <- new.env()
  suppressMessages({
    trace("interim_check", where = ns, print = FALSE,
          tracer = bquote(if (any(state$arm_status == "recruiting"))
            assign("interims", c(get("interims", envir = .(rec)), current_time),
                   envir = .(rec))))
    trace("slice_arm_data_at_time", where = ns, print = FALSE,
          tracer = bquote({
            assign("slice_time", calendar_time, envir = .(rec))
            assign("enroll", registry_df$enroll_time[!is.na(registry_df$id)], envir = .(rec))
          }))
  })
  on.exit(suppressMessages({
    untrace("interim_check", where = ns)
    untrace("slice_arm_data_at_time", where = ns)
  }), add = TRUE)

  rows <- vector("list", length(seeds))
  interim_list <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    rec$interims <- numeric(0); rec$slice_time <- NA_real_; rec$enroll <- numeric(0)
    set.seed(seeds[i])
    res <- do.call(evolveTrial::run_simulation_pure, args)
    stopped <- (res$PET_Futility + res$PET_Efficacy) > 0
    en <- rec$enroll
    interim_list[[i]] <- rec$interims
    rows[[i]] <- data.frame(
      seed = seeds[i],
      N = res$Exp_N, final_time = res$Exp_Time, stopped = stopped,
      final_eff = res$Pr_Final_Efficacy, final_fut = res$Pr_Final_Futility,
      slice_time = rec$slice_time,
      n_enrolled = length(en),
      last_enroll = if (length(en)) max(en) else 0,
      n_enroll_ge_cut = sum(en >= args$max_trial_time),
      n_enroll_before_final = sum(en < res$Exp_Time),
      n_interims = length(rec$interims),
      max_interim = if (length(rec$interims)) max(rec$interims) else -Inf,
      n_interims_before_final = sum(rec$interims < res$Exp_Time),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  attr(out, "interims") <- interim_list
  out
}

amd8_theta_t50 <- list(total_n = 65, eff = 0.894653709079939, fut = 0.380812975336976,
                       ev = 14, beat = 2.09782153331422, margin = 0.475391735149199,
                       prior_strength = 0.451397409221281, n_intervals = 4)
