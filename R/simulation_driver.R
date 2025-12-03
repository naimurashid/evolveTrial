
#' Run a full set of evolveTrial simulations for a design specification
#'
#' `run_simulation_pure()` is the workhorse simulator for evolveTrial.  It
#' enrols patients according to the supplied accrual plan, applies interim
#' gating/decision rules, optionally rebalances interval cut points, and
#' carries each replicate forward to the final analysis.  Operating
#' characteristics are returned summarised per arm.
#'
#' @param num_simulations Number of Monte Carlo replicates to run.
#' @param arm_names Character vector naming the trial arms.
#' @param reference_arm_name Character scalar naming the control/reference arm.
#' @param compare_arms_option Logical; `TRUE` evaluates vs-reference logic,
#'   `FALSE` evaluates arms independently against historical control targets.
#' @param weibull_shape_true_arms Named numeric vector of Weibull shape
#'   parameters under the truth.
#' @param weibull_median_true_arms Named numeric vector of true median PFS
#'   (months) for each arm.
#' @param null_median_arms Named numeric vector of null (historical control)
#'   medians used for single-arm evaluations.
#' @param futility_median_arms Named numeric vector of futility medians for
#'   single-arm logic.
#' @param interval_cutpoints_sim Numeric vector of interval boundaries (months)
#'   for piecewise exponential modelling.
#' @param max_follow_up_sim Maximum administrative follow-up time (months).
#' @param censor_max_time_sim Upper bound for random censoring draws (months).
#' @param prior_alpha_params_model Numeric vector of Gamma prior shape
#'   parameters for the piecewise exponential hazards.
#' @param prior_beta_params_model Numeric vector of Gamma prior rate parameters
#'   for the piecewise exponential hazards.
#' @param num_posterior_draws Number of posterior draws used for final analyses.
#' @param num_posterior_draws_interim Optional integer overriding the number of
#'   posterior draws used at interim looks.
#' @param cohort_size_per_arm Size of each enrolment batch per arm (typically 1).
#' @param max_total_patients_per_arm Named integer vector of per-arm sample size
#'   caps.
#' @param min_patients_for_analysis Minimum number of patients required to
#'   evaluate an arm in the single-arm path. If not specified, defaults to 0,
#'   allowing interim analyses to occur even with very few patients. Set this
#'   to a higher value to prevent interim analyses until a certain number of
#'   patients have been enrolled.
#' @param efficacy_stopping_rule_hc Logical; enable interim efficacy checks for
#'   the historical-control path.
#' @param efficacy_threshold_current_prob_hc **DEPRECATED**. Use
#'   \code{efficacy_threshold_hc_prob} instead. Interim success probability
#'   threshold for single-arm logic.
#' @param posterior_futility_threshold_hc **DEPRECATED**. Use
#'   \code{futility_threshold_hc_prob} instead. Interim futility probability
#'   threshold for single-arm logic.
#' @param efficacy_threshold_hc_prob Interim success probability threshold for
#'   single-arm logic (preferred harmonized name).
#' @param futility_threshold_hc_prob Interim futility probability threshold for
#'   single-arm logic (preferred harmonized name).
#' @param futility_stopping_rule_hc Logical; enable interim futility checks for
#'   the single-arm path.
#' @param efficacy_stopping_rule_vs_ref Logical; enable interim efficacy checks
#'   for the vs-reference path.
#' @param futility_stopping_rule_vs_ref Logical; enable interim futility checks
#'   for the vs-reference path.
#' @param efficacy_threshold_vs_ref_prob Posterior superiority threshold for
#'   vs-reference decisions.
#' @param futility_threshold_vs_ref_prob Posterior inferiority threshold for
#'   vs-reference decisions.
#' @param compare_arms_futility_margin Absolute median difference used when
#'   defining vs-reference futility.
#' @param compare_arms_hr_margin Optional hazard-ratio margin used when
#'   `use_ph_model_vs_ref = TRUE`.
#' @param use_ph_model_vs_ref Logical; use the proportional-hazards joint model
#'   for vs-reference comparisons.
#' @param ph_loghr_prior_mean Mean of the normal prior on the log hazard ratio
#'   (PH model).
#' @param ph_loghr_prior_sd Standard deviation of the normal prior on the log
#'   hazard ratio.
#' @param median_pfs_success_threshold_arms Named numeric vector of median PFS
#'   thresholds for declaring final success per arm.
#' @param final_success_posterior_prob_threshold Posterior probability threshold
#'   for final success declarations.
#' @param median_pfs_futility_threshold_arms Named numeric vector of median PFS
#'   futility thresholds for final analyses.
#' @param final_futility_posterior_prob_threshold Posterior probability
#'   threshold for final futility declarations.
#' @param overall_accrual_rate Expected accrual rate (patients per month).
#' @param randomization_probs Named numeric vector of randomisation probabilities.
#' @param min_follow_up_at_final Additional follow-up (months) required after
#'   last enrolment before the final analysis.
#' @param min_events_for_analysis **DEPRECATED**. Use \code{min_events_hc}
#'   instead. Minimum events required for interim review (global gate).
#' @param min_median_followup **DEPRECATED**. Use \code{min_median_followup_hc}
#'   instead. Minimum median follow-up required for interim review (global gate).
#' @param min_events_hc Minimum events required for single-arm interim review
#'   (preferred harmonized name).
#' @param min_median_followup_hc Minimum median follow-up required for single-arm
#'   interim review (preferred harmonized name).
#' @param interim_calendar_beat Calendar spacing (months) between scheduled
#'   interim looks when person-time milestones are not used.
#' @param diagnostics Logical; if `TRUE` prints interim diagnostic messages.
#' @param pred_success_pp_threshold_hc Predictive probability threshold for
#'   interim success in the single-arm predictive look (if enabled).
#' @param pred_futility_pp_threshold_hc Predictive probability threshold for
#'   interim futility in the single-arm predictive look (if enabled).
#' @param num_posterior_draws_pred Number of posterior draws used inside
#'   predictive probability calculations.
#' @param predictive_fast Logical; switch to the analytic predictive
#'   approximations.
#' @param min_events_per_arm Optional per-arm minimum event gate for vs-reference.
#' @param min_median_followup_per_arm Optional per-arm minimum median follow-up
#'   gate for vs-reference.
#' @param min_person_time_frac_per_arm Optional per-arm proportion of planned
#'   person-time required before evaluating vs-reference decisions.
#' @param person_time_milestones Optional numeric vector (fractions of total
#'   planned person-time) triggering interim looks.
#' @param latest_calendar_look Backstop calendar time for person-time schedules.
#' @param rebalance_after_events Optional integer; when non-`NULL` the piecewise
#'   cut points are re-estimated once that number of events has accrued.
#' @param parallel_replicates Logical; if `TRUE`, distribute Monte Carlo
#'   replicates across a parallel cluster.
#' @param num_workers Optional integer specifying the number of workers when
#'   `parallel_replicates = TRUE`. Defaults to `parallel::detectCores() - 1`.
#' @param cluster_type Type of parallel cluster to spawn when distributing
#'   replicates. One of `"auto"` (default, uses FORK on Unix, PSOCK on Windows),
#'   `"PSOCK"`, or `"FORK"`. FORK clusters are faster on Linux/macOS as they
#'   share memory and don't require package loading in workers.
#' @param cluster Optional pre-existing parallel cluster to reuse. When provided,
#'   the cluster is used for parallel execution and NOT stopped on exit. This
#'   enables cluster pooling for repeated calls (e.g., in Bayesian optimization).
#'   Create with `evolveTrial::create_simulation_cluster()` or `parallel::makeCluster()`.
#' @param progress Logical; show the simulation progress bar when running
#'   sequentially. Automatically disabled for parallel replicate execution.
#'
#' @return A data frame with one row per arm and columns summarising operating
#'   characteristics such as type I error / power, PETs, final decision
#'   probabilities, and expected sample size.
#'
#' @details When \code{parallel_replicates = TRUE}, results will vary based on
#'   \code{num_workers} due to different random number stream partitioning.
#'   For exact reproducibility across runs, use \code{parallel_replicates = FALSE}.
#'
#'   During package development with \code{devtools::load_all()}, parallel workers
#'   will load the installed package version, not the development code. For
#'   testing development changes, either use \code{parallel_replicates = FALSE} or
#'   reinstall the package with \code{devtools::install()}.
#' @export
run_simulation_pure <- function(
    num_simulations,
    arm_names,
    reference_arm_name,
    compare_arms_option,
    weibull_shape_true_arms,
    weibull_median_true_arms,
    null_median_arms = NULL,
    futility_median_arms = NULL,
    interval_cutpoints_sim,
    max_follow_up_sim,
    censor_max_time_sim,
    prior_alpha_params_model,
    prior_beta_params_model,
    num_posterior_draws,
    num_posterior_draws_interim = NULL,
    cohort_size_per_arm,
    max_total_patients_per_arm,
    min_patients_for_analysis = NULL,
    # stopping rule params
    efficacy_stopping_rule_hc = FALSE,
    efficacy_threshold_current_prob_hc = NULL,
    posterior_futility_threshold_hc = NULL,
    # NEW harmonized parameter names (preferred)
    efficacy_threshold_hc_prob = NULL,
    futility_threshold_hc_prob = NULL,
    futility_stopping_rule_hc = FALSE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE,
    efficacy_threshold_vs_ref_prob = NULL,
    futility_threshold_vs_ref_prob = NULL,
    compare_arms_futility_margin = 0,
    compare_arms_hr_margin = NULL,
    use_ph_model_vs_ref = FALSE,
    ph_loghr_prior_mean = 0,
    ph_loghr_prior_sd = 1,
    # final analysis params
    median_pfs_success_threshold_arms = NULL,
    final_success_posterior_prob_threshold = 0.85,
    median_pfs_futility_threshold_arms = NULL,
    final_futility_posterior_prob_threshold = 0.85,
    # accrual & randomization
    overall_accrual_rate,
    randomization_probs,
    min_follow_up_at_final = 0,
    # legacy info gates + calendar-beat interims
    min_events_for_analysis = NULL,
    min_median_followup = NULL,
    # NEW harmonized parameter names (preferred)
    min_events_hc = NULL,
    min_median_followup_hc = NULL,
    interim_calendar_beat = 2,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100,
    predictive_fast = FALSE,
    # NEW: per-arm info gates
    min_events_per_arm = NULL,
    min_median_followup_per_arm = NULL,
    min_person_time_frac_per_arm = 0,
    # NEW: person-time milestones
    person_time_milestones = NULL,
    latest_calendar_look = Inf,
    # optional interval rebalancing
    rebalance_after_events = NULL,
    # replicate-level parallelisation
    parallel_replicates = FALSE,
    num_workers = NULL,
    cluster_type = c("auto", "PSOCK", "FORK"),
    cluster = NULL,
    progress = interactive()
) {
  cluster_type <- match.arg(cluster_type)

  # Auto-detect optimal cluster type: FORK on Unix (faster), PSOCK on Windows
  if (cluster_type == "auto") {
    cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  }

  # ---- PARAMETER DEPRECATION HANDLING ----------------------------------------
  # Map old parameter names to new harmonized names with warnings






  # Default to 0 if NULL
  if (is.null(min_events_hc)) min_events_hc <- 0


  # Default to 0 if NULL
  if (is.null(min_median_followup_hc)) min_median_followup_hc <- 0

  # min_patients_for_analysis default
  if (is.null(min_patients_for_analysis)) min_patients_for_analysis <- 0

  # ---------------------------------------------------------------------------

  weibull_scale_true_arms <- sapply(arm_names, function(arm) {
    calculate_weibull_scale(weibull_median_true_arms[arm], weibull_shape_true_arms[arm])
  })
  names(weibull_scale_true_arms) <- arm_names
  
  results_data <- data.frame(
    Arm_Name = arm_names,
    True_Median = round(weibull_median_true_arms, 2),
    Type_I_Error_or_Power = 0,
    PET_Efficacy = 0,
    PET_Futility = 0,
    Pr_Reach_Max_N = 0,
    Pr_Final_Efficacy = 0,
    Pr_Final_Futility = 0,
    Pr_Final_Inconclusive = 0,
    Exp_N = 0,
    Exp_Events = 0,
    stringsAsFactors = FALSE
  )
  
  show_progress <- isTRUE(progress) && !isTRUE(parallel_replicates)
  if (show_progress) {
    pb <- progress::progress_bar$new(
      format = "  Sims [:bar] :percent in :elapsed",
      total = num_simulations, clear = FALSE, width = 60
    )
  }
  
  # total max person-time across arms (for milestone schedule)
  max_PT_per_arm <- setNames(as.numeric(max_total_patients_per_arm) * max_follow_up_sim,
                             names(max_total_patients_per_arm))
  total_max_PT <- sum(max_PT_per_arm)

  # Build absolute PT milestones (if any)
  pt_targets_abs <- NULL
  if (!is.null(person_time_milestones)) {
    stopifnot(all(person_time_milestones > 0 & person_time_milestones <= 1))
    pt_targets_abs <- sort(unique(person_time_milestones)) * total_max_PT
  }

  num_draws_interim <- if (is.null(num_posterior_draws_interim)) {
    num_posterior_draws
  } else {
    num_posterior_draws_interim
  }

  # pack args for interim_check
  args <- list(
    arm_names = arm_names,
    reference_arm_name = reference_arm_name,
    compare_arms_option = compare_arms_option,
    interval_cutpoints_sim = interval_cutpoints_sim,
    max_follow_up_sim = max_follow_up_sim,
    randomization_probs = randomization_probs,
    prior_alpha_params_model = prior_alpha_params_model,
    prior_beta_params_model  = prior_beta_params_model,
    num_posterior_draws = num_posterior_draws,
    num_posterior_draws_interim = num_draws_interim,
    min_patients_for_analysis = min_patients_for_analysis,
    min_events_hc = min_events_hc,  # Use new harmonized name
    min_median_followup_hc = min_median_followup_hc,  # Use new harmonized name
    null_median_arms = null_median_arms,
    futility_median_arms = futility_median_arms,
    efficacy_threshold_hc_prob = efficacy_threshold_hc_prob,  # Use new harmonized name
    futility_threshold_hc_prob = futility_threshold_hc_prob,  # Use new harmonized name
    efficacy_threshold_vs_ref_prob = efficacy_threshold_vs_ref_prob,
    futility_threshold_vs_ref_prob = futility_threshold_vs_ref_prob,
    compare_arms_futility_margin = compare_arms_futility_margin,
    compare_arms_hr_margin = compare_arms_hr_margin,
    use_ph_model_vs_ref = use_ph_model_vs_ref,
    ph_loghr_prior_mean = ph_loghr_prior_mean,
    ph_loghr_prior_sd = ph_loghr_prior_sd,
    max_total_patients_per_arm = max_total_patients_per_arm,
    median_pfs_success_threshold_arms = median_pfs_success_threshold_arms,
    final_success_posterior_prob_threshold = final_success_posterior_prob_threshold,
    final_futility_posterior_prob_threshold = final_futility_posterior_prob_threshold,
    censor_max_time_sim = censor_max_time_sim,
    predictive_fast = predictive_fast,
    num_posterior_draws_pred = num_posterior_draws_pred,
    predictive_boundary_band  = 0.12,
    run_pred_every_k_beats    = 2,
    pred_success_pp_threshold_hc  = pred_success_pp_threshold_hc,
    pred_futility_pp_threshold_hc = pred_futility_pp_threshold_hc,
    # per-arm gates + PT metadata
    min_events_per_arm = min_events_per_arm,
    min_median_followup_per_arm = min_median_followup_per_arm,
    min_person_time_frac_per_arm = min_person_time_frac_per_arm,
    max_PT_per_arm = max_PT_per_arm,
    rebalance_after_events = rebalance_after_events
  )
  
  num_intervals <- length(interval_cutpoints_sim) - 1
  interval_lengths_base <- diff(interval_cutpoints_sim)

  sum_final_n       <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_events  <- setNames(numeric(length(arm_names)), arm_names)
  sum_stop_efficacy <- setNames(numeric(length(arm_names)), arm_names)
  sum_stop_futility <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_efficacy <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_futility <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_inconclusive <- setNames(numeric(length(arm_names)), arm_names)

  tick_fun <- if (show_progress) function() pb$tick() else function() invisible(NULL)

  base_args_for_interim <- args

  simulate_chunk <- function(sim_indices, seed = NULL, tick = function() {}) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    chunk_sum_final_n <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_final_events <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_stop_eff <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_stop_fut <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_final_eff <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_final_fut <- setNames(numeric(length(arm_names)), arm_names)
    chunk_sum_final_inc <- setNames(numeric(length(arm_names)), arm_names)

    # Lightweight helper for rebalancing: count events without interval metrics
    extract_events_fast <- function(registry_df, calendar_time, max_follow_up) {
      reg <- get_active_registry(registry_df)
      if (nrow(reg) == 0) return(list(count = 0L, times = numeric(0)))

      time_since_enroll <- pmax(0, calendar_time - reg$enroll_time)
      time_available    <- pmin(time_since_enroll, max_follow_up)
      active <- time_available > 0
      if (!any(active)) return(list(count = 0L, times = numeric(0)))

      time_available <- time_available[active]
      te <- reg$true_event_time[active]
      rc <- reg$random_censor_time[active]
      observed_time <- pmin(te, rc, time_available)
      events_mask <- te <= pmin(rc, time_available)

      list(
        count = sum(events_mask),
        times = if (any(events_mask)) observed_time[events_mask] else numeric(0)
      )
    }

    for (sim_idx in sim_indices) {
      args_local <- base_args_for_interim
      # Reset cutpoints to baseline for each simulation (handles rebalancing correctly)
      interval_cutpoints_current <- interval_cutpoints_sim
      args_local$interval_cutpoints_sim <- interval_cutpoints_current
      interval_lengths <- interval_lengths_base
      rebalance_threshold <- rebalance_after_events
      rebalance_done <- is.null(rebalance_threshold)
      last_enroll_time <- 0

      state <- make_state(arm_names, max_total_patients_per_arm)
      current_time <- 0.0
      next_calendar_look <- interim_calendar_beat
      next_pt_idx <- 1L
      patient_id <- 0L

      is_eligible <- function(st) {
        which((st$arm_status == "recruiting") & (st$enrolled_counts < max_total_patients_per_arm))
      }

      while (length(is_eligible(state)) > 0) {
        interarrival <- rexp(1, rate = overall_accrual_rate)
        current_time <- current_time + interarrival

        if (!rebalance_done && !is.null(rebalance_threshold)) {
          total_events <- 0L
          event_times_all <- numeric(0)
          for (arm in arm_names) {
            ev_res <- extract_events_fast(
              state$registries[[arm]], current_time, max_follow_up_sim
            )
            total_events <- total_events + ev_res$count
            if (ev_res$count > 0) {
              event_times_all <- c(event_times_all, ev_res$times)
            }
          }
          if (total_events >= rebalance_threshold) {
            new_cuts <- rebalance_cutpoints_from_events(
              event_times_all, max_follow_up_sim, num_intervals
            )
            if (!is.null(new_cuts)) {
              interval_cutpoints_current <- new_cuts
              args_local$interval_cutpoints_sim <- new_cuts
              interval_lengths <- diff(new_cuts)
              rebalance_done <- TRUE
              if (diagnostics) {
                message(sprintf(
                  "Rebalanced interval cutpoints at t=%.2f using %d observed events",
                  current_time, total_events
                ))
              }
            } else {
              rebalance_done <- TRUE
            }
          }
        }

        if (!is.null(pt_targets_abs)) {
          total_PT_now <- cum_person_time_all_arms(state, current_time, max_follow_up_sim,
                                                   interval_cutpoints_current, arm_names)
          if (is.null(state$.__backstop_fired__)) state$.__backstop_fired__ <- FALSE

          while (!is.null(pt_targets_abs) &&
                 next_pt_idx <= length(pt_targets_abs) &&
                 (total_PT_now >= pt_targets_abs[next_pt_idx] ||
                  (!state$.__backstop_fired__ && is.finite(latest_calendar_look) && current_time >= latest_calendar_look))) {

            state <- interim_check(state, current_time, args_local, diagnostics = diagnostics)

            if (total_PT_now >= pt_targets_abs[next_pt_idx]) {
              next_pt_idx <- next_pt_idx + 1L
            }
            if (!state$.__backstop_fired__ && is.finite(latest_calendar_look) && current_time >= latest_calendar_look) {
              state$.__backstop_fired__ <- TRUE
            }

            if (all(state$arm_status != "recruiting")) break
            total_PT_now <- cum_person_time_all_arms(state, current_time, max_follow_up_sim,
                                                     interval_cutpoints_current, arm_names)
          }
        } else {
          if (current_time >= next_calendar_look) {
            state <- interim_check(state, current_time, args_local, diagnostics = diagnostics)
            next_calendar_look <- next_calendar_look + interim_calendar_beat
          }
        }

        if (all(state$arm_status != "recruiting")) break

        elig_idx <- is_eligible(state)
        if (length(elig_idx) == 0) break

        if (compare_arms_option) {
          exp_idx <- which(arm_names != reference_arm_name)
          exp_active <- exp_idx[state$arm_status[arm_names[exp_idx]] == "recruiting"]
          exp_active <- exp_active[state$enrolled_counts[arm_names[exp_active]] <
                                     max_total_patients_per_arm[arm_names[exp_active]]]
          if (length(exp_active) == 0) break
        }

        elig_arms <- arm_names[elig_idx]

        probs <- randomization_probs[elig_arms]
        probs <- probs / sum(probs)
        chosen_arm <- sample(elig_arms, size = 1, prob = probs)

        patient_id <- patient_id + 1L
        t_event_true <- rweibull(1, shape = weibull_shape_true_arms[chosen_arm],
                                 scale = weibull_scale_true_arms[chosen_arm])
        t_random_censor <- runif(1, min = 0, max = censor_max_time_sim)

        # PERFORMANCE: Use indexed assignment instead of rbind (O(1) vs O(n) per patient)
        idx <- state$registry_row_idx[[chosen_arm]]
        state$registries[[chosen_arm]][idx, "id"] <- patient_id
        state$registries[[chosen_arm]][idx, "enroll_time"] <- current_time
        state$registries[[chosen_arm]][idx, "true_event_time"] <- t_event_true
        state$registries[[chosen_arm]][idx, "random_censor_time"] <- t_random_censor
        state$registry_row_idx[[chosen_arm]] <- idx + 1L
        state$enrolled_counts[chosen_arm] <- state$enrolled_counts[chosen_arm] + 1L
        last_enroll_time <- current_time
      }

      state <- interim_check(state, current_time, args_local, diagnostics = diagnostics)

      final_time <- last_enroll_time + min_follow_up_at_final

      if (!identical(interval_cutpoints_current, interval_cutpoints_sim)) {
        interval_lengths <- diff(interval_cutpoints_current)
      }

      ref_slice_final <- NULL
      if (compare_arms_option) {
        ref_slice_final <- slice_arm_data_at_time(
          state$registries[[reference_arm_name]], final_time,
          max_follow_up_sim, interval_cutpoints_current
        )
      }

      final_eff_vec <- setNames(integer(length(arm_names)), arm_names)
      final_fut_vec <- setNames(integer(length(arm_names)), arm_names)
      final_inc_vec <- setNames(integer(length(arm_names)), arm_names)
      # PERFORMANCE: Cache arm slices to avoid redundant slice_arm_data_at_time calls
      cached_arm_slices <- list()

      for (arm in arm_names) {
        if (state$arm_status[arm] != "recruiting") next

        arm_slice <- slice_arm_data_at_time(state$registries[[arm]], final_time,
                                            max_follow_up_sim, interval_cutpoints_current)
        # Cache for reuse when computing final events
        cached_arm_slices[[arm]] <- arm_slice

        if (!compare_arms_option) {
          post_arm <- draw_posterior_hazard_samples(
            num_intervals = num_intervals,
            events_per_interval = arm_slice$metrics$events_per_interval,
            person_time_per_interval = arm_slice$metrics$person_time_per_interval,
            prior_alpha_params = prior_alpha_params_model,
            prior_beta_params  = prior_beta_params_model,
            num_samples = num_posterior_draws
          )
          # PERFORMANCE: Use C++ matrix version instead of apply() for 20-30x speedup
          med_arm <- calculate_median_survival_matrix_cpp(post_arm, interval_lengths)
          # Final success check
          success_thr <- if (!is.null(median_pfs_success_threshold_arms) && arm %in% names(median_pfs_success_threshold_arms)) {
            median_pfs_success_threshold_arms[arm]
          } else {
            Inf  # Impossible to cross if not specified
          }
          p_eff <- if (is.finite(success_thr)) mean(med_arm > success_thr) else 0

          if (!is.na(p_eff) && p_eff >= final_success_posterior_prob_threshold) {
            final_eff_vec[arm] <- 1L
          } else {
            # Final futility check
            futility_thr <- if (!is.null(median_pfs_futility_threshold_arms) && arm %in% names(median_pfs_futility_threshold_arms)) {
              median_pfs_futility_threshold_arms[arm]
            } else {
              -Inf  # Impossible to cross if not specified
            }
            p_fut <- if (is.finite(futility_thr)) mean(med_arm < futility_thr) else 0

            if (!is.na(p_fut) && p_fut >= final_futility_posterior_prob_threshold) {
              final_fut_vec[arm] <- 1L
            } else {
              final_inc_vec[arm] <- 1L
            }
          }
        } else {
          if (arm == reference_arm_name) {
            final_inc_vec[arm] <- 1L
          } else {
            med_samples <- sample_vs_ref_medians(
              slCtrl = ref_slice_final,
              slTrt = arm_slice,
              args = args_local,
              num_samples = num_posterior_draws
            )
            if (isTRUE(args_local$use_ph_model_vs_ref) && !is.null(med_samples$logHR)) {
              log_margin <- 0
              if (!is.null(args_local$compare_arms_hr_margin)) {
                log_margin <- log1p(max(0, args_local$compare_arms_hr_margin))
              }
              p_eff_ref <- mean(med_samples$logHR < 0)
              p_fut_ref <- mean(med_samples$logHR >= log_margin)
            } else {
              margin_abs <- coalesce_num(compare_arms_futility_margin, 0)
              pr <- final_vsref_probs_abs(med_samples$medTrt, med_samples$medCtrl, margin_abs)
              p_eff_ref <- pr$p_eff_ref
              p_fut_ref <- pr$p_fut_ref
            }

            if (p_eff_ref >= efficacy_threshold_vs_ref_prob) {
              final_eff_vec[arm] <- 1L
            } else if (p_fut_ref >= futility_threshold_vs_ref_prob) {
              final_fut_vec[arm] <- 1L
            } else {
              final_inc_vec[arm] <- 1L
            }
          }
        }
      }

      final_n_vec <- setNames(numeric(length(arm_names)), arm_names)
      final_events_vec <- setNames(numeric(length(arm_names)), arm_names)
      for (arm in arm_names) {
        final_n_vec[arm] <- if (is.na(state$sim_final_n_current_run[arm])) {
          state$enrolled_counts[arm]
        } else {
          state$sim_final_n_current_run[arm]
        }
        # PERFORMANCE: Use cached slice if available, otherwise compute
        if (arm %in% names(cached_arm_slices)) {
          arm_slice_for_events <- cached_arm_slices[[arm]]
        } else {
          arm_slice_for_events <- slice_arm_data_at_time(
            state$registries[[arm]], final_time,
            max_follow_up_sim, interval_cutpoints_current
          )
        }
        final_events_vec[arm] <- arm_slice_for_events$metrics$events_total
      }

      chunk_sum_final_n <- chunk_sum_final_n + final_n_vec
      chunk_sum_final_events <- chunk_sum_final_events + final_events_vec
      chunk_sum_stop_eff <- chunk_sum_stop_eff + state$stop_efficacy_per_sim_row
      chunk_sum_stop_fut <- chunk_sum_stop_fut + state$stop_futility_per_sim_row
      chunk_sum_final_eff <- chunk_sum_final_eff + final_eff_vec
      chunk_sum_final_fut <- chunk_sum_final_fut + final_fut_vec
      chunk_sum_final_inc <- chunk_sum_final_inc + final_inc_vec

      tick()
    }

    list(
      sum_final_n = chunk_sum_final_n,
      sum_final_events = chunk_sum_final_events,
      sum_stop_eff = chunk_sum_stop_eff,
      sum_stop_fut = chunk_sum_stop_fut,
      sum_final_eff = chunk_sum_final_eff,
      sum_final_fut = chunk_sum_final_fut,
      sum_final_inc = chunk_sum_final_inc,
      n_sims = length(sim_indices)
    )
  }

  aggregate_results <- function(partials) {
    for (res in partials) {
      sum_final_n       <<- sum_final_n + res$sum_final_n
      # Defensive: handle workers running old code without sum_final_events
      if (!is.null(res$sum_final_events)) {
        sum_final_events  <<- sum_final_events + res$sum_final_events
      }
      sum_stop_efficacy <<- sum_stop_efficacy + res$sum_stop_eff
      sum_stop_futility <<- sum_stop_futility + res$sum_stop_fut
      sum_final_efficacy <<- sum_final_efficacy + res$sum_final_eff
      sum_final_futility <<- sum_final_futility + res$sum_final_fut
      sum_final_inconclusive <<- sum_final_inconclusive + res$sum_final_inc
    }
  }

  chunk_results <- list()

  if (isTRUE(parallel_replicates) && num_simulations > 1) {
    # Determine number of workers
    # If external cluster is provided, use its size; otherwise use num_workers/detectCores
    workers <- if (!is.null(cluster)) {
      # Use the size of the supplied cluster
      length(cluster)
    } else if (is.null(num_workers)) {
      max(1L, parallel::detectCores() - 1L)
    } else {
      as.integer(num_workers)
    }
    workers <- max(1L, min(workers, num_simulations))

    # Sequential threshold: parallel overhead exceeds benefit for small rep counts

    # Based on profiling: cluster spawn/teardown ~5-10 sec, simulation ~1-2 ms/rep
    # Threshold of 100 ensures parallel overhead is worthwhile
    seq_threshold <- getOption("evolveTrial.sequential_threshold", 100L)

    # Check if parallelization is worthwhile
    if (num_simulations < seq_threshold) {
      if (interactive()) {
        message("Few simulations (", num_simulations,
                ") below threshold (", seq_threshold,
                "); running sequentially for efficiency")
      }
      chunk_results <- list(simulate_chunk(seq_len(num_simulations), seed = NULL, tick = tick_fun))
    } else {
      chunks <- parallel::splitIndices(num_simulations, workers)
      chunks <- chunks[lengths(chunks) > 0]

      pkg_name <- utils::packageName()
      if (is.null(pkg_name) || identical(pkg_name, "")) {
        pkg_name <- "evolveTrial"
      }

      # Cluster pooling: reuse external cluster if provided
      using_external_cluster <- !is.null(cluster)

      if (using_external_cluster) {
        cl <- cluster
        # Do not stop external cluster on exit - caller is responsible

        # Ensure package is loaded on external PSOCK clusters
        # (FORK clusters inherit parent environment, so no loading needed)
        # Check cluster type by examining the first node's class
        cluster_class <- class(cl[[1]])[1]
        is_psock <- grepl("SOCK", cluster_class, ignore.case = TRUE)

        if (is_psock) {
          # Check if package is already loaded on workers
          pkg_loaded <- tryCatch({
            test_result <- parallel::clusterCall(cl, function(pkg) {
              isNamespaceLoaded(pkg)
            }, pkg_name)
            all(unlist(test_result))
          }, error = function(e) FALSE)

          if (!pkg_loaded) {
            parallel::clusterCall(cl, function(pkg) {
              suppressPackageStartupMessages(require(pkg, character.only = TRUE))
              NULL
            }, pkg_name)
          }
        }
      } else {
        # Diagnostic for development workflow
        if (interactive()) {
          pkg_version <- tryCatch(
            as.character(utils::packageVersion(pkg_name)),
            error = function(e) "not installed"
          )
          message("Parallel execution: ", workers, " workers (", cluster_type,
                  ") will load ", pkg_name, " version ", pkg_version)
        }

        cl <- parallel::makeCluster(workers, type = cluster_type)
        on.exit(parallel::stopCluster(cl), add = TRUE)

        # Only load package for PSOCK clusters (FORK inherits parent environment)
        if (cluster_type == "PSOCK") {
          parallel::clusterCall(cl, function(pkg) {
            suppressPackageStartupMessages(require(pkg, character.only = TRUE))
            NULL
          }, pkg_name)
        }
      }

      seed_list <- replicate(length(chunks), sample.int(.Machine$integer.max, 1L), simplify = TRUE)

      chunk_results <- tryCatch(
        parallel::parLapply(
          cl,
          seq_along(chunks),
          function(idx, chunk_indices, seed_vals) {
            simulate_chunk(chunk_indices[[idx]], seed = seed_vals[[idx]], tick = function() {})
          },
          chunk_indices = chunks,
          seed_vals = as.list(seed_list)
        ),
        error = function(e) {
          stop("Parallel simulation failed: ", e$message, call. = FALSE)
        }
      )
    }
  } else {
    chunk_results <- list(simulate_chunk(seq_len(num_simulations), seed = NULL, tick = tick_fun))
  }

  aggregate_results(chunk_results)

  stopifnot(num_simulations > 0)
  inv_num_sims <- 1 / num_simulations
  
  for (j in seq_along(arm_names)) {
    arm <- arm_names[j]
    early_eff <- sum_stop_efficacy[arm] * inv_num_sims
    early_fut <- sum_stop_futility[arm] * inv_num_sims
    final_eff <- sum_final_efficacy[arm] * inv_num_sims
    final_fut <- sum_final_futility[arm] * inv_num_sims
    final_inc <- sum_final_inconclusive[arm] * inv_num_sims
    results_data$Type_I_Error_or_Power[j]   <- early_eff + final_eff
    results_data$PET_Efficacy[j]            <- early_eff
    results_data$PET_Futility[j]            <- early_fut
    results_data$Pr_Reach_Max_N[j]          <- 1 - early_eff - early_fut
    results_data$Pr_Final_Efficacy[j]       <- final_eff
    results_data$Pr_Final_Futility[j]       <- final_fut
    results_data$Pr_Final_Inconclusive[j]   <- final_inc
    results_data$Exp_N[j]                   <- sum_final_n[arm] * inv_num_sims
    results_data$Exp_Events[j]              <- sum_final_events[arm] * inv_num_sims
  }

  results_data
}


#' Build scenario overrides from a grid of design choices
#'
#' Expands a named list of options into a list of scenario override lists that
#' can be merged with `base_args` prior to simulation.
#'
#' @param choices Named list where each element is either a vector of scalar
#'   values or a list whose elements are per-arm vectors.
#'
#' @return A list of scenario override lists.  The underlying Cartesian grid is
#'   attached as an attribute named `"grid"`.
#'
#' @examples
#' scenarios_from_grid(list(
#'   max_total_patients_per_arm = list(
#'     c(Doublet = 60, Triplet = 60),
#'     c(Doublet = 60, Triplet = 70)
#'   ),
#'   compare_arms_futility_margin = c(0.3, 0.4)
#' ))
#'
#' @export
# ---------- 1) Build scenarios from a grid of options ----------
# Pass a named list where each element is either:
#   - a vector of scalar options (e.g., c(60, 70)), OR
#   - a list of per-arm options (e.g., list(c(Arm1=6,Arm2=6,Arm3=9), c(...)))
# Returns: a list of "override" lists, one per scenario.
scenarios_from_grid <- function(choices) {
  stopifnot(is.list(choices), length(choices) > 0)
  keys <- names(choices)
  # For list-valued choices (e.g., per-arm vectors), expand over their indices
  grid_input <- lapply(choices, function(v) if (is.list(v)) seq_along(v) else v)
  raw_grid <- do.call(expand.grid, c(grid_input, stringsAsFactors = FALSE))
  # Turn each grid row into an override list, replacing list indices by actual list elements
  scen_list <- lapply(seq_len(nrow(raw_grid)), function(i) {
    row <- raw_grid[i, , drop = FALSE]
    over <- setNames(vector("list", length(keys)), keys)
    for (k in keys) {
      v <- choices[[k]]
      if (is.list(v)) {
        over[[k]] <- v[[ as.integer(row[[k]]) ]]
      } else {
        over[[k]] <- row[[k]]
      }
    }
    over
  })
  # Attach a compact label for readability (optional)
  attr(scen_list, "grid") <- raw_grid
  scen_list
}

#' Evaluate a design across multiple scenarios
#'
#' Merges each scenario override list with `base_args`, runs
#' `run_simulation_pure()` for every scenario, and binds the results into a
#' single table.
#'
#' @param base_args Named list of arguments accepted by `run_simulation_pure()`.
#' @param scens List of scenario override lists, typically from
#'   `scenarios_from_grid()`.
#' @param parallel Logical; if `TRUE` uses `parallel::mclapply()` to distribute
#'   scenarios across cores.
#' @param seed Optional integer seed passed to `set.seed()` before simulations.
#'
#' @return A data.table/data.frame containing the combined operating
#'   characteristic summaries.  A `scenario` column identifies the originating
#'   scenario index.
#'
#' @examples
#' \dontrun{
#' base_args <- list(
#'   num_simulations = 200,
#'   arm_names = c("Doublet", "Triplet"),
#'   reference_arm_name = "Doublet",
#'   compare_arms_option = TRUE,
#'   weibull_shape_true_arms = c(Doublet = 1.2, Triplet = 1.2),
#'   weibull_median_true_arms = c(Doublet = 6, Triplet = 6),
#'   null_median_arms = c(Doublet = 6, Triplet = 6),
#'   futility_median_arms = c(Doublet = 6, Triplet = 6),
#'   interval_cutpoints_sim = seq(0, 24, by = 3),
#'   max_follow_up_sim = 24,
#'   censor_max_time_sim = 24,
#'   prior_alpha_params_model = rep(0.5, 8),
#'   prior_beta_params_model = rep(0.5, 8),
#'   num_posterior_draws = 400,
#'   cohort_size_per_arm = 1,
#'   max_total_patients_per_arm = c(Doublet = 60, Triplet = 60),
#'   min_patients_for_analysis = 10,
#'   efficacy_stopping_rule_hc = TRUE,
#'   efficacy_threshold_current_prob_hc = 0.95,
#'   posterior_futility_threshold_hc = 0.8,
#'   futility_stopping_rule_hc = TRUE,
#'   efficacy_threshold_vs_ref_prob = 0.98,
#'   futility_threshold_vs_ref_prob = 0.6,
#'   compare_arms_futility_margin = 0.4,
#'   overall_accrual_rate = 3,
#'   randomization_probs = c(Doublet = 1/3, Triplet = 2/3),
#'   min_follow_up_at_final = 0,
#'   min_events_for_analysis = 0,
#'   min_median_followup = 0,
#'   interim_calendar_beat = 3,
#'   pred_success_pp_threshold_hc = 1,
#'   pred_futility_pp_threshold_hc = 0,
#'   num_posterior_draws_pred = 200
#' )
#'
#' scens <- scenarios_from_grid(list(
#'   weibull_median_true_arms = list(
#'     c(Doublet = 6, Triplet = 6),
#'     c(Doublet = 6, Triplet = 9)
#'   )
#' ))
#'
#' run_scenarios(base_args, scens, parallel = FALSE, seed = 123)
#' }
#'
#' @export
run_scenarios <- function(base_args, scens, parallel = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # keep only args that run_simulation_pure actually accepts
  rs_formals <- names(formals(run_simulation_pure))
  
  run_one <- function(i) {
    over   <- scens[[i]]
    args_i <- utils::modifyList(base_args, over, keep.null = TRUE)
    args_i <- args_i[intersect(names(args_i), rs_formals)]
    res <- do.call(run_simulation_pure, args_i)
    res$scenario <- i
    res
  }
  
  if (isTRUE(parallel)) {
    cores <- max(1L, parallel::detectCores() - 1L)
    out <- parallel::mclapply(seq_along(scens), run_one, mc.cores = cores)
  } else {
    out <- lapply(seq_along(scens), run_one)
  }
  
  # validate all returned items are tabular
  ok_types <- vapply(out, function(x) is.data.frame(x) || is.list(x) || data.table::is.data.table(x), logical(1))
  if (!all(ok_types)) stop("run_scenarios: one or more scenarios did not return a tabular result.")
  
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

# helper to build new interval cutpoints once enough events accrued
rebalance_cutpoints_from_events <- function(event_times, max_follow_up, num_intervals) {
  if (length(event_times) < max(2, num_intervals)) {
    return(NULL)
  }
  event_times <- pmin(pmax(event_times, 0), max_follow_up)
  probs <- seq(0, 1, length.out = num_intervals + 1)
  q_vals <- as.numeric(stats::quantile(event_times, probs = probs, na.rm = TRUE, type = 7))
  q_vals[1] <- 0
  q_vals[length(q_vals)] <- max_follow_up
  for (i in 2:length(q_vals)) {
    if (!is.finite(q_vals[i]) || q_vals[i] <= q_vals[i - 1]) {
      q_vals[i] <- min(max_follow_up, q_vals[i - 1] + 1e-6)
    }
  }
  if (any(diff(q_vals) <= 0)) {
    return(NULL)
  }
  q_vals
}


#' Create a reusable evolveTrial simulation cluster
#'
#' Creates a parallel cluster optimized for evolveTrial simulations. The cluster
#' can be reused across multiple calls to `run_simulation_pure()` by passing it
#' via the `cluster` parameter, avoiding the overhead of repeated cluster
#' creation/destruction.
#'
#' @param workers Integer; number of worker processes. Defaults to
#'   `parallel::detectCores() - 1`.
#' @param cluster_type One of `"auto"` (default), `"FORK"`, or `"PSOCK"`. Auto
#'   selects FORK on Unix systems (faster) and PSOCK on Windows.
#'
#' @return A parallel cluster object that can be passed to `run_simulation_pure()`
#'   via the `cluster` parameter.
#'
#' @details This function is particularly useful for Bayesian optimization
#'   workflows where `run_simulation_pure()` is called hundreds of times. By
#'   creating the cluster once and reusing it, you can eliminate the 5-10 second
#'   cluster spawn/teardown overhead per call.
#'
#'   FORK clusters (default on Unix) are significantly faster because:
#'   - Worker processes inherit the parent environment (no package loading)
#'   - Data is shared via copy-on-write (no serialization overhead)
#'
#'   Remember to call `release_cluster()` when done to free resources.
#'
#' @examples
#' \dontrun{
#' # Create cluster once
#' cl <- create_simulation_cluster(workers = 8)
#'
#' # Use in repeated BO evaluations
#' for (i in 1:100) {
#'   result <- run_simulation_pure(
#'     num_simulations = 500,
#'     ...,
#'     parallel_replicates = TRUE,
#'     cluster = cl
#'   )
#' }
#'
#' # Clean up
#' release_cluster(cl)
#' }
#'
#' @seealso [release_cluster()], [run_simulation_pure()]
#' @export
create_simulation_cluster <- function(workers = NULL, cluster_type = c("auto", "FORK", "PSOCK")) {
  cluster_type <- match.arg(cluster_type)

  # Auto-detect optimal cluster type
  if (cluster_type == "auto") {
    cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  }

  # Default workers
  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores() - 1L)
  }
  workers <- as.integer(workers)

  # Create cluster
  cl <- parallel::makeCluster(workers, type = cluster_type)

  # For PSOCK clusters, pre-load evolveTrial package in workers
  if (cluster_type == "PSOCK") {
    pkg_name <- utils::packageName()
    if (is.null(pkg_name) || identical(pkg_name, "")) {
      pkg_name <- "evolveTrial"
    }
    parallel::clusterCall(cl, function(pkg) {
      suppressPackageStartupMessages(require(pkg, character.only = TRUE))
      NULL
    }, pkg_name)
  }

  # Add class for identification
  class(cl) <- c("evolveTrial_cluster", class(cl))

  message("Created ", cluster_type, " cluster with ", workers, " workers")
  cl
}


#' Release an evolveTrial simulation cluster
#'
#' Stops and releases resources for a cluster created by
#' `create_simulation_cluster()`.
#'
#' @param cluster A cluster object created by `create_simulation_cluster()` or
#'   `parallel::makeCluster()`.
#'
#' @return NULL invisibly.
#'
#' @seealso [create_simulation_cluster()]
#' @export
release_cluster <- function(cluster) {
  if (!is.null(cluster)) {
    tryCatch(
      parallel::stopCluster(cluster),
      error = function(e) {
        warning("Failed to stop cluster: ", e$message)
      }
    )
  }
  invisible(NULL)
}
