# simulate_chunk_cpp.R -- R side of the C++ chunk engine (engine C++ port, step 3).
#
# run_simulation_pure() packs the closure's free variables once per frame with .pack_chunk_cfg(), defines
# simulate_chunk_cpp() next to simulate_chunk_r(), and makes simulate_chunk() a dispatcher: the C++ body
# (src/simulate_chunk.cpp) when getOption("evolveTrial.use_cpp_engine", TRUE) is TRUE and .chunk_cpp_supported() holds,
# otherwise the R body. simulate_chunk_cpp() itself never falls back: outside the supported set it stops.
#
# The supported set is the frozen single-arm branch set of the R01 to R04 acceptance designs
# (dev/plans/2026-11-03-step3-brief.md, Section 2; inventory dev/perf/evidence/chunk_body_inventory.md, Section 2).

# Item labels, in the order of the brief. Each check returns TRUE when the item holds; an error counts as a violation.
.chunk_cpp_static_checks <- function(env) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  is_num1 <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  is_fin1 <- function(x) is_num1(x) && is.finite(x)
  args <- g("base_args_for_interim")
  arm_names <- g("arm_names")
  arm1 <- function() arm_names[[1L]]
  named_finite_has_arm <- function(x) {
    is.numeric(x) && !is.null(names(x)) && all(arm_names %in% names(x)) && all(is.finite(x[arm_names]))
  }
  list(
    # 1. HC path (inventory 2.10, 2.20, 2.22, 2.35)
    "1 compare_arms_option not TRUE" = function() !isTRUE(g("compare_arms_option")) && !isTRUE(args$compare_arms_option),
    # 2. single arm
    "2 one arm" = function() is.character(arm_names) && length(arm_names) == 1L && !is.na(arm_names) &&
      identical(args$arm_names, arm_names),
    # 3. names(max_total_patients_per_arm) == arm_names (P19, 6.10), so IL:203-204 takes the named branch (2.50)
    "3 max_total_patients_per_arm named by arm_names" = function() {
      m <- g("max_total_patients_per_arm")
      is.numeric(m) && identical(names(m), arm_names) && all(is.finite(m)) && all(m >= 0) &&
        identical(args$max_total_patients_per_arm, m)
    },
    # 4. looks_stop_at_final (2.5, 2.13, I7)
    "4 looks_stop_at_final TRUE" = function() isTRUE(g("looks_stop_at_final")),
    # 5. no person-time milestones (2.7)
    "5 pt_targets_abs NULL" = function() is.null(g("pt_targets_abs")),
    # 6. dropout off (2.11, P3)
    "6 random_dropout FALSE" = function() identical(g("random_dropout"), FALSE),
    # 7. finite max_trial_time (2.14, 2.18)
    "7 max_trial_time finite" = function() is_fin1(g("max_trial_time")) && identical(args$max_trial_time, g("max_trial_time")),
    # 8. rebalancing threshold non-NULL (2.4)
    "8 rebalance_after_events non-NULL" = function() is_num1(g("rebalance_after_events")),
    # 9. no interim efficacy stopping (2.53)
    "9 efficacy_stopping_rule_hc not TRUE" = function() !isTRUE(args$efficacy_stopping_rule_hc) &&
      !isTRUE(g("efficacy_stopping_rule_hc")),
    # 10. median thresholds named, containing the arm, finite: at the final from the frame (2.23, 2.24, 2.26), and at the
    #     interim as IL:57-61 resolves them from args (2.55). args (SD:377-428) has no median_pfs_futility_threshold_arms
    #     element, so the interim futility threshold resolves to args$futility_median_arms (or the success vector).
    "10 success and futility median thresholds named and finite" = function() {
      s <- g("median_pfs_success_threshold_arms"); f <- g("median_pfs_futility_threshold_arms")
      iv <- .chunk_cpp_interim_thresholds(args)
      named_finite_has_arm(s) && named_finite_has_arm(f) &&
        named_finite_has_arm(iv$success) && named_finite_has_arm(iv$futility)
    },
    # 12. no progress bar (P13)
    "12 show_progress FALSE" = function() identical(g("show_progress"), FALSE),
    # 13. other single-arm branches the packed configuration can reach, pinned to the frozen branch or value
    "13a min_person_time_frac_per_arm is the frozen unnamed 0 (2.49)" = function() {
      x <- args$min_person_time_frac_per_arm
      is.numeric(x) && length(x) == 1L && is.null(names(x)) && identical(as.numeric(x), 0)
    },
    "13b diagnostics FALSE (3.5)" = function() isFALSE(g("diagnostics")),
    "13c interim futility_threshold_hc_prob one finite number (2.54)" = function() {
      x <- args$futility_threshold_hc_prob
      is.numeric(x) && length(x) == 1L && is.finite(x)
    },
    "13d final posterior probability thresholds finite (2.25, 2.27)" = function() {
      is_fin1(g("final_success_posterior_prob_threshold")) && is_fin1(g("final_futility_posterior_prob_threshold"))
    },
    "13e posterior draw counts at least 1 (2.56, 6.8)" = function() {
      nd <- g("num_posterior_draws"); ndi <- args$num_posterior_draws_interim %||% args$num_posterior_draws
      is_fin1(nd) && nd >= 1 && is_fin1(ndi) && identical(args$num_posterior_draws, nd)
    },
    "13f HC gate values are scalars (2.48)" = function() {
      ok1 <- function(x) is.null(x) || (is.numeric(x) && length(x) == 1L)
      ok1(args$min_patients_for_analysis) && ok1(args$min_events_hc) && ok1(args$min_median_followup_hc)
    },
    "13g well-formed design inputs" = function() {
      cuts <- g("interval_cutpoints_sim"); K <- length(cuts) - 1L
      pa <- g("prior_alpha_params_model"); pb <- g("prior_beta_params_model")
      rp <- g("randomization_probs"); sh <- g("weibull_shape_true_arms"); sc <- g("weibull_scale_true_arms")
      is.numeric(cuts) && K >= 1L && all(is.finite(cuts)) && !is.unsorted(cuts) &&
        identical(args$interval_cutpoints_sim, cuts) &&
        is.numeric(pa) && is.numeric(pb) && length(pa) == K && length(pb) == K && all(is.finite(pa)) && all(is.finite(pb)) &&
        identical(args$prior_alpha_params_model, pa) && identical(args$prior_beta_params_model, pb) &&
        is.numeric(rp) && all(is.finite(rp[arm_names])) && all(rp[arm_names] > 0) && identical(args$randomization_probs, rp) &&
        all(is.finite(sh[arm_names])) && all(sh[arm_names] > 0) && all(is.finite(sc[arm_names])) && all(sc[arm_names] > 0) &&
        is_fin1(g("overall_accrual_rate")) && g("overall_accrual_rate") > 0 &&
        is_fin1(g("interim_calendar_beat")) && g("interim_calendar_beat") > 0 &&
        is_fin1(g("min_follow_up_at_final")) &&
        is_fin1(g("max_follow_up_sim")) && g("max_follow_up_sim") > 0 &&
        identical(args$max_follow_up_sim, g("max_follow_up_sim")) &&
        is_fin1(g("censor_max_time_sim")) && g("censor_max_time_sim") > 0
    }
  )
}

# IL:57-59 verbatim: the vectors calculate_current_probs_hc() indexes by arm name at an interim.
.chunk_cpp_interim_thresholds <- function(args) {
  success_thr_vec <- args$median_pfs_success_threshold_arms %||% args$null_median_arms
  futility_thr_vec <- args$median_pfs_futility_threshold_arms %||%
    args$futility_median_arms %||% success_thr_vec
  list(success = success_thr_vec, futility = futility_thr_vec)
}

# The labels of the static items that fail for a frame (character(0) when all hold).
.chunk_cpp_violations <- function(env) {
  checks <- .chunk_cpp_static_checks(env)
  ok <- vapply(checks, function(f) isTRUE(tryCatch(f(), error = function(e) FALSE)), logical(1))
  names(checks)[!ok]
}

# Item 11: EVOLVETRIAL_USE_CPP in the set the posterior dispatchers treat as TRUE (P14, inventory 2.69). Read at dispatch.
.chunk_cpp_env_ok <- function() .use_cpp_posterior()

# The supported-set predicate: every static item holds for the frame and EVOLVETRIAL_USE_CPP is truthy now.
.chunk_cpp_supported <- function(env) {
  length(.chunk_cpp_violations(env)) == 0L && isTRUE(.chunk_cpp_env_ok())
}

# Called by simulate_chunk_cpp() before anything else (before set.seed), so a refusal leaves the stream untouched.
.chunk_cpp_assert_supported <- function(env, cfg) {
  v <- .chunk_cpp_violations(env)
  if (!isTRUE(.chunk_cpp_env_ok())) v <- c(v, "11 EVOLVETRIAL_USE_CPP is TRUE, 1 or YES")
  if (length(v) > 0L || !isTRUE(cfg$supported)) {
    stop("simulate_chunk_cpp: configuration outside the supported set of the C++ chunk engine (it never falls back; ",
         "use simulate_chunk or simulate_chunk_r). Failing items: ",
         if (length(v)) paste(v, collapse = "; ") else "packed configuration marked unsupported", call. = FALSE)
  }
  invisible(TRUE)
}

# Packs the closure's free variables (inventory Section 1) into the list src/simulate_chunk.cpp reads. Built once per
# frame. Draws nothing. Values the R body reads through args_local (the HC interim, IL:171-263) are taken from
# base_args_for_interim with the same expressions; values the R body reads from the frame are taken from the frame.
# Outside the supported set it returns list(supported = FALSE, violations = ...) and never errors.
.pack_chunk_cfg <- function(env) {
  force(env$looks_stop_at_final)   # plan I7: the lazy default is forced here, once, in the master
  v <- .chunk_cpp_violations(env)
  if (length(v) > 0L) return(list(supported = FALSE, violations = v))
  out <- tryCatch(.pack_chunk_cfg_supported(env), error = function(e) {
    list(supported = FALSE, violations = paste("packing failed:", conditionMessage(e)))
  })
  out
}

.pack_chunk_cfg_supported <- function(env) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  args <- g("base_args_for_interim")
  arm_names <- g("arm_names")
  mtpa <- g("max_total_patients_per_arm")
  cuts <- g("interval_cutpoints_sim")

  # HC gate thresholds, exactly as IL:182-216 forms them
  per_arm <- lapply(arm_names, function(arm) {
    min_pt_frac_vec <- resolve_gate_vec(
      args$min_person_time_frac_per_arm,
      target_arms = arm,
      all_arm_names = args$arm_names,
      randomization_probs = args$randomization_probs,
      default = 0,
      scale = TRUE
    )
    m <- args$max_total_patients_per_arm
    max_follow <- coalesce_num(args$max_trial_time, coalesce_num(args$max_follow_up_sim, 0))
    list(
      min_pat = as.numeric(coalesce_num(args$min_patients_for_analysis, 0)),
      min_events = as.numeric(coalesce_num(args$min_events_hc, 0)),
      min_median_fu = as.numeric(coalesce_num(args$min_median_followup_hc, 0)),
      min_pt_frac = as.numeric(min_pt_frac_vec[[arm]] %||% 0),
      max_pt_arm = as.numeric(coalesce_num(m[[arm]], 0) * max_follow)
    )
  })
  pa <- function(field) vapply(per_arm, function(x) x[[field]], numeric(1))

  # interim draw count, IL:45-46
  num_draws_interim <- args$num_posterior_draws_interim %||% args$num_posterior_draws
  num_draws_interim <- max(1L, as.integer(num_draws_interim))

  list(
    supported = TRUE,
    arm_names = arm_names,
    mtpa = as.numeric(mtpa),
    capacity = as.integer(ceiling(mtpa * 1.2)),
    rand_probs = as.numeric(g("randomization_probs")[arm_names]),
    shape = as.numeric(g("weibull_shape_true_arms")[arm_names]),
    scale = as.numeric(g("weibull_scale_true_arms")[arm_names]),
    accrual_rate = as.numeric(g("overall_accrual_rate")),
    beat = as.numeric(g("interim_calendar_beat")),
    max_trial_time = as.numeric(g("max_trial_time")),
    min_follow_up_at_final = as.numeric(g("min_follow_up_at_final")),
    max_follow_up_sim = as.numeric(g("max_follow_up_sim")),
    censor_max_time_sim = as.numeric(g("censor_max_time_sim")),
    random_dropout = g("random_dropout"),
    looks_stop_at_final = g("looks_stop_at_final"),
    K = as.integer(g("num_intervals")),
    cut_base = as.numeric(cuts),
    len_base = as.numeric(g("interval_lengths_base")),
    rebalance_threshold = as.numeric(g("rebalance_after_events")),
    prior_alpha = as.numeric(g("prior_alpha_params_model")),
    prior_beta = as.numeric(g("prior_beta_params_model")),
    draws_interim = num_draws_interim,
    draws_final = as.integer(g("num_posterior_draws")),
    min_pat = pa("min_pat"),
    min_events = pa("min_events"),
    min_median_fu = pa("min_median_fu"),
    min_pt_frac = pa("min_pt_frac"),
    max_pt_arm = pa("max_pt_arm"),
    futility_prob_interim = as.numeric(args$futility_threshold_hc_prob),
    # interim futility threshold, IL:60-61: futility_thr_vec[[arm_name]] %||% futility_thr_vec[1]
    futility_thr_interim = vapply(arm_names, function(arm) {
      fv <- .chunk_cpp_interim_thresholds(args)$futility
      as.numeric(fv[[arm]] %||% fv[1])
    }, numeric(1), USE.NAMES = FALSE),
    # final thresholds, SD:770-772 and SD:781-783 (finite in the supported set, so the Inf / -Inf fallbacks are not taken)
    success_thr = as.numeric(g("median_pfs_success_threshold_arms")[arm_names]),
    futility_thr = as.numeric(g("median_pfs_futility_threshold_arms")[arm_names]),
    final_success_prob = as.numeric(g("final_success_posterior_prob_threshold")),
    final_futility_prob = as.numeric(g("final_futility_posterior_prob_threshold")),
    collect_raw_n = isTRUE(g("collect_raw_n")),
    return_variance = isTRUE(g("return_variance"))
  )
}
