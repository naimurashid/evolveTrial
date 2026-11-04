# simulate_chunk_cpp.R -- R side of the C++ chunk engine (engine C++ port, steps 3 and 4).
#
# run_simulation_pure() packs the closure's free variables once per frame with .pack_chunk_cfg(), defines
# simulate_chunk_cpp() next to simulate_chunk_r(), and makes simulate_chunk() a dispatcher: the C++ body
# (src/simulate_chunk.cpp) when getOption("evolveTrial.use_cpp_engine", TRUE) is TRUE and .chunk_cpp_supported() holds,
# otherwise the R body. simulate_chunk_cpp() itself never falls back: outside the supported set it stops.
#
# The supported set is the union of two frozen branch sets:
# - single-arm (HC analysis): the R01 to R04 acceptance designs (dev/plans/2026-11-03-step3-brief.md, Section 2);
# - between-arm (vs-reference, PH model): the R23 to R28 acceptance designs (dev/plans/2026-11-03-step4-brief.md,
#   Section 3).
# A frame is checked against the common items plus the HC items when compare_arms_option is not TRUE in the frame, and
# against the common items plus the vs-ref items when it is (inventory dev/perf/evidence/chunk_body_inventory.md,
# Section 2). The HC list is exactly the step 3 list (same labels, same checks, same order).

# Items shared by both paths (step 3 items 4 to 8, 12, 13b and 13g). Item 11 (EVOLVETRIAL_USE_CPP) is read at dispatch
# by .chunk_cpp_env_ok(). Item 13e has a path-specific form, see the two lists below.
.chunk_cpp_common_checks <- function(env) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  is_num1 <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  is_fin1 <- function(x) is_num1(x) && is.finite(x)
  args <- g("base_args_for_interim")
  arm_names <- g("arm_names")
  list(
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
    # 12. no progress bar (P13)
    "12 show_progress FALSE" = function() identical(g("show_progress"), FALSE),
    "13b diagnostics FALSE (3.5)" = function() isFALSE(g("diagnostics")),
    # 13g. looked up by arm name for every arm in arm_names (one arm on the HC path, two on the vs-ref path). The prior
    #     shapes and rates (prior_alpha_params_model, prior_beta_params_model) are the only prior parameters either C++
    #     engine reads that are not clamped: they must be finite and > 0 (Codex review P2: negative shapes made the exact
    #     engine report inconclusive where the R body and phase 1 report efficacy). The log-HR prior sd is clamped to
    #     max(1e-6, sd) by the R body and by .pack_vsref() alike, and the prior mean is a location; V10 keeps both finite.
    "13g well-formed design inputs" = function() {
      cuts <- g("interval_cutpoints_sim"); K <- length(cuts) - 1L
      pa <- g("prior_alpha_params_model"); pb <- g("prior_beta_params_model")
      rp <- g("randomization_probs"); sh <- g("weibull_shape_true_arms"); sc <- g("weibull_scale_true_arms")
      is.numeric(cuts) && K >= 1L && all(is.finite(cuts)) && !is.unsorted(cuts) &&
        identical(args$interval_cutpoints_sim, cuts) &&
        is.numeric(pa) && is.numeric(pb) && length(pa) == K && length(pb) == K && all(is.finite(pa)) && all(is.finite(pb)) &&
        all(pa > 0) && all(pb > 0) &&
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

# The single-arm (HC) items of step 3: 1 to 3, 9, 10, 13a, 13c to 13f.
.chunk_cpp_hc_checks <- function(env) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  is_num1 <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  is_fin1 <- function(x) is_num1(x) && is.finite(x)
  args <- g("base_args_for_interim")
  arm_names <- g("arm_names")
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
    # 13. other single-arm branches the packed configuration can reach, pinned to the frozen branch or value
    "13a min_person_time_frac_per_arm is the frozen unnamed 0 (2.49)" = function() {
      x <- args$min_person_time_frac_per_arm
      is.numeric(x) && length(x) == 1L && is.null(names(x)) && identical(as.numeric(x), 0)
    },
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
    }
  )
}

# The between-arm (vs-reference, PH model) items of step 4, V1 to V11, plus the vs-ref form of 13e.
.chunk_cpp_vsref_checks <- function(env) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  is_num1 <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  is_fin1 <- function(x) is_num1(x) && is.finite(x)
  null_or_fin1 <- function(x) is.null(x) || (is.numeric(x) && length(x) == 1L && is.null(names(x)) && is.finite(x))
  args <- g("base_args_for_interim")
  arm_names <- g("arm_names")
  list(
    # V1. vs-ref path, in the frame (SD:561, 640, 737, 758) and in args (IL:78) (2.35)
    "V1 compare_arms_option TRUE in the frame and in args" = function() {
      isTRUE(g("compare_arms_option")) && isTRUE(args$compare_arms_option)
    },
    # V2. PH model at the interim (IL:10, PH:96) and at the final (SD:804) (2.29, 2.42, 2.61)
    "V2 use_ph_model_vs_ref TRUE" = function() isTRUE(args$use_ph_model_vs_ref),
    # V3. two arms, the reference one of them, exactly one experimental arm; the interim reads args$arm_names and
    #     args$reference_arm_name (IL:85-90, SM:154-162), the body the frame copies (SD:562, 738, 795). Refuses P6.
    "V3 two arms, one reference, one experimental" = function() {
      ref <- g("reference_arm_name")
      is.character(arm_names) && length(arm_names) == 2L && !anyNA(arm_names) && !anyDuplicated(arm_names) &&
        identical(args$arm_names, arm_names) &&
        is.character(ref) && length(ref) == 1L && !is.na(ref) && ref %in% arm_names &&
        identical(args$reference_arm_name, ref) && sum(arm_names != ref) == 1L
    },
    # V4. names(max_total_patients_per_arm) == arm_names (P19, 6.10): is_eligible is positional (SD:538-540), the
    #     vs-ref breaks (SD:564-565, 643-644) and the gates (SM:228-236) index by name
    "V4 max_total_patients_per_arm named by arm_names" = function() {
      m <- g("max_total_patients_per_arm")
      is.numeric(m) && identical(names(m), arm_names) && all(is.finite(m)) && all(m >= 0) &&
        identical(args$max_total_patients_per_arm, m)
    },
    # V5. no interim efficacy stopping (IL:147-149, 2.44)
    "V5 efficacy_stopping_rule_vs_ref not TRUE" = function() {
      !isTRUE(args$efficacy_stopping_rule_vs_ref) && !isTRUE(g("efficacy_stopping_rule_vs_ref"))
    },
    # V6. futility HR margin set, efficacy HR margin NULL, both read through args (IL:13-28, SD:806-815) (2.30, 2.31, 2.43)
    "V6 compare_arms_hr_margin one finite number, compare_arms_hr_eff_margin NULL" = function() {
      is_fin1(args$compare_arms_hr_margin) && is.null(args$compare_arms_hr_eff_margin)
    },
    # V7. the probability thresholds: the interim reads args (IL:142-159), the final the frame (SD:825, 827) (2.32, 2.45)
    "V7 vs-ref probability thresholds one finite number" = function() {
      is_fin1(args$efficacy_threshold_vs_ref_prob) && is_fin1(args$futility_threshold_vs_ref_prob) &&
        is_fin1(g("efficacy_threshold_vs_ref_prob")) && is_fin1(g("futility_threshold_vs_ref_prob"))
    },
    # V8. no target-arm override (IL:105-109, 2.38). args$ partial matching is what IL:105 reads, so both are checked.
    "V8 compare_arms_target_arm absent from args" = function() {
      !("compare_arms_target_arm" %in% names(args)) && is.null(args$compare_arms_target_arm)
    },
    # V9. vs-ref gate values are NULL or one unnamed finite number, so resolve_gate_vec takes its NULL or "recycle first
    #     element" branch (GD:23-41) and then the proportional scaling (GD:44-68) (2.60, 6.7)
    "V9 vs-ref gate values NULL or one unnamed finite number" = function() {
      null_or_fin1(args$min_events_per_arm) && null_or_fin1(args$min_median_followup_per_arm) &&
        null_or_fin1(args$min_person_time_frac_per_arm)
    },
    # V10. log-HR prior (PH:113-114, 2.63). A longer vector would make coalesce_num's is.na() error.
    "V10 ph_loghr_prior_mean and ph_loghr_prior_sd NULL or one finite number" = function() {
      null_or_fin1(args$ph_loghr_prior_mean) && null_or_fin1(args$ph_loghr_prior_sd)
    },
    # 13e (vs-ref form). The vs-ref interim and final both draw num_posterior_draws (IL:7, SD:802, 6.8); the interim
    #     count is not read on this path.
    "13e num_posterior_draws at least 1 (6.8)" = function() {
      nd <- g("num_posterior_draws")
      is_fin1(nd) && nd >= 1 && identical(args$num_posterior_draws, nd)
    },
    # V11. other vs-ref branches the packed configuration can reach, pinned to the frozen branch or value
    "V11a num_posterior_draws integral (PH:124, 127)" = function() {
      nd <- g("num_posterior_draws")
      is_fin1(nd) && nd == floor(nd) && nd <= .Machine$integer.max
    }
  )
}

# The static items for a frame: the step 3 list on the HC path, the common items plus V1 to V11 on the vs-ref path.
.chunk_cpp_static_checks <- function(env) {
  vsref <- tryCatch(isTRUE(get("compare_arms_option", envir = env, inherits = FALSE)), error = function(e) FALSE)
  cm <- .chunk_cpp_common_checks(env)
  if (vsref) return(c(cm, .chunk_cpp_vsref_checks(env)))
  hc <- .chunk_cpp_hc_checks(env)
  # step 3 order: 1, 2, 3, 4 to 8, 9, 10, 12, 13a, 13b, 13c, 13d, 13e, 13f, 13g
  c(hc[1:3], cm[1:5], hc[4:5], cm[6], hc[6], cm[7], hc[7:10], cm[8])
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
# frame. Draws nothing. Values the R body reads through args_local (the interim, IL:1-263, and on the vs-ref path the
# final's sampler and margins) are taken from base_args_for_interim with the same expressions; values the R body reads
# from the frame are taken from the frame. Outside the supported set it returns list(supported = FALSE, violations = ...)
# and never errors.
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
  vsref <- isTRUE(g("compare_arms_option"))

  common <- list(
    supported = TRUE,
    mode = if (vsref) "vsref" else "hc",
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
    draws_final = as.integer(g("num_posterior_draws")),
    collect_raw_n = isTRUE(g("collect_raw_n")),
    return_variance = isTRUE(g("return_variance"))
  )
  c(common, if (vsref) .pack_vsref(env, args, arm_names) else .pack_hc(env, args, arm_names))
}

# HC interim gates, rule and thresholds (step 3), exactly as IL:43-66 and IL:182-253 form them.
.pack_hc <- function(env, args, arm_names) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
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
    draws_interim = num_draws_interim,
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
    final_futility_prob = as.numeric(g("final_futility_posterior_prob_threshold"))
  )
}

# vs-ref interim gates (SM:153-258), PH sampler constants (PH:105-135, PD:98-118), margins and thresholds (IL:1-169,
# SD:794-833). Every constant the R body computes from the frame or args with a fixed expression is computed here with the
# same expression, once per frame, so the C++ receives the identical double.
.pack_vsref <- function(env, args, arm_names) {
  g <- function(name) get(name, envir = env, inherits = FALSE)
  ref_name <- g("reference_arm_name")
  trt_name <- arm_names[arm_names != ref_name]

  # gates_pass_for_both_arms (SM:153-258) with trt_name supplied by interim_check_vs_ref (IL:127)
  ctrl_name <- args$reference_arm_name
  gate_arms <- c(ctrl_name, trt_name)
  min_ev_vec <- resolve_gate_vec(args$min_events_per_arm, target_arms = gate_arms, all_arm_names = args$arm_names,
                                 randomization_probs = args$randomization_probs, default = 0, scale = TRUE)
  min_mfu_vec <- resolve_gate_vec(args$min_median_followup_per_arm, target_arms = gate_arms,
                                  all_arm_names = args$arm_names, randomization_probs = args$randomization_probs,
                                  default = 0, scale = FALSE)
  min_pt_vec <- resolve_gate_vec(args$min_person_time_frac_per_arm, target_arms = gate_arms,
                                 all_arm_names = args$arm_names, randomization_probs = args$randomization_probs,
                                 default = 0, scale = TRUE)
  mtn <- args$max_total_patients_per_arm
  max_follow <- coalesce_num(args$max_trial_time, coalesce_num(args$max_follow_up_sim, 0))

  # sample_vs_ref_medians_ph (PH:113-114) and the ph_beta_mode_var dispatcher's defaults (PD:98-99)
  ph_formals <- formals(ph_beta_mode_var)

  # IL:11-28 (interim, from args) and SD:805-815 (final, from args_local, the same list)
  log_margin_interim <- 0
  if (!is.null(args$compare_arms_hr_margin)) log_margin_interim <- log1p(-min(max(0, args$compare_arms_hr_margin), 0.95))
  log_eff_interim <- 0
  if (!is.null(args$compare_arms_hr_eff_margin)) log_eff_interim <- log1p(-min(max(0, args$compare_arms_hr_eff_margin), 0.95))
  args_local <- args
  log_margin_final <- 0
  if (!is.null(args_local$compare_arms_hr_margin)) {
    log_margin_final <- log1p(-min(max(0, args_local$compare_arms_hr_margin), 0.95))
  }
  log_eff_final <- 0
  if (!is.null(args_local$compare_arms_hr_eff_margin)) {
    log_eff_final <- log1p(-min(max(0, args_local$compare_arms_hr_eff_margin), 0.95))
  }

  list(
    ref_idx = match(ref_name, arm_names) - 1L,          # 0-based
    trt_idx = match(trt_name, arm_names) - 1L,
    vr_min_ev = as.numeric(c(min_ev_vec[[ctrl_name]] %||% 0, min_ev_vec[[trt_name]] %||% 0)),       # ctrl, trt
    vr_min_mfu = as.numeric(c(min_mfu_vec[[ctrl_name]] %||% 0, min_mfu_vec[[trt_name]] %||% 0)),
    vr_min_pt = as.numeric(c(min_pt_vec[[ctrl_name]] %||% 0, min_pt_vec[[trt_name]] %||% 0)),
    vr_max_pt = as.numeric(c(coalesce_num(mtn[[ctrl_name]], 0) * max_follow, coalesce_num(mtn[[trt_name]], 0) * max_follow)),
    draws_vsref_interim = as.integer(args$num_posterior_draws),   # IL:7
    ph_mu = as.numeric(coalesce_num(args$ph_loghr_prior_mean, 0)),
    ph_sigma = as.numeric(max(1e-6, coalesce_num(args$ph_loghr_prior_sd, 1))),
    ph_tol = as.numeric(eval(ph_formals$tol)),
    ph_max_iter = as.integer(eval(ph_formals$max_iter)),
    log_margin_interim = as.numeric(log_margin_interim),
    log_eff_interim = as.numeric(log_eff_interim),
    log_margin_final = as.numeric(log_margin_final),
    log_eff_final = as.numeric(log_eff_final),
    vr_fut_prob_interim = as.numeric(args$futility_threshold_vs_ref_prob),   # IL:158-159
    vr_eff_prob_final = as.numeric(g("efficacy_threshold_vs_ref_prob")),     # SD:825
    vr_fut_prob_final = as.numeric(g("futility_threshold_vs_ref_prob"))      # SD:827
  )
}

# Test hook (step 4): sample_vs_ref_medians_ph through the C++ port, for a unit test against the R sampler on slices.
# slCtrl and slTrt are slice_arm_data_at_time() results. Returns list(logHR, medCtrl, medTrt) like the R function.
.chunk_ph_sample_r <- function(slCtrl, slTrt, args, num_samples) {
  ph_formals <- formals(ph_beta_mode_var)
  .chunk_ph_sample(
    as.integer(slCtrl$metrics$events_per_interval), as.numeric(slCtrl$metrics$person_time_per_interval),
    as.integer(slTrt$metrics$events_per_interval), as.numeric(slTrt$metrics$person_time_per_interval),
    as.numeric(args$prior_alpha_params_model), as.numeric(args$prior_beta_params_model),
    as.numeric(diff(args$interval_cutpoints_sim)),
    as.numeric(coalesce_num(args$ph_loghr_prior_mean, 0)), as.numeric(max(1e-6, coalesce_num(args$ph_loghr_prior_sd, 1))),
    as.numeric(eval(ph_formals$tol)), as.integer(eval(ph_formals$max_iter)), as.integer(num_samples)
  )
}
