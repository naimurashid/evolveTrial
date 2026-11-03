# simulate_chunk body inventory for the bit-identical C++ port

Worktree `~/research/evolveTrial-perf`, branch perf/engine-cpp, HEAD c2a88d5. Read-only survey. All
file:line references are to this worktree. `SD` = `R/simulation_driver.R`, `IL` = `R/interim_logic.R`,
`SM` = `R/state_management.R`, `PH` = `R/posterior_helpers.R`, `PD` = `R/posterior_cpp_dispatchers.R`,
`PS` = `src/posterior_sampling.cpp`, `GD` = `R/gate_diagnostics.R`, `DG` = `R/data_generation.R`.

Frozen settings come from `dev/perf/evidence/args_R01.txt` (single-arm, "R01") and `args_R23.txt`
(between-arm, "R23"), plus `acceptance_set.md` for settings the capture does not show. One setting is
NOT in the args files and decides a major branch: `looks_stop_at_final`. It is not a base_args entry. The
default is `getOption("evolveTrial.looks_stop_at_final", FALSE)` (SD:236), and the certify phase sets the
option to TRUE (`acceptance_set.md:63`, `scripts/capture.R:27`). So the frozen designs run with TRUE.

Empirical checks I ran (R 4.6.0, scratchpad scripts `chk.R`, `chk2.R`, `chk3.R`, nothing written to the
repo) are marked "checked". Statements about R internals that I did not execute are marked "from R source
memory, verify by a bit-compare test".

Tags used in section 2: FROZEN-R01, FROZEN-R23 (exercised by that design), BOTH, NOT EXERCISED (by
the frozen designs).

---------------------------------------------------------------------------------------------------

## 1. Closure free variables

`simulate_chunk(sim_indices, seed = NULL, tick = function() {})` is defined at SD:468. Its enclosing
frame is the `run_simulation_pure` call frame (SD:186). Because `parLapply` serialises the anonymous
function (SD:1042-1044) with its enclosure, the whole frame travels to the worker. Names that are not
in the frame resolve in the evolveTrial namespace.

### 1A. Variables read from the enclosing frame

| # | Name | Read at | Set at | Type and frozen value |
|---|---|---|---|---|
| 1.1 | `arm_names` | 472-478, 483, 492, 532, 563, 575, 642, 750, 836-838, 868, 876 | formal SD:187 | character; R01 "Experimental", R23 c("Doublet","Triplet"). Plain unnamed. |
| 1.2 | `collect_raw_n` | 482 | SD:445 (`isTRUE(return_percentiles) \|\| isTRUE(return_variance)`) | logical; TRUE in both (run_scenarios called with return_variance = TRUE) |
| 1.3 | `return_variance` | 491, 498 | formal SD:268, overwritten by run_scenarios at SD:1286 | logical TRUE |
| 1.4 | `interval_cutpoints_sim` | 525, 732 | formal SD:195 | double vector, length K+1. R01 c(0,3,...,24) (K=8). R23 args file shows c(0,2,...,24) (K=12), see 6.20. |
| 1.5 | `interval_lengths_base` | 527 | SD:431 `diff(interval_cutpoints_sim)` | double vector length K |
| 1.6 | `rebalance_after_events` (copied to `rebalance_threshold`) | 528 | formal SD:257 | NULL or numeric scalar; frozen 50 (double) |
| 1.7 | `base_args_for_interim` | 523 | SD:466, copy of `args` built at SD:377-428 | list of 46 named entries, see 1B |
| 1.8 | `max_total_patients_per_arm` | 532 (make_state), 539, 564-565, 643-644 | formal SD:203 | named vector; R01 c(Experimental=32L) integer, R23 c(Doublet=52, Triplet=51) double |
| 1.9 | `interim_calendar_beat` | 534, 631, 690, 711 | formal SD:243 | double scalar; R01 1.28536033253331, R23 1.00580743406046 |
| 1.10 | `overall_accrual_rate` | 543 | formal SD:229 | double scalar 2 |
| 1.11 | `looks_stop_at_final` | 552, 561, 681 | formal SD:236, validated SD:277-280 | logical; TRUE (frozen, via option) |
| 1.12 | `max_trial_time` | 552, 687-695, 701-704, 726-727 | formal SD:234 | double; 72 (default Inf) |
| 1.13 | `compare_arms_option` | 561, 640, 737, 758 | formal SD:190 | logical; R01 FALSE, R23 TRUE |
| 1.14 | `reference_arm_name` | 562, 563, 642, 738-740, 795 | formal SD:191 | NULL (R01) or "Doublet" (R23) |
| 1.15 | `pt_targets_abs` | 605, 606, 610-617, 625-627 | SD:364-368 | NULL in both (person_time_milestones NULL). When not NULL: `sort(unique(milestones)) * total_max_PT`. |
| 1.16 | `latest_calendar_look` | 613, 620 | formal SD:255 | Inf default (only read inside the pt-milestone branch) |
| 1.17 | `max_follow_up_sim` | 577, 740, 754, 849 | formal SD:197 | double 120 |
| 1.18 | `num_intervals` | 586 | SD:430 `length(interval_cutpoints_sim) - 1` | double (not integer); R01 8, R23 12 (capture) |
| 1.19 | `diagnostics` | 593, 615, 630, 682, 692, 714 | formal SD:244 | logical FALSE |
| 1.20 | `randomization_probs` | 650 | formal SD:230 | named double; R01 c(Experimental=1), R23 c(0.5, 0.5) |
| 1.21 | `weibull_shape_true_arms` | 655 | formal SD:191 | named double 1.3 |
| 1.22 | `weibull_scale_true_arms` | 656 | SD:329-332 `sapply(arm_names, calculate_weibull_scale(median[arm], shape[arm]))`, names reset at 332 | named double; `median / (log(2)^(1/shape))` (DG:50-55) |
| 1.23 | `censor_max_time_sim` | 657 | formal SD:197 | double 120 |
| 1.24 | `random_dropout` | 658 | formal SD:271, validated SD:274-276 | logical FALSE (frozen) |
| 1.25 | `min_follow_up_at_final` | 702, 704, 727, 729 | formal SD:231 | double 6 |
| 1.26 | `prior_alpha_params_model`, `prior_beta_params_model` | 763-764 | formals SD:198-199 | double vectors length K (0.2167... R01, 0.3051... R23) |
| 1.27 | `num_posterior_draws` | 765, 802 | formal SD:200 | double; R01 5000, R23 800 |
| 1.28 | `median_pfs_success_threshold_arms` | 770-771 | formal SD:224 | named double; R01 c(Experimental=3.8), R23 c(6,6) |
| 1.29 | `final_success_posterior_prob_threshold` | 777 | formal SD:225 | double; R01 0.893541169485304 |
| 1.30 | `median_pfs_futility_threshold_arms` | 781-782 | formal SD:226 | named double; R01 4.5273483339221 |
| 1.31 | `final_futility_posterior_prob_threshold` | 788 | formal SD:227 | double; R01 0.345920625284874 |
| 1.32 | `compare_arms_futility_margin` | 819 | formal SD:217 | double, default 0 (not in either base_args). Read only on the non-PH final path. |
| 1.33 | `efficacy_threshold_vs_ref_prob`, `futility_threshold_vs_ref_prob` | 825, 827 | formals SD:215-216 | double scalars; R23 0.950199307715971 / 0.372152871977103. R01 passes 0.98 / 0.6 but never reads them (HC final). |

Arguments of the closure itself: `sim_indices` (integer vector, used for `length()` at 483, 492, 498, 894
and for `which(sim_indices == sim_idx)` at 867, 875), `seed`, `tick`.

Variables NOT read by the closure although they look related: `args` (the closure reads the copy
`base_args_for_interim`), `total_max_PT`, `max_PT_per_arm` (only inside `args`), `num_draws_interim`
(only inside `args`), `cohort_size_per_arm`, `min_events_hc` and `min_median_followup_hc` (only inside
`args`), `show_progress`, `pb`, `tick_fun`, `all_final_*_raw`, `sum_*` accumulators.

### 1B. `args` / `base_args_for_interim` fields actually read downstream

`args_local <- base_args_for_interim` (SD:523) is passed to `interim_check` and to `sample_vs_ref_medians`
at the final step. Fields set at SD:377-428 and the ones overwritten per trial (`interval_cutpoints_sim`
at 526 and 590).

Fields read on the between-arm path (`interim_check_vs_ref`, IL:84-169, plus callees):
- `compare_arms_option` (IL:78), `reference_arm_name` (IL:85; SM:154, 157), `arm_names` (IL:90; SM:159, 162, 165)
- `compare_arms_target_arm` (IL:105). **Not a member of `args`**, so always NULL, the override branch is dead.
- `max_follow_up_sim` (IL:112, 125), `interval_cutpoints_sim` (IL:112, 125; PH:106; IL:44; PD:173, 244, 277)
- `use_ph_model_vs_ref` (IL:10, 117; PH:96), `num_posterior_draws` (IL:7, 118)
- gates (SM:182-211, 225-242): `min_events_per_arm`, `min_median_followup_per_arm`, `min_person_time_frac_per_arm`,
  `randomization_probs`, `max_total_patients_per_arm`, `max_trial_time`, `max_follow_up_sim`
- `compare_arms_hr_margin` (IL:13, 18), `compare_arms_hr_eff_margin` (IL:26-27), `compare_arms_futility_margin`
  (IL:35, 143; only non-PH and diagnostics)
- `efficacy_stopping_rule_vs_ref` (IL:147), `efficacy_threshold_vs_ref_prob` (IL:142, 148, 149),
  `futility_threshold_vs_ref_prob` (IL:143, 158, 159)
- `prior_alpha_params_model`, `prior_beta_params_model` (PH:111-112), `ph_loghr_prior_mean`, `ph_loghr_prior_sd` (PH:113-114)

Fields read on the single-arm path (`interim_check_hc`, IL:171-263, plus `calculate_current_probs_hc`, IL:43-66):
- `arm_names` (IL:172, 192, 206-207), `max_follow_up_sim`, `interval_cutpoints_sim` (IL:44, 176-177)
- `min_patients_for_analysis` (IL:182), `min_events_hc` (IL:183), `min_median_followup_hc` (IL:184)
- `min_person_time_frac_per_arm` (IL:188), `randomization_probs` (IL:191), `max_total_patients_per_arm` (IL:201-206), `max_trial_time` (IL:202)
- `num_posterior_draws_interim` (IL:45; falls back to `num_posterior_draws`), `prior_alpha_params_model`, `prior_beta_params_model` (IL:51-52)
- `median_pfs_success_threshold_arms` (IL:57), `null_median_arms` (IL:57, fallback), `median_pfs_futility_threshold_arms` (IL:58), `futility_median_arms` (IL:59, fallback)
- `efficacy_stopping_rule_hc` (IL:240), `efficacy_threshold_hc_prob` (IL:230, 241-243), `futility_threshold_hc_prob` (IL:234, 252-253)

Fields read at the final analysis through `args_local` (SD:798-808, 813): `use_ph_model_vs_ref`,
`compare_arms_hr_margin`, `compare_arms_hr_eff_margin`, and inside `sample_vs_ref_medians` the
interval cutpoints (rebalanced), priors, `ph_loghr_prior_*`.

Fields packed into `args` and never read anywhere downstream (dead): `predictive_fast`,
`num_posterior_draws_pred`, `predictive_boundary_band`, `run_pred_every_k_beats`,
`pred_success_pp_threshold_hc`, `pred_futility_pp_threshold_hc`, `censor_max_time_sim`,
`rebalance_after_events` (the closure uses its own copy), `max_PT_per_arm`,
`final_success_posterior_prob_threshold`, `final_futility_posterior_prob_threshold`,
`futility_stopping_rule_hc`, `futility_stopping_rule_vs_ref`, `median_pfs_*` on the vs-ref path.
(Verified by grep over SD, IL, SM, PH, PD, GD. Only `R/predictive_probabilities.R` reads the pred fields and nothing in the chunk path calls it.)

### 1C. Package-level functions the closure calls (namespace lookups at call time)

`get_active_registry` (SM:108), `make_state` (SM:70), `slice_arm_data_at_time` (SM:115),
`calculate_interval_metrics_fast` (SM:12), `interim_check` (IL:77), `rebalance_cutpoints_from_events` (SD:1331),
`cum_person_time_all_arms` (DG:26, NOT EXERCISED), `draw_posterior_hazard_samples` (PD:29),
`calculate_median_survival_matrix_cpp` (PS:151), `sample_vs_ref_medians` (PH:93), `coalesce_num` (SM:265),
`final_vsref_probs_abs` (PH:60, NOT EXERCISED), `%||%` (`R/globals.R:4`).
Because they are looked up at call time, `evolveTrial:::interim_check` and `run_scenarios` can be replaced
in the namespace. The certify phase does this (`acceptance_set.md:~62`: `verify_m4_install` wraps
`run_scenarios`, and wraps `interim_check` "for a probe only"). A C++ chunk bypasses any such wrapper.

Count: 33 rows in 1A, plus 1B (3 field groups, 8 dead-field groups) and 1C (12 functions).

---------------------------------------------------------------------------------------------------

## 2. Branches and options that depend on configuration

Notation: "R01" is the single-arm design (HC path), "R23" the between-arm design (vs-ref, PH).

### 2A. Inside `simulate_chunk` (SD:468-898)

| # | Line | Condition | R01 | R23 | Tag |
|---|---|---|---|---|---|
| 2.1 | 469 | `!is.null(seed)` | parallel path: seed from `sample.int` (SD:1036) so non-NULL | same | BOTH. Sequential path (`num_simulations < 100`, SD:968-974, or `parallel_replicates = FALSE`, SD:1054) passes `seed = NULL`: the chunk then continues the MASTER stream (see 3.4). |
| 2.2 | 482 | `collect_raw_n` | TRUE | TRUE | BOTH |
| 2.3 | 491, 498 | `isTRUE(return_variance)` | TRUE | TRUE | BOTH |
| 2.4 | 529, 572 | `rebalance_done <- is.null(rebalance_threshold)`; `!rebalance_done && !is.null(...)` | threshold 50, so checked each arrival until done | same | BOTH. NULL threshold not exercised. |
| 2.5 | 552 | `looks_stop_at_final && current_time >= max_trial_time` | TRUE branch active | active | BOTH |
| 2.6 | 561 | `looks_stop_at_final && isTRUE(compare_arms_option)` | FALSE (single-arm) | TRUE: phantom-arrival break 561-570 active | R23 only |
| 2.7 | 605 | `!is.null(pt_targets_abs)` | FALSE, so else branch 628-633 (calendar beats) | same | NOT EXERCISED: 605-627 (person-time milestones, `cum_person_time_all_arms`, `.__backstop_fired__`) |
| 2.8 | 629 | `current_time >= next_calendar_look` | exercised | exercised | BOTH |
| 2.9 | 635, 638 | `all(arm_status != "recruiting")`; `length(elig_idx) == 0` | R01: arm can stop, break fires | R23: reference never stops (see 6.9), `all(...)` never TRUE | R01 mainly |
| 2.10 | 640 | `isTRUE(compare_arms_option)` | FALSE | TRUE (645 break when no experimental arm can enroll) | split |
| 2.11 | 658 | `!random_dropout` | TRUE: `t_random_censor <- Inf` AFTER the uniform at 657 is drawn | same | BOTH. The `random_dropout = TRUE` stored value (the uniform) is NOT EXERCISED. |
| 2.12 | 665 | `idx > registry_capacity` | stop() | stop() | error path only |
| 2.13 | 681 | `!looks_stop_at_final` | FALSE, so the else block 697-716 runs | same | NOT EXERCISED: 681-696 (look at end of accrual plus beats to `max_trial_time`) |
| 2.14 | 701 | `is.finite(max_trial_time)` | TRUE (72) | TRUE | `last_enroll_time + F` with infinite `max_trial_time` NOT EXERCISED (704) |
| 2.15 | 710 | `is.finite(t_star) && any(arm_status == "recruiting")` | exercised | exercised | BOTH |
| 2.16 | 713 | `post_enroll_time <= current_time` then `next` | exercised when an arrival-driven look was the last event | rare, see 6.12 | BOTH |
| 2.17 | 723 | `all_stopped_early && any(!is.na(stop_time))` | reachable (all arms stopped) | never (reference stays "recruiting") | R01 |
| 2.18 | 726 | `is.finite(max_trial_time)` (final time) | TRUE | TRUE | the Inf branch (729) NOT EXERCISED |
| 2.19 | 732 | `!identical(interval_cutpoints_current, interval_cutpoints_sim)` | TRUE once rebalanced | same | BOTH. Redundant: line 591 already set `interval_lengths`. |
| 2.20 | 737 | `isTRUE(compare_arms_option)` (reference slice) | FALSE | TRUE | split |
| 2.21 | 751 | `state$arm_status[arm] != "recruiting"` skip | exercised | exercised (reference counts as recruiting) | BOTH |
| 2.22 | 758 | `!isTRUE(compare_arms_option)` (HC final vs vs-ref final) | HC final 759-793 | vs-ref final 794-833 | split |
| 2.23 | 770 | `!is.null(median_pfs_success_threshold_arms) && arm %in% names(...)` | TRUE | n/a | R01. The `Inf` fallback (773) NOT EXERCISED. |
| 2.24 | 775 | `is.finite(success_thr)` | TRUE | n/a | R01. The `0` fallback NOT EXERCISED. |
| 2.25 | 777 | `!is.na(p_eff) && p_eff >= final_success_posterior_prob_threshold` | exercised | n/a | R01 |
| 2.26 | 781, 786 | futility analogs (`-Inf` fallback, `0` fallback) | TRUE branches | n/a | R01. Fallbacks NOT EXERCISED. |
| 2.27 | 788 | `!is.na(p_fut) && p_fut >= final_futility_posterior_prob_threshold` | exercised | n/a | R01 |
| 2.28 | 795 | `arm == reference_arm_name` (reference arm marked inconclusive) | n/a | exercised | R23 |
| 2.29 | 804 | `isTRUE(args_local$use_ph_model_vs_ref) && !is.null(med_samples$logHR)` | n/a | TRUE | R23. Non-PH branch 818-823 (uses `final_vsref_probs_abs`) NOT EXERCISED. |
| 2.30 | 806 | `!is.null(args_local$compare_arms_hr_margin)` | n/a | TRUE (0.0337897681037559) | R23. NULL (log_margin 0) NOT EXERCISED. |
| 2.31 | 813 | `!is.null(args_local$compare_arms_hr_eff_margin)` | n/a | FALSE (NULL), so `log_eff <- 0` | R23. Non-NULL NOT EXERCISED. |
| 2.32 | 825, 827 | `p_eff_ref >= efficacy_threshold_vs_ref_prob`, `else if p_fut_ref >= futility_threshold_vs_ref_prob` | n/a | exercised. A NULL threshold would make `if` error ("argument is of length zero"). | R23 |
| 2.33 | 845 | `arm %in% names(cached_arm_slices)` | recruiting arms cached, stopped arms re-sliced | same | BOTH |
| 2.34 | 866, 874 | `!is.null(chunk_final_n_raw)`, `!is.null(chunk_final_events_raw)` | non-NULL | non-NULL | BOTH. NULL variants NOT EXERCISED. |

### 2B. `interim_check` and helpers (IL, SM, PH, PD)

| # | Line | Condition | R01 | R23 | Tag |
|---|---|---|---|---|---|
| 2.35 | IL:78 | `isTRUE(args$compare_arms_option)` dispatch | HC | vs-ref | split |
| 2.36 | IL:86 | reference arm missing: stop | n/a | not hit | error path |
| 2.37 | IL:91, 100 | no experimental arms / none recruiting: return state | n/a | 100 reachable | R23 |
| 2.38 | IL:106-109 | `compare_arms_target_arm` override | n/a | field absent, `target_override` NULL, `eval_arms <- active_experimental` | dead |
| 2.39 | IL:117 | `!use_ph && length(eval_arms) > 1` (precompute control posteriors) | n/a | FALSE | NOT EXERCISED. If entered it draws gamma RNG whose result the C++ dispatcher ignores (PD:169-191 never touches `ctrl_cache`). Stream-shifting if ever enabled. |
| 2.40 | IL:122 | `state$arm_status[trt_name] != "recruiting"` next | exercised | exercised | R23 |
| 2.41 | IL:127 | `gates_pass_for_both_arms(...)` | n/a | exercised | R23 |
| 2.42 | IL:10, 33 | PH branch of `calculate_current_probs_vs_ref` | n/a | PH (`logHR` present) | non-PH branch 33-39 NOT EXERCISED |
| 2.43 | IL:13, 26 | `compare_arms_hr_margin` / `compare_arms_hr_eff_margin` non-NULL | n/a | margin set, eff margin NULL | as 2.30, 2.31 |
| 2.44 | IL:147-149 | `isTRUE(efficacy_stopping_rule_vs_ref) && !is.null(thr) && is.finite(thr) && pr_eff >= thr` | n/a | rule FALSE, so efficacy never stops at an interim | R23. TRUE case NOT EXERCISED. |
| 2.45 | IL:158-159 | `!is.null(futility_threshold_vs_ref_prob) && is.finite(..) && pr_fut >= thr` | n/a | exercised. **The flag `futility_stopping_rule_vs_ref` (TRUE in R23) is never consulted.** | R23 |
| 2.46 | IL:173 | `arm_status != "recruiting"` next | exercised | n/a | R01 |
| 2.47 | IL:179 | `n_pat > 0` for median | exercised | n/a | R01 |
| 2.48 | IL:182-184 | `coalesce_num` of gates | min_patients 0, min_events_hc 3, min_median_fu_hc 0 | n/a. Note R23 min_patients_for_analysis = 10 and min_events_hc = 11 are IGNORED on the vs-ref path. | R01 |
| 2.49 | IL:187-195 | `resolve_gate_vec(min_person_time_frac_per_arm, scale = TRUE)` | raw 0, single arm: branch "length(raw) == length(all_arm_names)" (GD:31-35); scaling by `probs/max_prob` = 1 | n/a | R01. Value 0 so gate never binds. |
| 2.50 | IL:203-215 | `mtpa` named / positional / first-element fallbacks | named branch (204) | n/a | R01. Other two NOT EXERCISED. |
| 2.51 | IL:216 | `maxPT_arm > 0` | TRUE | n/a | R01 |
| 2.52 | IL:218-223 | gate: `n_pat < min_pat \|\| events_total < min_events \|\| median_fu < min_median_fu \|\| pt_frac < min_pt_frac` | exercised | n/a | R01 |
| 2.53 | IL:240-243 | `isTRUE(efficacy_stopping_rule_hc) && ...` | FALSE (futility_only), no interim efficacy stop | n/a | R01. TRUE NOT EXERCISED. |
| 2.54 | IL:251-253 | futility stop. **`futility_stopping_rule_hc` (TRUE) is never consulted.** | exercised | n/a | R01 |
| 2.55 | IL:57-61 | `success_thr_vec <- median_pfs_success_threshold_arms %\|\|% null_median_arms`; `futility_thr_vec <- median_pfs_futility_threshold_arms %\|\|% futility_median_arms %\|\|% success_thr_vec`; `[[arm_name]] %\|\|% [1]` | first operand non-NULL, name present. Note: `[[` with an absent name on an atomic vector ERRORS, so the `%\|\|% [1]` fallback never helps. | n/a | R01 |
| 2.56 | IL:45-46 | `num_posterior_draws_interim %\|\|% num_posterior_draws`, `max(1L, as.integer(.))` | `num_posterior_draws_interim` is `num_draws_interim` (SD:370-374) = 5000 | n/a | R01 |
| 2.57 | SM:154, 159, 164-178 | argument validation; `trt_name` supplied | not hit | passes | R23 |
| 2.58 | SM:226, 228, 231 | `coalesce_num(max_trial_time, ...)`; null/names checks | `max_trial_time` 72 | 72 | R23 |
| 2.59 | SM:238-239 | `maxPT > 0 ? pt/maxPT : 0` | n/a | exercised | R23 |
| 2.60 | GD:23-77 | `resolve_gate_vec` branches | see 2.49 | raw scalar, two target arms: "recycle first element" branch (GD:40-41), then proportional scaling (GD:44-68). `min_events_per_arm` 11L -> 11*1, `min_median_followup_per_arm` 3 (scale = FALSE), `min_person_time_frac_per_arm` 0 | split. Named-raw and positional-raw branches NOT EXERCISED. |
| 2.61 | PH:96 / PH:93-103 | `isTRUE(args$use_ph_model_vs_ref)` | n/a | PH | independent model (PD:169-196) NOT EXERCISED |
| 2.62 | PH:95 | `num_samples %\|\|% args$num_posterior_draws` | passed explicitly | passed explicitly | (800) |
| 2.63 | PH:113-114 | `coalesce_num(ph_loghr_prior_mean, 0)`, `max(1e-6, coalesce_num(ph_loghr_prior_sd, 1))` | n/a | 0 and 1 | R23 |
| 2.64 | SM:14, 25, 33-37 | empty patient data; `any(event_status == 1L)`; `length(idx) > 0` | exercised | exercised | BOTH |
| 2.65 | SM:118, 144 | empty registry early return; `nrow(pd) > 0` for median | exercised early in the trial | exercised | BOTH |
| 2.66 | SD:1332 | `length(event_times) < max(2, num_intervals)` returns NULL | rebalance_done set TRUE (SD:600) | same | BOTH |
| 2.67 | SD:1341, 1345 | non-increasing quantiles bump by 1e-6; `any(diff(q_vals) <= 0)` returns NULL | exercised rarely | same | BOTH |

### 2C. Hidden configuration reads (outside the argument list)

| # | Where | What | Frozen value | Note |
|---|---|---|---|---|
| 2.68 | SD:236 | `getOption("evolveTrial.looks_stop_at_final", FALSE)` | TRUE | Evaluated lazily when the promise is first forced, in the master. Read once per `run_simulation_pure` call, not per trial. |
| 2.69 | PD:12 | `Sys.getenv("EVOLVETRIAL_USE_CPP", "TRUE")`, read on EVERY call of `draw_posterior_hazard_samples` (PD:37), `ph_beta_mode_var` (PD:100) | unset, so TRUE (C++ branches) | Values outside TRUE/1/YES switch to the R fallbacks (`draw_posterior_hazard_samples_r`: R `rgamma` per column; `ph_beta_mode_var_r`: R `sum`). They are bit-different from the C++ ones for `ph_beta_mode_var` (long double vs Armadillo sums, see 4.30). The port must read this once and REFUSE (fall back to the R closure) when it is not truthy. |
| 2.70 | SD:965 | `getOption("evolveTrial.sequential_threshold", 100L)` | 100 | outside the chunk, decides sequential versus parallel (and thus `seed`) |
| 2.71 | SD:979-982 | `utils::packageName()` | irrelevant | outside the chunk |

### 2D. Arguments the frozen base_args carry that never reach the chunk

`run_scenarios` keeps only formals of `run_simulation_pure` (SD:1274-1287). Dropped from both args files:
`min_patients_per_arm`, `enable_efficacy_stopping`, `enable_futility_stopping`, `efficacy_threshold`,
`futility_threshold`, `modeling_horizon`. `cohort_size_per_arm` (1L) is a formal but is never read by any
code (SD:34, 202, 1236 only, 6.5). `progress` defaults to `interactive()` (SD:263), FALSE under Rscript.

Count: 71 numbered items (2.1 to 2.71) in 2A to 2C, plus 2D (6 dropped args, 1 unused formal).

---------------------------------------------------------------------------------------------------

## 3. Side effects other than the return value

| # | Where | Effect | Frozen? |
|---|---|---|---|
| 3.1 | SD:882 `tick()` | once per trial, at the end. Parallel path passes `function() {}` (SD:1043). Sequential path passes `tick_fun` (SD:442, 1054, 974), which is `pb$tick()` only if `progress` is TRUE and not parallel (SD:350), otherwise `invisible(NULL)`. | no-op in frozen runs |
| 3.2 | SD:470 `set.seed(seed)` | resets the process RNG (Mersenne-Twister, Inversion, Rejection defaults) from an integer in `[1, .Machine$integer.max]`. Worker-local in the parallel path. Also materialises `.Random.seed` in the worker's global env. | yes |
| 3.3 | Every RNG call | R-level `rexp`, `sample`, `rweibull`, `runif`, `rnorm`, `rgamma` and the C++ `R::rgamma` (inside Rcpp's `RNGScope` of the generated wrapper) read and write `.Random.seed` (GetRNGstate/PutRNGstate). A C++ chunk must take the state once at entry and write it back at exit to leave `.Random.seed` exactly as R would. | yes |
| 3.4 | Sequential path (seed NULL) | the chunk consumes the CALLER's stream. `run_scenarios` calls `set.seed(seed)` once (SD:1272) and then runs scenarios in order on that stream (SD:1305, `lapply`). In the sequential path scenario 2's draws depend on exactly how many uniforms scenario 1 consumed, and anything after the call (e.g. BATON) sees the advanced stream. In the parallel path the master consumes only one `sample.int` per chunk (SD:1036) and the worker's RNG is discarded. | parallel in frozen runs (num_simulations >= 100); sequential path NOT EXERCISED but must be bit-exact if allowed |
| 3.5 | `message()` when `diagnostics` is TRUE | SD:594-597 (rebalance), IL:92, 101, 130, 140-143, 230-235, SM:247-251. Frozen `diagnostics = FALSE`. | NOT EXERCISED |
| 3.6 | `warning()` / `Rcpp::warning` | PS:54 "Posterior beta parameter ... non-positive. Setting to 1e-6" (draw_posterior_hazard_samples_cpp); PD:74 R analog; R's `rgamma` / `rnorm` "NAs produced" warnings if a parameter is NaN. None expected (prior rate > 0, person-time >= 0). | not triggered |
| 3.7 | `stop()` | SD:665-670 registry overflow; IL:86-87 missing reference arm; SM:155-156, 160-161, 167, 173, 176, 229-233; PS:40 length mismatch; `FixupProb` errors in `sample()` ("NA in probability vector", "too few positive probabilities"); `if (NA)` errors when a probability is NaN; `if (logical(0))` at SD:825 for a NULL threshold. | error paths |
| 3.8 | `<<-` assignments | **None inside `simulate_chunk`.** The only `<<-` are in `aggregate_results` (SD:906-941), which runs on the master after the chunks return. | n/a |
| 3.9 | `state$.__backstop_fired__` (SD:608) | adds a list element to `state`; local to the trial. | NOT EXERCISED |
| 3.10 | Environment writes | none besides `.Random.seed`. No files, no options changes, no `Sys.setenv`. | n/a |
| 3.11 | Attributes | names only. `setNames` on accumulators (472-478), `names(tmp) <- arm_names` (484, 494). No class, dim or other attributes on the returned objects. | n/a |
| 3.12 | `parallel::parLapply` serialisation | the whole frame, including `pb` if created, is serialised to the worker; no semantic effect on results. | n/a |

Count: 12 items.

---------------------------------------------------------------------------------------------------

## 4. Floating-point operations that decide a branch or feed an output

Legend: LD = R accumulates in `long double` (80-bit x87 on x86-64) then casts to double. "int" means
exact integer arithmetic. All frozen runs are on x86-64 Linux with R's default build flags
(`R CMD config CXXFLAGS`: `-g -O2 ...`, no `-march`, no `-ffast-math`, so no FMA contraction; the
CPU has FMA, so a port built with `-march=native` or `-mfma` WOULD change results).

### 4A. RNG-fed values and timeline (SD)

| # | Line | Operation | Notes |
|---|---|---|---|
| 4.1 | 329-331, DG:50-55 | `weibull_scale = median / (log(2)^(1/shape))` | libm `log`, R's `^` (`R_pow`: `x == 1 \|\| y == 0` shortcuts, `y == 2` uses `x*x`, else libm `pow`). Compute once in R and pass the scale vector in. |
| 4.2 | 543-544 | `current_time <- current_time + rexp(1, rate)` | R's `rexp(n, rate)` evaluates `1/rate` in R then C `scale * exp_rand()`. In C++: `R::rexp(1.0/rate)`, NOT `exp_rand()/rate`. Accumulation by repeated addition. |
| 4.3 | 552 | `current_time >= max_trial_time` | exact comparison, boundary excluded from enrolling |
| 4.4 | 650-651 | `probs / sum(probs)` | `sum` of doubles is LD (checked: sum(dbl) equals long-double accumulation in 50,000 of 50,000 random cases and differs from plain double accumulation in 28,433 of 50,000). Then a double divide. |
| 4.5 | 652 | `sample(elig_arms, size = 1, prob = probs)` | `sample.int(n, 1, FALSE, prob)` -> `.Internal(sample())`: `FixupProb` divides again by a plain-double in-order sum of the (already normalised) p; then `revsort` (descending, heap-based, NOT stable); then cumulative sum of the sorted p; then ONE `unif_rand()` and `rU <= p[j]` scan. From R source memory, verify by test. **Checked:** with p = (0.5, 0.5) the draw `u <= 0.5` selects the LAST arm ("b"), i.e. tie order after `revsort` is reversed. **Checked:** `sample(x, 1, prob = 1)` (single eligible arm) still consumes exactly one uniform. Port: call R's `revsort` (R_ext/Utils.h) and replicate the two normalisations. |
| 4.6 | 655-656 | `rweibull(1, shape, scale)` | libR: `scale * pow(-log(unif_rand()), 1./shape)`. One uniform. Returns 0 if `scale == 0`, NaN for non-finite/non-positive parameters. |
| 4.7 | 657 | `runif(1, 0, censor_max_time_sim)` | libR: rejection loop until `0 < u < 1`, then `a + (b-a)*u`. Always drawn, also with `random_dropout = FALSE` (658 only overwrites the stored value with Inf). If `b` is non-finite, R returns NaN and draws nothing (from the earlier inventory, `rng_driver_inventory.md`). |
| 4.8 | 629, 631 | `current_time >= next_calendar_look`; `next_calendar_look <- next_calendar_look + beat` | repeated addition from `interim_calendar_beat` (initial value, 534), one beat per look. If an arrival gap exceeds a beat, `next_calendar_look` lags and each following arrival fires a look until it catches up. |
| 4.9 | 584 | `total_events >= rebalance_threshold` | `total_events` is an integer (0L + integer counts); threshold double |
| 4.10 | 505-519 (`extract_events_fast`) | `pmax(0, t - enroll)`, `pmin(., max_follow_up)`, `> 0`, `pmin(te, rc, ta)`, `te <= pmin(rc, ta)`, `sum(events_mask)` | the mask sum is an integer sum of a logical. `rc` is Inf under `random_dropout = FALSE`. Events returned as `observed_time[mask]`, which equals `te` exactly. |
| 4.11 | 1335 | `pmin(pmax(event_times, 0), max_follow_up)` | exact |
| 4.12 | 1336 | `probs <- seq(0, 1, length.out = num_intervals + 1)` | checked: equals `i * (1/n)` (and `0 + i*by`, `by = 1/n`) bit-for-bit for n in {3,7,8,9,12,13,15,16}, but NOT `i/n` for n = 7, 9, 12, 13. Use `i * (1.0/n)`. |
| 4.13 | 1337 | `quantile(event_times, probs, na.rm = TRUE, type = 7)` | R 4.6.0 `quantile.default`: `probs <- pmax(0, pmin(1, probs))`; `index <- 1 + max(n-1, 0) * probs`; `lo <- floor(index * (1 + fuzz))` with `fuzz = 0` for type 7; `hi <- ceiling(index)`; order statistics from a partial sort; `qs <- x[lo]`; where `index > lo & x[hi] != qs`: `qs[i] <- (1-h)*qs[i] + h*x[hi[i]]`, `h <- (index - lo)[i]`. Values depend only on order statistics, not on the order of `event_times`. |
| 4.14 | 1338-1339 | `q_vals[1] <- 0`, `q_vals[K+1] <- max_follow_up` | after rebalancing the LAST cutpoint becomes `max_follow_up_sim` (120), not the modelling horizon (24), see 6.15. |
| 4.15 | 1340-1344 | `!is.finite(q) \|\| q <= q_prev` then `min(max_follow_up, q_prev + 1e-6)` | sequential, uses the already-bumped previous value |
| 4.16 | 1345 | `any(diff(q_vals) <= 0)` returns NULL | exact |
| 4.17 | 591, 733 | `diff(new_cuts)` | `x[i+1] - x[i]` |
| 4.18 | 732 | `identical(interval_cutpoints_current, interval_cutpoints_sim)` | exact bit identity of double vectors |

### 4B. Slicing and metrics (SM)

| # | Line | Operation | Notes |
|---|---|---|---|
| 4.19 | SM:130-136 | `ts <- pmax(0, cal - enroll)`; `ta <- pmin(ts, max_fu)`; `obs <- pmin(te, rc, ta)`; `status <- as.integer(te <= pmin(rc, ta))`; `keep <- ta > 0` | patients enrolled exactly at `cal` (ta = 0) drop out; ties `te == ta` count as events |
| 4.20 | SM:27-33 | `findInterval(event_times, cutpoints, left.open = FALSE, rightmost.closed = FALSE)`; drop 0s; `pmin(idx, num_intervals)`; `tabulate` | `idx` = number of cutpoints `<= x`. Events at or beyond the LAST cutpoint all land in the last interval. Events before cutpoint 1 (<= 0) cannot occur. |
| 4.21 | SM:45-59 | for each interval: `at_risk <- obs > lo` (strict); `contrib <- pmin(obs[at_risk], hi) - lo`; `contrib[contrib < 0] <- 0`; `person_time[i] <- sum(contrib)` | `sum` is LD. Person-time is TRUNCATED at the last cutpoint while events past the last cutpoint are still counted (6.14). |
| 4.22 | SM:142 | `events_total <- sum(events_per_interval)` | integer sum (empty-registry branch returns `0L`, with double zero vectors at SM:122-123) |
| 4.23 | SM:143 | `person_time_total <- sum(person_time_per_interval)` | LD, over K already-LD-rounded interval sums (two-level summation) |
| 4.24 | SM:144 | `stats::median(pd$observed_time)` | R 4.6.0 `median.default`: odd n -> `sort(x, partial = half)[half]`; even n -> `mean(sort(x, partial = half + 0:1)[half + 0:1])`. `mean` of two doubles is R's two-pass long double mean: `s = (a+b)/2` (LD), `t = (a-s) + (b-s)` (LD), `s += t/2`, cast to double. Checked: equals plain `(a+b)/2` in 20,000 of 20,000 random pairs, but implement the LD formula. |
| 4.25 | IL:178-180 | `sum(patient_data$event_status)` (integer), `median` as 4.24, `sum(person_time_per_interval)` (LD) | recomputed in HC; same values as the slice metrics |
| 4.26 | IL:201-216 | `maxPT_arm <- coalesce_num(mtpa[[arm]], 0) * max_follow`; `pt_frac <- pt_total / maxPT_arm`; `max_follow = coalesce_num(max_trial_time, coalesce_num(max_follow_up_sim, 0))` | note `max_trial_time` (72) not `max_follow_up_sim` (120) is the denominator; Inf `max_trial_time` gives frac 0 |
| 4.27 | IL:218-221 | `n_pat < min_pat \|\| events_total < min_events \|\| median_fu < min_median_fu \|\| pt_frac < min_pt_frac` | integer-vs-double compares; gate fails when strictly below |
| 4.28 | SM:235-242, GD:44-68 | `maxPT_C <- mtn[[ctrl]] * max_follow`; `fracC <- ptC / maxPT_C`; `passC <- evC >= min_ev && mfuC >= min_mfu && fracC >= min_pt` | `min_ev_vec` = raw value `* (probs/max_prob)` (GD:44-68), e.g. 11 * (0.5/0.5) = 11 exactly for R23; for general probs this is a double multiply and `evC >= 11*x` compares an integer count with a double. |

### 4C. Posterior draws and medians

| # | Line | Operation | Notes |
|---|---|---|---|
| 4.29 | IL:47-54, PD:39-46, PS:48-65 | `draw_posterior_hazard_samples_cpp`: for each interval j (outer), for each draw i (inner): `R::rgamma(prior_alpha[j] + events[j], 1.0/(prior_beta[j] + person_time[j]))` | column-major draw order, `posterior_alpha_j = alpha + (double)events`, scale computed once per column. Same draws as the R fallback `rgamma(n, shape, rate)` (R forms `scale = 1/rate` first), which is why the stream is stable across both. |
| 4.30 | PH:116 -> PD:98 -> PS:196-261 | `ph_beta_mode_var_cpp` (Newton solver, live) | uses **Armadillo `sum`** (`sum(E_C)`, `sum(PT_C)`, `sum((alpha + conv_to<vec>::from(E_C+E_T)) % g)`, etc.). Armadillo's accumulators are plain double with a two-accumulator unrolled loop, NOT R's long double `sum`. The R fallback `ph_beta_mode_var_r` (PD:123-155) uses LD sums and would differ. The port must call the same C++ function (or replicate Armadillo's exact summation order), not re-implement it with a naive loop. Iteration stops at `abs(step) < 1e-6` or 50 iterations, or when `beta_new` is non-finite; variance `min(1e6, max(1e-6, -1/hess))` if `hess < 0` else 1. |
| 4.31 | PH:124 | `rnorm(num_samples, mean = beta_params$mean, sd = beta_params$sd)` | R `rnorm` = `mu + sigma * norm_rand()` (inversion; two uniforms per normal: `BIG = 134217728`). Returns `mu` WITHOUT drawing if `sigma == 0` or `mu` non-finite. `sd >= 1e-3` here. All `num_samples` normals are drawn BEFORE any gamma. |
| 4.32 | PH:125 | `shapes <- alpha_prior + E_C + E_T` | `(alpha + E_C) + E_T`, exact for the integers involved |
| 4.33 | PH:127-133 | for each draw i: `exp_beta <- exp(beta_draws[i])`; `rates <- beta_prior + PT_C + exp_beta * PT_T`; `lambda <- rgamma(K, shape = shapes, rate = rates)`; medians for `lambda` and `lambda * exp_beta` | evaluation order `(beta_prior + PT_C) + (exp_beta * PT_T)`; R forms `scale = 1/rates` then draws K gammas in interval order. **`med_ctrl`/`med_trt` are never used on the frozen PH path** (final: SD:804-817 uses only `logHR`; interim: IL:29-32 only `logHR`), but the gamma draws must still be made to keep the stream aligned (6.1). |
| 4.34 | PH:16-57, PS:98-135 | piecewise median: `cum_end = cum + lambda*dt`; `if cum_end >= log(2)`: `t + (log(2) - cum)/lambda`; else `cum <- cum_end`, `t <- t + dt`; tail `t + (log(2) - cum)/lambda_last` if `lambda_last > 0` else Inf | R uses `log(2)` (libm), C++ uses the constant 0.69314718055994530942; same double. No FMA permitted (`lambda*dt` added to `cum`). |
| 4.35 | IL:56, SD:768 | `calculate_median_survival_matrix_cpp(lam, interval_lengths)` | row-wise `_impl` as 4.34 |
| 4.36 | IL:63-64, SD:775, 786 | `mean(med_draws > thr)`, `mean(med_draws < thr)` | **`mean()` of a LOGICAL vector in R is count divided by n in long double, then cast** (checked: matches `(double)((long double)c / (long double)n)` in 100,000 of 100,000 random (c, n), but differs from plain `c/n` in 12 of 100,000). Compare `NA` cases: NaN medians make `mean` NA. |
| 4.37 | IL:30-31, SD:816-817 | `mean(logHR < log_eff)`, `mean(logHR >= log_margin)` | same as 4.36 |
| 4.38 | IL:18, 27, SD:808, 814 | `log1p(-min(max(0, margin), 0.95))` | libm `log1p` |
| 4.39 | IL:149, 159, SD:777, 788, 825, 827 | `pr >= threshold` | exact double comparison; the mean from 4.36 may sit within 1 ulp of a threshold only by coincidence |
| 4.40 | SD:775 | `p_eff <- if (is.finite(success_thr)) mean(...) else 0`, then `!is.na(p_eff) && p_eff >= thr` | |

### 4D. Final time and accumulation

| # | Line | Operation | Notes |
|---|---|---|---|
| 4.41 | 701-715 | `t_star <- min(last_enroll_time + F, max_trial_time)`; `post_enroll_time <- last_enroll_time`; loop: `post_enroll_time <- post_enroll_time + beat`; `>= t_star` break; `<= current_time` next | **repeated addition from `last_enroll_time`**, NOT `last_enroll_time + j*beat` (the comment at 698 says j*beat). Different rounding after many beats. |
| 4.42 | 723-730 | `max(state$stop_time, na.rm = TRUE)` or `min(last_enroll_time + F, max_trial_time)` | exact |
| 4.43 | 856-863 | `chunk_sum_* <- chunk_sum_* + vec`; `chunk_sum_final_time <- chunk_sum_final_time + final_time` | per-trial sequential double additions in trial order, starting from 0. The chunk partial sum is later added to the global accumulator (SD:906-918) and scaled by `1/num_simulations` (SD:1060-1081, multiplication by the reciprocal, not division). The grouping depends on the chunk partition (`splitIndices`). |
| 4.44 | 336 | `True_Median = round(weibull_median_true_arms, 2)` | outside the chunk |
| 4.45 | 359-361, 367 | `max_PT_per_arm`, `total_max_PT <- sum(...)` (LD), `pt_targets_abs` | computed always, used only with milestones (NOT EXERCISED) |
| 4.46 | DG:37 | `cum_person_time_arm`: `sum(observed_time)` (LD) | NOT EXERCISED |
| 4.47 | SD:1086-1101 | `stats::var(raw) / num_simulations` (R `var` is a two-pass LD algorithm) | outside the chunk but consumes the chunk's `final_*_raw` vectors, so their exact double values matter |

Count: 47 items.

---------------------------------------------------------------------------------------------------

## 5. Return value structure (SD:885-898)

An unnamed-call `list(...)` of exactly 12 elements in this order:

| # | Name | Type, names, length |
|---|---|---|
| 5.1 | `sum_final_n` | double vector, names `arm_names`, length n_arms. Built as `numeric + final_n_vec` (472, 856). |
| 5.2 | `sum_final_events` | double vector, names `arm_names` |
| 5.3 | `sum_stop_eff` | double vector, names `arm_names` (double + integer `stop_efficacy_per_sim_row`) |
| 5.4 | `sum_stop_fut` | double vector, names `arm_names` |
| 5.5 | `sum_final_eff` | double vector, names `arm_names` (double + integer `final_eff_vec`) |
| 5.6 | `sum_final_fut` | double vector, names `arm_names` |
| 5.7 | `sum_final_inc` | double vector, names `arm_names` |
| 5.8 | `sum_final_time` | double scalar, NO names (0 + unnamed `final_time`) |
| 5.9 | `n_sims` | **integer** scalar (`length(sim_indices)`) |
| 5.10 | `final_n_raw` | NULL, or a list of length n_arms with names `arm_names`; each element a double vector of length `length(sim_indices)`, NO names, zero-initialised, filled by position `which(sim_indices == sim_idx)`. Frozen: non-NULL. |
| 5.11 | `final_events_raw` | NULL, or same shape as 5.10. Frozen: non-NULL (return_variance TRUE). |
| 5.12 | `final_time_raw` | NULL, or unnamed double vector of length `length(sim_indices)`. Frozen: non-NULL. |

Notes:
- `list(final_n_raw = NULL, ...)` keeps a NULL ELEMENT in place, so the list always has length 12.
- Names on the sum vectors come from the first operand of `+` (the `setNames(numeric(...), arm_names)` accumulator).
  With a single arm they are length-1 named vectors, not scalars.
- `aggregate_results` (SD:901-944) depends on structure: it tests `!is.null(res$final_n_raw) && !is.null(all_final_n_raw)`, advances `raw_offset` only inside the `final_n_raw` block (SD:922-929), and the `final_events_raw` / `final_time_raw` blocks (SD:931-942) index backwards from that offset. A chunk that returns `final_events_raw` without `final_n_raw` would corrupt the placement. The port must return all three together or none.
- Elements 5.1 to 5.8 enter `sum_*` accumulators and are later multiplied by `inv_num_sims` (SD:1060-1081).

Count: 12 elements.

---------------------------------------------------------------------------------------------------

## 6. Surprising or risky items

### 6A. RNG alignment and hidden dependencies

- 6.1 **Dead computation that still consumes RNG.** On the frozen PH path the `rgamma` draws at PH:130 feed only `med_ctrl` / `med_trt`, which are never used for the decision (4.33). The port may skip the median computation but MUST make the 800 x K gamma draws, after the 800 normals, in that order. Skipping them shifts the whole stream.
- 6.2 **Two definitions of `sample_vs_ref_medians_ph`.** PD:241-271 (dispatches to `sample_vs_ref_medians_ph_cpp` (PS:285), which uses Armadillo's `randn`, a different normal generator) and PH:105-135 (pure R loop). DESCRIPTION has no `Collate:` field, so files load alphabetically and `posterior_helpers.R` wins. Confirmed on the installed namespace: `deparse(evolveTrial:::sample_vs_ref_medians_ph)` shows the `coalesce_num` / R-loop version. The C++ one (PS:285 onward) is dead and does NOT match. Also `sample_vs_ref_medians_ph_r` (PD:276) is a dead third copy.
- 6.3 **`sample()` with a single eligible arm still draws one uniform** (checked), and with tied probabilities the tie order is reversed by `revsort` (checked: `u <= 0.5` picks the last of two equally likely arms). R23 randomises 0.5/0.5, so this decides every patient assignment.
- 6.4 **Per-patient draw order**: `rexp` (every arrival, including the one that breaks the loop) -> optional rebalance/look (posterior draws) -> `sample` -> `rweibull` -> `runif`. Draws made inside an interim look occur between the arrival `rexp` and the patient's `sample`. The `rexp` of an arrival that then triggers a break (635, 638, 645; and 552/566 before the break tests in P6 mode) is consumed with no patient enrolled.
- 6.5 **Unused formals/args**: `cohort_size_per_arm` is never read; `predictive_fast`, `num_posterior_draws_pred`, `pred_*`, `predictive_boundary_band`, `run_pred_every_k_beats` are packed (SD:408-413) and never read.
- 6.6 **Flags that look like switches but are not read**: `futility_stopping_rule_hc` (SD:212, 424), `futility_stopping_rule_vs_ref` (SD:214, 427). Futility stopping is governed only by the threshold being non-NULL and finite (IL:158-159, 251-253). `efficacy_stopping_rule_*` ARE read (IL:147, 240) but only at interims, not at the final.
- 6.7 **Gate asymmetry**: HC uses `min_patients_for_analysis`, `min_events_hc`, `min_median_followup_hc`, `min_person_time_frac_per_arm`; vs-ref uses `min_events_per_arm`, `min_median_followup_per_arm`, `min_person_time_frac_per_arm` only. R23's `min_patients_for_analysis = 10` and `min_events_hc = 11` are ignored; R01's `min_events_per_arm = 3L` is ignored (it is also equal to `min_events_hc` = 3).
- 6.8 **Draw counts**: vs-ref interim and final both use `args$num_posterior_draws` (IL:7, SD:802), never `num_posterior_draws_interim`; the HC interim uses `num_draws_interim` (IL:45) but the HC final uses `num_posterior_draws` (SD:765). In the frozen designs `num_draws_interim == num_posterior_draws`, so this does not show.

### 6B. Control flow subtleties

- 6.9 In vs-ref mode the reference arm's status never leaves "recruiting" (only experimental arms are updated, IL:150, 160). Consequently `all(state$arm_status != "recruiting")` (624, 635, 722) is never TRUE, `all_stopped_early` is always FALSE (so `final_time = min(last_enroll_time + F, max_trial_time)` for R23), and accrual ends by the explicit "no experimental arm can enroll" breaks (566, 645). The reference arm counts as `final_inc` at the final (796).
- 6.10 `is_eligible` (538-540) compares `enrolled_counts < max_total_patients_per_arm` POSITIONALLY, while other code indexes `max_total_patients_per_arm` by name (565, 644). Silent wrongness if the two vectors are in different arm orders. `which()` returns a named integer; `arm_names[elig_idx]` uses positions.
- 6.11 `is_eligible` is a closure created inside the per-trial loop (538). `extract_events_fast` is created once per chunk (501). Neither carries state.
- 6.12 In P6 mode, `current_time` after the accrual loop is either `last_enroll_time` (breaks at 553-555, 567-568 reset it) or the time of a later non-enrolling arrival (breaks at 635, 638, 645 do not reset it). Grid points `<= current_time` are skipped (713) so look times stay increasing. The phantom-arrival reset at 561-570 exists to keep this clean for vs-ref.
- 6.13 Rebalancing happens only inside the accrual loop, before the look of the same iteration, using the new arrival's time; the arrival's patient is not yet enrolled. `rebalance_done` stays FALSE for the whole trial if fewer than 50 events ever occur, and then `extract_events_fast` re-slices every arm on EVERY arrival (an O(N^2) hot spot, but also the reason the port must track events exactly the same way). Once `total_events >= 50` the attempt is made once, and `rebalance_done` becomes TRUE whether or not new cutpoints were produced (SD:588-602).
- 6.14 **Person-time/event asymmetry**: person-time per interval is clipped at the interval's upper cutpoint (SM:54), so total person-time stops at the last cutpoint (24 before rebalancing), yet events observed after the last cutpoint are counted in the LAST interval (SM:33). The gate `person_time_total` and the posterior use these truncated values.
- 6.15 After rebalancing the cutpoints span `[0, max_follow_up_sim]` = [0, 120] (SD:1338-1339), not the modelling horizon. The count of intervals K is unchanged (`num_intervals` from SD:430 is used, and priors keep length K).
- 6.16 Per-trial reset: `args_local`, `interval_cutpoints_current`, `interval_lengths`, `rebalance_threshold`, `rebalance_done`, `last_enroll_time`, `state`, `current_time`, `next_calendar_look`, `next_pt_idx`, `patient_id`, `final_*_vec`, `cached_arm_slices` are all re-initialised each trial (523-536, 744-748). **Nothing carries across trials inside a chunk except the RNG stream and the `chunk_*` accumulators/raw vectors (472-498).**
- 6.17 Registries are pre-allocated at `ceiling(1.2 * max_total_patients_per_arm)` rows (SM:74) with `id = NA`; "active" rows are the non-NA-id rows (SM:108-112), i.e. the first `registry_row_idx - 1` rows. `id` is `patient_id`, global across arms in a trial. Only the `NA` pattern is used.
- 6.18 The slice excludes patients with zero follow-up (`keep`, SM:136); `n_pat` for the gate is the count of kept patients, not enrolled patients.
- 6.19 Redundant code: SD:732-734 (line 591 already set `interval_lengths`); `cached_arm_slices` for the events; `elig_idx` computed at 637 and `elig_arms` at 648.

### 6C. Capture and environment discrepancies to be aware of

- 6.20 `args_R23.txt` shows `theta$n_intervals = 15L` but `interval_cutpoints_sim = c(0, 2, ..., 24)` (K = 12) and priors of length 12, while `acceptance_set.md` states K = 15 (`seq(0, 24, length.out = 16)`). The port and its tests must not hard-code K. Worth a check of which one the production wrapper builds.
- 6.21 The captured `num_workers = 1L`, `num_simulations = 300L`; production uses 19 workers x 10,000 trials (`acceptance_set.md:~77`). The chunk partition (`splitIndices`) and per-chunk seeds come from `sample.int` on the master stream (SD:1036) and decide the result; the chunk body itself is partition-agnostic.
- 6.22 `random_dropout = FALSE` is injected by the `verify_m4` wrapper and `min_follow_up_at_final = 6` by the P6 wrapper (`acceptance_set.md:~62-63`); the args files record the post-wrapper values ("after P6 F=6, random_dropout=FALSE").
- 6.23 Platform: long double is 80-bit on x86-64. On aarch64 (128-bit long double) R's own sums, means and the logical `mean` would differ, so bit-identity claims are x86-64 specific. The port must use `long double` for these reductions and must be compiled without `-ffast-math`, `-Ofast`, `-march=native`/`-mfma` (FMA contraction), and with `-ffp-contract=off` to be safe.
- 6.24 `ph_beta_mode_var_cpp` is the one live numerical routine that uses Armadillo reductions instead of R's long double ones (4.30). Calling it unchanged from the port is the safe choice, and it is already exported C++.

Count: 24 items.

---------------------------------------------------------------------------------------------------

## 7. Item counts

| Section | Items |
|---|---|
| 1 Closure free variables | 33 table rows (1.1-1.33) + field-read groups (3 read, 1 dead) + 12 package functions |
| 2 Branches and options | 71 numbered rows (2.1-2.71) + 6 dropped args + 1 unused formal |
| 3 Side effects | 12 |
| 4 Floating-point | 47 |
| 5 Return value | 12 elements |
| 6 Surprises | 24 |
