# P3-2: acceptance frames for gates B, S and the G3 check (plan dev/plans/2026-11-04-p3-2-plan.md, Section 5). Built from
# the R01, R03 and R23 captures with the P3-1 helpers (p31_frames.R: p31_frame_inputs, p31_chunk_frame), same spec format:
# design D, base_args modifications mod (NULL set explicitly), chronology lsaf (TRUE = P6 on, FALSE = P6 off), scen_mod.
# New field same_as: a modification list whose R-body result must equal this frame's bitwise (rule TRUE with a NULL or
# infinite probability threshold must be the rule FALSE frame). Every frame runs with return_variance = TRUE.
source(file.path(normalizePath("~/research/evolveTrial-merge"), "dev/perf/phase3/p31_frames.R"))

p32_frame_specs <- function() {
  effon <- list(efficacy_stopping_rule_hc = TRUE)                       # R01 threshold 0.894, R03 0.974 (captured)
  vson <- list(efficacy_stopping_rule_vs_ref = TRUE)                    # R23 threshold 0.95 (captured), margin 0.0338
  list(
    # ---- single-arm (HC), R01: interim success value 3.8 (median_pfs_success_threshold_arms), futility 4.53
    C01_R01_on_eff        = list(D = "R01", lsaf = TRUE,  mod = effon),
    C02_R01_off_eff       = list(D = "R01", lsaf = FALSE, mod = effon),
    # both tails cross at one look: Pr(med > 3.8) >= 0.3 and Pr(med < 4.53) >= 0.3 hold together when the posterior
    # median sits between the two values; efficacy must win
    C03_R01_on_compete    = list(D = "R01", lsaf = TRUE,  mod = c(effon, list(efficacy_threshold_hc_prob = 0.3,
                                                                                futility_threshold_hc_prob = 0.3))),
    C04_R01_off_compete   = list(D = "R01", lsaf = FALSE, mod = c(effon, list(efficacy_threshold_hc_prob = 0.3,
                                                                                futility_threshold_hc_prob = 0.3))),
    # rule TRUE with a NULL / Inf probability threshold: the R body never stops for efficacy, equal to rule FALSE
    C05_R01_off_effNULL   = list(D = "R01", lsaf = FALSE, mod = c(effon, list(efficacy_threshold_hc_prob = NULL)),
                                 same_as = list(efficacy_stopping_rule_hc = FALSE, efficacy_threshold_hc_prob = NULL)),
    C06_R01_on_effInf     = list(D = "R01", lsaf = TRUE,  mod = c(effon, list(efficacy_threshold_hc_prob = Inf)),
                                 same_as = list(efficacy_stopping_rule_hc = FALSE, efficacy_threshold_hc_prob = Inf)),
    # interim and final success values differ: median_pfs_success_threshold_arms NULL, so the interim reads
    # null_median_arms (scenario value, here 5.0) and the final falls back to Inf (no final efficacy)
    C07_R01_off_thrdiff   = list(D = "R01", lsaf = FALSE, mod = c(effon, list(median_pfs_success_threshold_arms = NULL)),
                                 scen_mod = list(null_median_arms = c(Experimental = 5.0))),
    C08_R01_on_thrdiff    = list(D = "R01", lsaf = TRUE,  mod = c(effon, list(median_pfs_success_threshold_arms = NULL)),
                                 scen_mod = list(null_median_arms = c(Experimental = 5.0))),
    # efficacy on with NULL interim futility probability (13c): the efficacy check still runs
    C09_R01_off_eff_futNULL = list(D = "R01", lsaf = FALSE, mod = c(effon, list(futility_threshold_hc_prob = NULL))),
    C10_R01_on_eff_futNULL  = list(D = "R01", lsaf = TRUE,  mod = c(effon, list(futility_threshold_hc_prob = NULL))),
    # an infinite interim success value: mean(med > Inf) = 0, so no efficacy stop although the rule is on
    C11_R01_off_succInf   = list(D = "R01", lsaf = FALSE, mod = c(effon, list(median_pfs_success_threshold_arms = c(Experimental = Inf)))),
    # ---- single-arm, R03 (efficacy threshold 0.974, min_events_hc 11)
    C12_R03_on_eff        = list(D = "R03", lsaf = TRUE,  mod = effon),
    C13_R03_off_eff       = list(D = "R03", lsaf = FALSE, mod = c(effon, list(efficacy_threshold_hc_prob = 0.8))),
    # ---- between-arm, R23 (PH model, efficacy threshold 0.95, futility 0.372)
    D01_R23_off_current   = list(D = "R23", lsaf = FALSE, mod = vson),
    D02_R23_on_current    = list(D = "R23", lsaf = TRUE,  mod = vson),
    D03_R23_off_symmetric = list(D = "R23", lsaf = FALSE, mod = c(vson, list(compare_arms_hr_margin = 0.1, compare_arms_hr_eff_margin = 0.1))),
    D04_R23_on_symmetric  = list(D = "R23", lsaf = TRUE,  mod = c(vson, list(compare_arms_hr_margin = 0.1, compare_arms_hr_eff_margin = 0.1))),
    D05_R23_off_nomargin  = list(D = "R23", lsaf = FALSE, mod = c(vson, list(compare_arms_hr_margin = NULL, compare_arms_hr_eff_margin = NULL))),
    D06_R23_on_zeromargin = list(D = "R23", lsaf = TRUE,  mod = c(vson, list(compare_arms_hr_margin = 0, compare_arms_hr_eff_margin = 0))),
    D07_R23_off_effonly   = list(D = "R23", lsaf = FALSE, mod = c(vson, list(compare_arms_hr_margin = NULL, compare_arms_hr_eff_margin = 0.1))),
    # 0.97 caps at 0.95 (log 0.05); a negative margin clamps to 0
    D08_R23_on_capped     = list(D = "R23", lsaf = TRUE,  mod = c(vson, list(compare_arms_hr_margin = 0.97, compare_arms_hr_eff_margin = 0.97))),
    D09_R23_off_negative  = list(D = "R23", lsaf = FALSE, mod = c(vson, list(compare_arms_hr_margin = -0.2, compare_arms_hr_eff_margin = -0.2))),
    # both tails cross at one look: with equal margins Pr(eff) + Pr(fut) = 1, so thresholds 0.45 let both hold
    D10_R23_off_compete   = list(D = "R23", lsaf = FALSE, mod = c(vson, list(compare_arms_hr_margin = 0.1, compare_arms_hr_eff_margin = 0.1,
                                                                             efficacy_threshold_vs_ref_prob = 0.45,
                                                                             futility_threshold_vs_ref_prob = 0.45))),
    # V6 without V5: the symmetric margins with interim efficacy off
    D11_R23_off_sym_noeff = list(D = "R23", lsaf = FALSE, mod = list(compare_arms_hr_margin = 0.1, compare_arms_hr_eff_margin = 0.1))
  )
}
