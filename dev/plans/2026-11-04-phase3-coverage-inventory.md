# Phase 3 coverage inventory: why the fast engines do not run everywhere

Worktree `~/research/evolveTrial-perf`, branch perf/engine-cpp, HEAD 37bb518. Read-only survey, 2026-11-04.
Line references are to HEAD 37bb518. `SD` = `R/simulation_driver.R`, `IL` = `R/interim_logic.R`,
`SM` = `R/state_management.R`, `PH` = `R/posterior_helpers.R`, `PD` = `R/posterior_cpp_dispatchers.R`,
`GD` = `R/gate_diagnostics.R`, `CC` = `R/simulate_chunk_cpp.R`, `XC` = `src/simulate_chunk_exact.cpp`.
The SD line numbers in `dev/perf/evidence/chunk_body_inventory.md` (HEAD c2a88d5) are 17 lower than here.

## 0. Short answer for the PI

Nothing intrinsic stops the fast engines from covering almost every configuration. Phase 1 was gated by bit
identity against the ten frozen designs (R01 to R04, R23 to R28), so it only ports and proves the branches those
designs exercise. Every other branch is refused by the predicate and runs the R body. That was a scope choice for
the gate, not a limit of the method.

Of the 32 predicate items, 18 pin a valid configuration that real or default calls can reach. Seventeen of those
need only trial-loop logic, and the phase 2 exact posterior applies to them unchanged. One (the non-PH between-arm
model) needs new posterior math. The rest are two path selectors, three draw-count checks that matter only to
phase 1, and nine error, degenerate, dead or deliberate-R items that should stay on the R body.

A call with every evolveTrial option unset never reaches a fast engine today. On the single-arm path it fails
items 4 (`looks_stop_at_final` defaults to FALSE), 6 (`random_dropout` defaults to TRUE), 7, 8, 10 and 13c, and on
the between-arm path also V2, V6 and V7 (Section 5). The paper's own wrappers fail only those
two items on futility-only designs, because they already set `max_trial_time = 72`, `rebalance_after_events = 50`,
named thresholds and the PH model. The atbo-longleaf certify and search pipelines reach the fast engines only
because their shims set both values (`verify_m4_install("off")`, `p6r_install("p6")`).

## 1. Where the predicate lives

- Common items: `CC:18-56` (`.chunk_cpp_common_checks`). Single-arm items: `CC:59-113` (`.chunk_cpp_hc_checks`).
  Between-arm items: `CC:117-186` (`.chunk_cpp_vsref_checks`). Assembly: `CC:189-196`. Item 11: `CC:213-219`.
  Packing and `cfg$supported`: `CC:236-260`.
- Dispatch: `SD:966-978` (phase 1, `simulate_chunk`) and `SD:990-999` (phase 2 default selection and fallback).
- The exact engine adds no item of its own. It reuses `.chunk_cpp_supported` and `.chunk_cpp_assert_supported`.
  Its C++ re-checks items 4, 6 and 7 as hard stops (`XC:214-216`) and the arm count per mode (`XC:191-198`).
- The two C++ engines carry separate copies of the trial loop (`src/simulate_chunk.cpp`, 837 lines,
  `src/simulate_chunk_exact.cpp`, 1387 lines). Every coverage item below therefore lands twice, unless phase 3
  first factors the trial logic into a shared header.

## 2. Item table

Columns: what the item pins, the R-body branch it excludes with file:line, the callers that reach the excluded
branch, and whether a call with all options unset hits it. "Paper" means `~/Downloads/adaptive-trial-bo-paper`
`{R,scripts}`. "atbo" means `~/Downloads/atbo-longleaf`. BATON's `R/` makes no evolveTrial call (BATON takes the
simulator as a user function), and no `~/research/trials-research/*/scripts` tree calls evolveTrial
(`jasa-baton-paper` is a symlink to the paper repo).

### 2a. Common items (both paths)

| Item | Pins | Excluded branch (R code) | Callers of the excluded branch | Default hits it? |
|---|---|---|---|---|
| 4 | `looks_stop_at_final` TRUE | FALSE branch, the 422eb40 look schedule: a look at the end of accrual, then beats to `max_trial_time` (`SD:698-713`). The TRUE-only caps at `SD:569`, `SD:578` are skipped. | Every paper wrapper (`R/warmstart_wrappers.R`, `R/wrapper_factory.R`, `R/warmstart_core.R`, `R/verification_helpers.R`) leaves the option unset. atbo sets TRUE only through `pilot/p6recert/R/p6r_protocol.R:44` and the shims that source it (e2c, retune, retune2a). | Yes. The default is `getOption("evolveTrial.looks_stop_at_final", FALSE)` (`SD:253`). |
| 5 | `pt_targets_abs` NULL | Person-time milestone looks and the `latest_calendar_look` backstop (`SD:622-645`, `cum_person_time_all_arms`) | `scripts/example_between_arm_comparison.R` (paper, and its atbo copy), one d1-fix test file in atbo docs | No. `person_time_milestones = NULL`. |
| 6 | `random_dropout` FALSE | Keeps the uniform as the stored censoring time (`SD:674-675`) | Every paper wrapper (the formal is not passed). atbo passes FALSE only through `pilot/verify/R/verify_m4_switch.R:42`. | Yes. The formal defaults to TRUE (`SD:288`). |
| 7 | `max_trial_time` finite | Infinite trial time: no post-enrollment looks under item 4 FALSE (`SD:704`), `t_star = last_enroll + F` (`SD:718-722`), final time (`SD:743-747`), person-time gate denominators become Inf, so the fraction is 0 (`IL:202`, `SM:226`) | None found. Paper wrappers set 72 (`row$max_trial_time %||% 72`). | Yes. Default `Inf` (`SD:250`). |
| 8 | `rebalance_after_events` non-NULL | No interval rebalancing: `rebalance_done` TRUE from the start, the block at `SD:589-619` never runs | None found. Paper wrappers set 50. | Yes. Default NULL (`SD:274`). |
| 12 | `show_progress` FALSE | Progress bar ticks per trial (`SD:367-370`, `SD:459`, tick at `SD:899`) | Any interactive call with `parallel_replicates = FALSE` | Yes in an interactive session (`progress = interactive()`, `parallel_replicates = FALSE`). No under Rscript. |
| 13b | `diagnostics` FALSE | Diagnostic `message()` calls (`SD:610`, `IL:101` and others) | None found that reach evolveTrial (the `baton-dosefinding` hits are its own code) | No |
| 13g | well-formed design inputs (finite sorted cutpoints, prior lengths K, positive rates, probs and Weibull parameters) | Malformed or degenerate inputs, where the R body errors or produces NaN (for example `runif(1, 0, Inf)`) | None | No |

### 2b. Single-arm (HC) items

| Item | Pins | Excluded branch (R code) | Callers of the excluded branch | Default hits it? |
|---|---|---|---|---|
| 1 | path selector, `compare_arms_option` not TRUE | none (selects the HC list) | n/a | n/a |
| 2 | exactly one arm | Several arms each against a historical control: arm loop at the final (`SD:768-810`) and interim (`IL:171-260`), arm draw over more than one eligible arm | None found | No (depends on `arm_names`) |
| 3 | `max_total_patients_per_arm` named by `arm_names` | Positional and first-element fallbacks (`IL:203-213`), positional `is_eligible` (`SD:555-557`) | None found in production. Any caller that passes a plain scalar such as 60 | Not a default (the argument is required), but a plain scalar hits it |
| 9 | `efficacy_stopping_rule_hc` not TRUE | Interim efficacy stop (`IL:238-249`) | Paper single-arm wrapper on bidirectional rows (`R/warmstart_wrappers.R:421`, 25 of 39 rows in `scenarios/cohortB_singlearm_production.csv`), `R/wrapper_factory.R:595`, `R/warmstart_core.R:993`, `R/verification_helpers.R:75`, `scripts/dgp_sensitivity_step1_verify.R:119`, `scripts/rule_variant_ablation.R:113` | No. Default FALSE. |
| 10 | success and futility median thresholds named and finite | Final fallbacks Inf and -Inf, giving p = 0 (`SD:787-806`), interim fallback to `null_median_arms` and `futility_median_arms` (`IL:57-61`) | None found. Paper wrappers set both. | Yes. `median_pfs_success_threshold_arms` and `median_pfs_futility_threshold_arms` default NULL (`SD:240-242`). |
| 13a | `min_person_time_frac_per_arm` unnamed 0 | A binding person-time fraction gate, named or positional (`IL:186-195`, `GD:19-77`, `IL:218-223`) | None found. Paper wrappers pass 0. | No. Default 0. |
| 13c | interim `futility_threshold_hc_prob` one finite number | NULL threshold: no interim futility stop (`IL:251-253`) | None found | Yes. Default NULL (`SD:228`). |
| 13d | final probability thresholds finite | NULL makes `if` error (`SD:794`, `SD:805`). Inf is a valid "never" setting. | None | No. Defaults 0.85. |
| 13e | posterior draw counts at least 1 | Zero or non-finite draw counts (`IL:45-46`) | None | No |
| 13f | HC gate values scalars | Vectors longer than one make `coalesce_num` error (`SM:265-266`, `IL:182-184`) | None | No |

### 2c. Between-arm (vs-reference) items

| Item | Pins | Excluded branch (R code) | Callers of the excluded branch | Default hits it? |
|---|---|---|---|---|
| V1 | path selector, `compare_arms_option` TRUE | none (selects the vs-ref list) | n/a | n/a |
| V2 | `use_ph_model_vs_ref` TRUE | Independent-arms model: median-difference probabilities at the interim (`IL:33-39`) and final (`SD:834-840`, `PH:60-65`), sampler `PD:169-196` | None found on the vs-ref path. The paper's `use_ph_model_vs_ref = FALSE` settings are all in single-arm base lists, where the flag is not read. | Yes for any vs-ref call. Default FALSE (`SD:237`). |
| V3 | two arms, one reference, one experimental | More than one experimental arm: per-arm loop (`IL:110-160`) and, without PH, the shared control cache with discarded draws (`IL:116-119`, breakage inventory P6) | None found (all paper `arm_names = c(...)` calls have two arms) | No |
| V4 | `max_total_patients_per_arm` named | Positional `is_eligible` (`SD:555-557`). The vs-ref gate stops on an unnamed vector (`SM:231-232`). | None found | Not a default |
| V5 | `efficacy_stopping_rule_vs_ref` not TRUE | Interim efficacy stop (`IL:146-149`) | Paper multi-arm wrapper on bidirectional rows (`R/warmstart_wrappers.R:789`, 28 of 35 rows in `scenarios/cohortB_multiarm_production.csv`), `R/wrapper_factory.R:492`, `R/warmstart_core.R:1182` | No. Default FALSE. |
| V6 | futility HR margin one finite number, efficacy HR margin NULL | Margin NULL, so `log_margin = 0` (`IL:12-18`, `SD:823-826`). Efficacy margin set (`IL:26-28`, `SD:830-832`). | Efficacy margin: `rule_variant = "symmetric"` (`R/warmstart_wrappers.R:740`, the Tolaney LOI design). Margin 0 is reached by `no_margin`, which is inside the set. | Yes. `compare_arms_hr_margin` defaults to NULL (`SD:235`). |
| V7 | vs-ref probability thresholds finite | NULL thresholds make the final `if` error (`SD:842-844`) | None | Yes, but the R body then errors, so it is not a coverage gap |
| V8 | `compare_arms_target_arm` absent | Target-arm override (`IL:105-109`) | None. `run_simulation_pure` never puts this element in `args` (`SD:393-440`), so the branch is dead from the public entry. | No |
| V9 | vs-ref gates NULL or one unnamed number | Named and positional gate vectors (`GD:27-38`) | None found. Paper wrappers pass scalars. | No |
| V10 | log-HR prior NULL or one finite number | Longer or non-finite values, where `coalesce_num` errors or the prior is degenerate (`PH:113-114`) | None | No. Defaults 0 and 1. |
| 13e (vs-ref) | `num_posterior_draws` at least 1 | Zero or non-finite draws (`SD:802`) | None | No |
| V11a | `num_posterior_draws` integral | Non-integer draw counts (`PH:124`, `PH:127`) | None | No |

### 2d. Dispatch-level items

| Item | Pins | Excluded branch | Callers | Default hits it? |
|---|---|---|---|---|
| 11 | `EVOLVETRIAL_USE_CPP` in TRUE, 1, YES | R posterior fallbacks (`PD:11-16`), bit-different from C++ in `ph_beta_mode_var` | None. No caller sets the variable. | No. Unset means TRUE. |
| cfg | `cfg$supported` (packing succeeded) | Catch-all for a packing error (`CC:242-245`) | follows from the other items | No |

Outside the predicate entirely, and not counted: the hybrid and seamless simulators (`compute_hybrid_oc_rcpp`,
`run_hybrid_simulations_rcpp`, used by `R/warmstart_wrappers.R:1393-1430` and `R/hybrid_dual_wrappers.R`), and the
binary simulator (`run_simulation_binary`, six paper scripts). See `2026-11-04-hybrid-seamless-speedup-note.md`.

## 3. Classification, effort and gate

The gate for every coverage item has two layers, as in phases 1 and 2. First, phase 1 bit identity: a new
acceptance design that exercises the branch, `identical()` against the R body on the same R stream (layers 0 and 1,
then layer 2 at production R). Second, for the exact engine, a fresh-run statistical equivalence check against the
R body (operating characteristics within Monte Carlo tolerance, plus the paired decision-disagreement diagnostic).
Items with new posterior math also need the P2-1 accuracy gate (absolute 1e-7 against a deterministic reference,
and within 4 Monte Carlo standard errors against 10^7 draws).

| Item | Exact posterior | Effort | Reason |
|---|---|---|---|
| 4 | applies directly | small | A second look schedule in the trial loop. No new posterior. The FALSE branch is the simpler of the two. |
| 6 | applies directly | small | The uniform is already drawn for P3. Store it instead of Inf. |
| 7 | applies directly | small | Three `is.finite` branches and an Inf gate denominator (the fraction becomes 0). |
| 8 | applies directly | small | Removes a branch. The trial runs with fixed cutpoints. |
| 9 | applies directly | small | The interim efficacy tail is the same gamma-sum tail the final already uses, with the success threshold. |
| V5 | applies directly | small | The interim efficacy tail is the normal tail already used at the final. |
| V6 | applies directly | small | Pack `log_margin = 0` when NULL and a finite `log_eff` when set. Both are thresholds on the same normal tail. |
| 10 | applies directly | small | Pack the fallback thresholds. A non-finite threshold gives p = 0 without a posterior call. The phase 1 port must still draw. |
| 13c | applies directly | small | A NULL threshold means no interim futility stop. |
| 12 | no math | small | Tick per block, or once at the end. Output unchanged. A UI decision for the PI. |
| 3, V4 | applies directly | small | Normalize the vector at pack time, mirroring `IL:203-213` and the positional `is_eligible`. |
| 13a | applies directly | small | The gate logic of `GD:19-77` and `IL:218-223`. No posterior. |
| V9 | applies directly | small | The named and positional branches of `GD:27-38` at pack time. |
| 2 | applies directly, per arm | medium | The trial loop over n arms, the arm draw over several eligible arms, per-arm final tallies. The return shape is already per arm. |
| V3 | applies directly, per pair (PH) | medium | Loop over experimental arms. Without PH, phase 1 must also reproduce the discarded control-cache draws. |
| 5 | applies directly | medium | Person-time milestone schedule, `cum_person_time_all_arms` and the backstop. No posterior. |
| V2 | needs new math | medium to large | P(med_T > med_C + m) and P(med_T <= med_C - m) under independent gamma posteriors on each arm's hazards. med_T > x exactly when the cumulative hazard H_T(x) < log 2. The probability is therefore a one-dimensional Stieltjes integral of S_T(x + m) against the distribution of med_C. Both functions are tails the phase 2 code already computes, but each integrand evaluation is a characteristic-function inversion, so it needs adaptive quadrature with error control. Note that the interim uses `diff_med > 0` (unbuffered) and `diff_med <= -margin`, while the final uses `med_T > med_C + margin`. Both forms must be reproduced. |
| 13e, 13e (vs-ref), V11a | irrelevant to exact (no draws) | small | Phase 1 needs them. The exact engine could drop them, except as input validation. |
| 13b | keep R body | n/a | Message side effects. No callers. |
| 13d, 13f, V7, V10, 13g | keep R body | n/a | Error or degenerate paths. The R body produces the reference error. 13d with Inf could be admitted, small. |
| V8 | keep R body | n/a | Dead from the public entry. |
| 11 | policy ruling | small | Phase 1 must keep the R fallback. The exact engine draws nothing, but it calls `ph_beta_mode_var_cpp`, so admitting it changes the last bits of the Laplace mode. A ruling decides whether the variable binds the exact engine. |
| cfg | follows from the others | n/a | Catch-all. |
| 1, V1 | path selectors | n/a | Not exclusions. |

## 4. Ordered phase 3 plan (caller impact per effort)

1. **Items 4 and 6** (small, small). These make every paper wrapper run fast on futility-only designs (14 of 39 single-arm and 7 of 35 multi-arm production rows), without the
   atbo shims. The phase 1 layer reproduces the 422eb40 look schedule and the 7565a17 dropout model bit for bit.
   This step does not decide which model is scientifically right. It only makes both fast.
2. **Items 7 and 8** (small). Infinite trial time and no rebalancing. With steps 1, 2 and items 10 and 13c from
   step 4, a bare single-arm `run_simulation_pure()` call is covered. A bare between-arm call also needs V2 and V6,
   and V7 errors in the R body.
3. **Items 9, V5 and V6** (small each). These cover the bidirectional designs (25 of 39 single-arm and 28 of 35
   multi-arm production rows) and the symmetric rule variant. The exact tails already exist.
4. **Items 10, 13c, 3, V4, 13a, V9** (small each). Defaults and input-shape cleanup. No known production caller,
   but each widens what a new user's call reaches.
5. **Item 12** (small, PI ruling on ticks). This covers interactive sequential calls.
6. **Shared trial-loop refactor** (medium). Worth doing before step 7, because items 2, V3 and 5 each add loop
   structure that otherwise lands twice.
7. **Items 2 and V3** (medium each). Multi-arm single-arm-vs-historical and more than one experimental arm. No
   caller found today. These matter for future platform designs.
8. **Item 5** (medium). Person-time milestones. One example script uses it.
9. **Item V2** (medium to large, new math). A default vs-ref call hits it, but no production caller uses it.
   It needs the P2-1 accuracy gate in addition to the two layers.
10. **Item 11** (policy ruling, then small). Optional.

Not scheduled: 13b, 13d, 13f, V7, V10, 13g and V8 stay on the R body. The hybrid, seamless and binary simulators
are separate engines and are outside this inventory.

## 5. Counts

- Predicate items: 32 (8 common, 10 single-arm, 12 between-arm, item 11, `cfg$supported`). The exact engine adds 0.
- Coverage targets: 18 (4, 5, 6, 7, 8, 12, 2, 3, 9, 10, 13a, 13c, V2, V3, V4, V5, V6, V9). 17 apply the exact
  posterior directly or need no posterior. 1 (V2) needs new math. 12 are small, 5 medium, 1 medium to large.
- Draw-count items relevant only to phase 1: 3 (13e, 13e vs-ref, V11a).
- Kept on the R body (error, degenerate, dead, deliberate or catch-all): 9 (13b, 13d, 13f, 13g, V7, V8, V10,
  11, cfg). Item 11 has an optional policy ruling.
- Path selectors: 2 (1, V1).
- Unclassified: 0.
- Default call (all options unset) fails: 4, 6, 7, 8, 10, 13c on the single-arm path, plus V2, V6 and V7 on the
  between-arm path (V7 then errors in the R body). In an interactive session it also fails 12. The paper wrappers
  fail only 4 and 6 on futility-only designs.

## 6. Where the brief conflicts with the code

1. "P6" is overloaded. In the phase 1 breakage inventory, P6 is the discarded control-cache draws
   (`IL:116-119`). The predicate labels `looks_stop_at_final` as item 4. The "(P6)" comments in `SD:563-578`
   and `SD:715` refer to the atbo-longleaf `p6r_install("p6")` protocol.
2. `~/research/atbo-longleaf` does not exist. I used `~/Downloads/atbo-longleaf` (HEAD 1416e3d6).
3. The predicate spans more than `CC:150-240`. The common and single-arm lists start at `CC:18`.
4. The exact engine has no extra supported-set conditions. `XC:191-216` only re-checks items 2, V3, 4, 6 and 7.
5. BATON's `R/` makes no evolveTrial call, and no trials-research script calls evolveTrial. The callers that
   matter are the paper wrappers and the atbo pipeline shims.

## 7. Merging main (4dc2f1f, amd8 chronology 56f0328) before phase 3

Sources: `git show 56f0328 -- R/simulation_driver.R`, main's `R/simulation_driver.R` at 4dc2f1f (cited as
`MAIN:`), and this branch at 37bb518. The merge base is 422eb40. Main does not contain c124dca to 58a4850, so main
has neither `looks_stop_at_final` nor `random_dropout`. Both branches call themselves 0.0.0.9002, which is a version
collision the merge must resolve.

### 7a. Unconditional or behind an option, and which frames change

56f0328 is unconditional. It edits the one R body with no option and no version switch. It adds four rules and
one assertion block.

1. Accrual closes at `max_trial_time`. An arrival at or after a finite cutoff enrolls nobody, fires no look, and
   sets `current_time <- max_trial_time` (`MAIN:523-529`).
2. No accrual-phase look at an arrival known to reach T_final. That is (a) `min_follow_up_at_final <= 0` and the
   arrival fills the last open slot, or (b) vs-ref mode with no eligible experimental arm and
   `current_time >= T_final` (`MAIN:531-551`, the `interim_ok_here` flag).
3. The look at the moment enrollment ends runs only if `current_time < T_final` (`MAIN:665-677`). Post-enrollment
   looks start at `current_time` and stop strictly before T_final, not before `max_trial_time` (`MAIN:681-690`).
4. Final time is T_final = min(last enrollment + F, `max_trial_time`) (`MAIN:694-702`, same formula as before).
5. Each trial errors if an enrollment, stop or final time exceeds a finite cutoff, or an enrollment follows the
   final time (`MAIN:704-718`).

Relative to 58a4850 with P6 off (the 422eb40 schedule, which is the default): amd8 changes results on a frame when
accrual reaches a finite cutoff, when any post-enrollment look would land in [T_final, `max_trial_time`), when
F <= 0 (the end-of-enrollment look and the last-slot look are dropped, even with `max_trial_time = Inf`), or in
vs-ref case (b). With the default F = 0 and a finite cutoff, that is nearly every trial that finishes enrollment
with an arm still recruiting. Unchanged frames are those with `max_trial_time = Inf` and F > 0, and those where the
cutoff never binds and no look reaches T_final (main's tests at `tests/testthat/test-chronology-amd8.R:66` and `:85`).

Relative to 58a4850 with P6 on: amd8 has no option, so the comparison depends on how the merge combines the two
branches. Taken as written, amd8 differs from the P6 chronology in four places, listed in 7b. On the frames the
predicate accepts (P6 on, finite cutoff, F = 6 in the certify shims) at least difference 2 below fires in almost
every trial, and difference 3 fires in almost every vs-ref trial.

### 7b. Subsumes or conflicts with P6

Neither. The two agree on the core: accrual closes at the cutoff, T_final has the same formula, and interims run
strictly before T_final. They conflict in four details.

1. **Clock after a cutoff.** P6 resets `current_time <- last_enroll_time` (`SD:569-571`), so post-enrollment looks
   at last enrollment + j * beat < T* still run. amd8 sets `current_time <- max_trial_time`, so no post-enrollment
   look runs. They differ whenever accrual is truncated by the cutoff.
2. **Look at the moment enrollment ends (j = 0).** P6 never runs it ("No look at j = 0", `SD:715-717`, following
   ADAPT Amendment 2). amd8 runs it whenever `current_time < T_final`, which holds for every F > 0. That is a
   second look at the last arrival time, with one more patient counted (zero follow-up), and in phase 1 it
   consumes posterior draws.
3. **Phantom vs-ref arrival.** P6 lets the arrival that would leave the loop through the vs-ref break enroll
   nobody and fire neither rebalancing nor a look (`SD:578-589`). amd8 still rebalances and looks at that arrival
   unless `current_time >= T_final` (`MAIN:546-551`), then the post-enrollment grid starts from the phantom time.
   In vs-ref mode the reference arm never stops (P19), so essentially every vs-ref trial ends this way.
4. **Last-slot rule for F <= 0.** amd8 suppresses the look at the arrival that fills the last slot. P6 has no such
   rule.

A PI ruling is needed on which chronology is the target: P6 (ADAPT Amendment 2, what the engines and the
atbo certify runs implement), amd8 (BATON-C freeze Amendment 8), or both as named modes.

### 7c. Silent divergence after a merge (the key hazard)

Yes, unless the resolution leaves the P6 branch of the R body textually as it is at 37bb518. The predicate admits
a frame on `looks_stop_at_final` TRUE and a finite cutoff (`CC:26`, `CC:30`). It says nothing about chronology
mode. Dispatch picks the C++ or exact engine with no runtime comparison to the R body (`SD:973`, `SD:990-999`).
The P6 chronology is hard-coded in the C++ at these lines:

- phase 1: cutoff with clock reset `src/simulate_chunk.cpp:616-619`, phantom vs-ref arrival `:621-626`,
  post-enrollment grid from last enrollment with no j = 0 look `:698-705`, final time `:723`;
- exact: cutoff `src/simulate_chunk_exact.cpp:881`, post-enrollment grid `:954-960`, final time `:979`.

Any resolution that puts amd8's rules 1 to 3 on the P6-on path (for example the cutoff block of `MAIN:523-529`
placed ahead of the P6 block, or `MAIN:675` replacing the absent j = 0 look) makes the R body diverge from both
engines on accepted frames, with no error at run time. The C++ also lacks amd8's assertions (`MAIN:704-718`).
`tests/testthat/test-engine-equivalence.R`, `test-simulate-chunk-cpp.R` and `test-simulate-chunk-cpp-vsref.R`
compare phase 1 against the R body and would fail, so the divergence shows up at test time if the suite runs
after the merge. The golden-output tests would still pass, because they pin the engines and not the R body.

A safe merge does one of two things. It applies amd8 only on the P6-off path, so the accepted frames keep the P6
R body. Or it adds a chronology-mode item to the predicate that refuses amd8 frames until the engines port them.

### 7d. Effort to port amd8 chronology

- **Phase 1 engine:** small to medium. Four local edits to the trial loop: clock to the cutoff, the conditional
  j = 0 look, look and rebalancing at the phantom arrival, and the last-slot test. The last-slot test uses
  `rep_len(unname(max_total_patients_per_arm), ...)`, which is positional, but the predicate already pins named
  vectors (items 3 and V4). The assertions are cheap. The bit-identical gate must reproduce the RNG order of the
  extra looks, so it needs new acceptance frames that bind the cutoff, use F = 0 and F > 0, and cover vs-ref
  phantom exits, then layers 0, 1 and 2.
- **Exact engine:** small to medium. The same loop edits in `src/simulate_chunk_exact.cpp`, with no new posterior
  math, then the fresh-run statistical-equivalence gate.
- Together this is about the size of item 4. If amd8 replaces the 422eb40 schedule as the P6-off default after the
  merge, then item 4 in Section 4, step 1, becomes "port amd8" rather than "port 422eb40". Doing the
  shared-trial-loop refactor (step 6) first would land each edit once instead of twice.
