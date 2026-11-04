# P3-2: bidirectional designs and the symmetric rule variant (items 9, V5, V6)

Date 2026-11-04. Worktree `~/research/evolveTrial-merge`, branch `perf/merge-main`, base 8e1430f (0.0.0.9004). Plan
`dev/plans/2026-11-04-p3-2-plan.md`. Nothing committed, pushed or installed into `~/R/library`. Builds (`p32_install.sh`):
`dev/perf/lib/p32_dev` (working tree before the test edits and the version bump, used by gate B, layer 1 and the gate S
pilot), `dev/perf/lib/p32_final` (final working tree, 0.0.0.9005, used by gate T, goldens, P6-on recheck, G3),
`dev/perf/lib/p32_base` (git archive of 8e1430f), `dev/perf/lib/p32_mut_{m,v}` (mutation builds). The C++ and R sources
are the same in p32_dev and p32_final; only tests, DESCRIPTION and NEWS changed between them. Scripts and logs
`dev/perf/phase3/p32_*`.

## Counts

| Item | Result |
|---|---|
| Gate B (`p32_gateB.log`, `p32_gateB.rds`, p32_dev) | 24 frames, 3,408 chunks, 10,712 trials, 0 mismatches, 0 chunks ending in an R-body error |
| Gate B, rule TRUE with NULL / Inf threshold equals rule FALSE (C05, C06) | 284 of 284 chunks identical (result and stream) |
| Gate B mutation check (`p32_mutation.sh`) | interim success value swapped for the final one in phase 1 (C07): 78 of 142 chunks mismatch; vs-ref interim efficacy check disabled (D10): 124 of 142 mismatch |
| Layer 1, old designs (`p32_layer1.log`, p32_dev) | 2,120 chunks, 0 mismatches |
| Gate T, full suite (`p32_tests_final.log`, p32_final) | 1,792 passed, **1 failed**, 37 skipped, 0 errors, 0 warnings (failure: Section 4, item 1) |
| Long tier, `test-engine-equivalence.R` (`p32_tests_long.log`) | 812 passed, 0 failed |
| Goldens (`p32_golden.log`, `p32_golden_identical.log`) | 18 of 18 reproduced bit for bit by p32_final; the 12 existing entries stored byte-identical to the pre-P3-2 file (`p32_engine_golden_before.rds`); 6 new entries |
| P6-on recheck vs 37bb518 (`p32_onpath.log`) | R body 15/15, cpp 15/15, default 12/15 (the same 3 dropout frames as P3-1); all three 15/15 identical to the P3-1 build |
| G3 worker count (`p32_g3.log`, `p32_g3_results.csv`) | C02 (R01 bidirectional, P6 off) and D03 (R23 symmetric): 20 runs, 0 failed |

Suite reconciliation against P3-1 (1,516 passed): `test-engine-equivalence.R` 340 to 499 (+159), `-cpp-vsref` 234 to
246 (+12), `-cpp` 158 to 183 (+25), `-exact` 282 to 363 (+81), `test-return-variance.R` 47 to 46 (the failure), total
1,516 + 277 = 1,793 = 1,792 passed + 1 failed. The 4 skips in `test-engine-equivalence.R` are the long-tier tests (P3-2
adds one). The P3-1 runner attached the package namespace; my first `p32_tests.R` run did not, so two
`test-hybrid-simulator.R` tests that call an internal function errored (also on the 8e1430f base build). The runner now
evaluates tests in a child of the namespace and that file passes 60 of 60.

## 1. What changed (diff vs 8e1430f)

| File | Lines (+/-) | Change |
|---|---|---|
| `R/simulate_chunk_cpp.R` | +57 / -13 | Item 9: rule TRUE admitted when `efficacy_threshold_hc_prob` is NULL or one number (a longer vector errors in the R body's `&&`, a non-numeric one has R's `is.finite` semantics); the rule must be the same in args and the frame. New `.chunk_cpp_eff_interim_hc()`: rule TRUE and threshold one finite number. Item 10: with that flag TRUE the interim success value (IL:57-60) must be numeric, name the arm and not be NA (+-Inf allowed). V5: rule TRUE admitted (V7 already requires a finite threshold). V6: each margin NULL or one finite number. Packer: `eff_interim`, `eff_prob_interim`, `success_thr_interim` (HC), `eff_interim`, `vr_eff_prob_interim` (vs-ref); values NA and unread when the flag is FALSE |
| `src/simulate_chunk.cpp` | +49 / -5 | `interim_hc`: `pr_eff = mean(med > success_thr_interim)` on the same draws, checked before futility, sets STOPPED_EFFICACY, `stop_eff`, `final_n_cur`, `stop_time`, `continue`; runs also when interim futility is NULL. `interim_vsref`: `pr_eff = mean(logHR < log_eff_interim)` before `ph_ge`. `unpack` refuses `eff_interim` TRUE with a non-finite probability or an NA success value |
| `src/simulate_chunk_exact.cpp` | +89 / -1 | Policy calls `hc_interim_begin / hc_interim_above / hc_interim_below` (McPost: one draw, two tallies; ExactPost: the gamma-sum tail `Pr(median > s)` the final uses, `0` for `s = Inf` and `1` for `s = -Inf` without a posterior call; TracePost: kind 6 for the efficacy tail). vs-ref interim: `post.ph_lt(log_eff_interim)` first. McPost's `ph_lt` and `ph_ge` tally the stored draws and ExactPost's evaluate the stored normal, so neither draws again. Same `unpack` refusals |
| `tests/testthat/test-engine-equivalence.R` | +153 / -10 | P3-2 frames, routine (R01, R23) and long (R01, R03, R23) tiers: cpp identical to the R body, default = exact, no NA, power = PET_Efficacy + Pr_Final_Efficacy; efficacy stops taken and `Var_PET_Efficacy = p(1-p)/n > 0` in both engines; rule TRUE with NULL / Inf / NA threshold identical to rule FALSE (R body, cpp, exact) and packed `eff_interim` FALSE; NULL or NA interim success value with the rule on stays on the R body and errors with its message; 6 golden entries; still-excluded list now 2, 5, 9 (non-numeric threshold), 13a, V2, V6 (two-element margin), V9; the fallback test uses 13a and a two-element margin instead of items 9 and V6 |
| `tests/testthat/test-simulate-chunk-cpp.R` | +32 / -7 | Item 9 refusals (length 2, non-numeric, args-only rule), 10e (NULL success value with the rule on), 10f (NA value); the "no HR margin" frame now runs on C++ and equals the R body, and a non-PH frame (V2) checks the fallback |
| `tests/testthat/test-simulate-chunk-cpp-vsref.R` | +8 / -4 | Item 9 label; V5 refusal (args-only rule), V6 refusals (NA, length 2, non-numeric) |
| `tests/testthat/test-simulate-chunk-exact.R` | +75 / 0 | Replay instantiation = R body on 8 efficacy frames x 2 scenarios (efficacy stops taken); M1 trace: kind-6 records use `success_thr_interim`, equal `.exact_pr_median_above`, precede every interim futility record on the same inputs, and a crossing is never followed by a futility record; PH interim `ph_lt` before `ph_ge` on the same posterior |
| `tests/testthat/fixtures/engine_golden.rds` | binary | 12 old entries copied verbatim, 6 added |
| DESCRIPTION, NEWS.md | +1/-1, +41 | 0.0.0.9005, NEWS section |

`src/symbols.rds` shows as deleted in `git status`; it was already deleted before P3-2 started (initial status).

## 2. Gate B detail

Frames from R01, R03, R23 (`p32_frames.R`, same spec format as `p31_frames.R`), masters 3210001 and 3210002 x workers
16 and 19 plus one long chunk per scenario, every frame with `return_variance = TRUE`. Per chunk: cpp = R body (result,
`.Random.seed`, error), the instrumented R body = R body, the exact file's replay instantiation = phase 1, no NA; for
C05 and C06 also R body (frame) = R body (rule FALSE). Branch counts (`p32_instrument.R`): stops are trials summed from
the R body's chunk results; crossings are looks (exit traces on `calculate_current_probs_hc` / `_vs_ref` with the R body's
guard); both-held is trials / looks where both tails crossed at one look; inert is looks with the rule TRUE and the
threshold NULL or Inf.

| Frame | Chunks | Trials | Mismatch | same_as | Stops eff / fut | Looks crossed eff / fut | Both held (trials / looks) | Inert looks |
|---|---|---|---|---|---|---|---|---|
| C01 R01 P6 on, eff 0.894 | 142 | 560 | 0 | - | 273 / 287 | 273 / 289 | 2 / 2 | 0 |
| C02 R01 P6 off | 142 | 560 | 0 | - | 273 / 287 | 273 / 289 | 2 / 2 | 0 |
| C03 R01 on, eff 0.3 fut 0.3 | 142 | 560 | 0 | - | 454 / 106 | 454 / 282 | 176 / 176 | 0 |
| C04 R01 off, eff 0.3 fut 0.3 | 142 | 560 | 0 | - | 454 / 106 | 454 / 282 | 176 / 176 | 0 |
| C05 R01 off, rule TRUE thr NULL | 142 | 560 | 0 | 142/142 | 0 / 326 | 0 / 326 | 0 / 0 | 3,262 |
| C06 R01 on, rule TRUE thr Inf | 142 | 560 | 0 | 142/142 | 0 / 329 | 0 / 329 | 0 / 0 | 2,967 |
| C07 R01 off, interim 5.0 vs final Inf | 142 | 560 | 0 | - | 237 / 318 | 237 / 318 | 0 / 0 | 0 |
| C08 R01 on, interim 5.0 vs final Inf | 142 | 560 | 0 | - | 252 / 306 | 252 / 306 | 0 / 0 | 0 |
| C09 R01 off, eff on, futility prob NULL | 142 | 560 | 0 | - | 344 / 0 | 344 / 0 | 0 / 0 | 0 |
| C10 R01 on, eff on, futility prob NULL | 142 | 560 | 0 | - | 342 / 0 | 342 / 0 | 0 / 0 | 0 |
| C11 R01 off, interim success Inf | 142 | 560 | 0 | - | 0 / 326 | 0 / 326 | 0 / 0 | 0 |
| C12 R03 on, eff 0.974 | 142 | 560 | 0 | - | 268 / 278 | 268 / 278 | 0 / 0 | 0 |
| C13 R03 off, eff 0.8 | 142 | 560 | 0 | - | 313 / 238 | 313 / 256 | 18 / 18 | 0 |
| D01 R23 off, current | 142 | 312 | 0 | - | 177 / 134 | 177 / 134 | 0 / 0 | 0 |
| D02 R23 on, current | 142 | 312 | 0 | - | 177 / 133 | 177 / 133 | 0 / 0 | 0 |
| D03 R23 off, symmetric 0.1 | 142 | 312 | 0 | - | 166 / 131 | 166 / 131 | 0 / 0 | 0 |
| D04 R23 on, symmetric 0.1 | 142 | 312 | 0 | - | 165 / 137 | 165 / 137 | 0 / 0 | 0 |
| D05 R23 off, no_margin (both NULL) | 142 | 312 | 0 | - | 181 / 126 | 181 / 126 | 0 / 0 | 0 |
| D06 R23 on, both margins 0 | 142 | 312 | 0 | - | 177 / 126 | 177 / 126 | 0 / 0 | 0 |
| D07 R23 off, eff margin only | 142 | 312 | 0 | - | 157 / 129 | 157 / 129 | 0 / 0 | 0 |
| D08 R23 on, margins 0.97 (capped) | 142 | 312 | 0 | - | 0 / 312 | 0 / 312 | 0 / 0 | 0 |
| D09 R23 off, margins -0.2 (clamped) | 142 | 312 | 0 | - | 181 / 126 | 181 / 126 | 0 / 0 | 0 |
| D10 R23 off, symmetric, thresholds 0.45 | 142 | 312 | 0 | - | 255 / 57 | 255 / 76 | 19 / 19 | 0 |
| D11 R23 off, symmetric, efficacy off (V6 alone) | 142 | 312 | 0 | - | 0 / 160 | 0 / 160 | 0 / 0 | 0 |

Consistency of the counts: in every frame futility stops = futility crossings minus both-held looks (289 - 2 = 287,
282 - 176 = 106, 256 - 18 = 238, 76 - 19 = 57), so whenever both tails crossed the arm stopped for efficacy. D08 (eff
margin capped at 0.95, log 0.05) never reaches the efficacy threshold, as expected; D09 gives the same counts as D05
(both clamp to margin 0). All efficacy-stop frames include null and alternative scenarios. Every branch the plan lists
was hit at least once. No branch had to stay excluded for RNG reasons.

## 3. Gate S preparation (not run)

Frames (`p32_frames.R`): C01 HC bidirectional P6 on, C02 HC bidirectional P6 off, D01 vs-ref bidirectional `current`
(P6 off, the default chronology), D03 vs-ref bidirectional `symmetric` (margins 0.1, P6 off), each under scenario 1
(null) and 2 (alternative): eight frame x scenario cells. Scripts: `p32_gateS.sh` (runs), `p32_gateS_run.R` (copy of
p31 run script; also checks `cfg$eff_interim` is TRUE), `p32_gateS_compare.R` (p31 test unchanged), `p32_gateS_pilot.sh`,
`p32_gateS_size.R`. `p32_gateS.sh` uses `LIB=${LIB:-dev/perf/lib/p32_final}`; point LIB at the build of the commit.

**Pilot** (`p32_gateS_pilot.log`, seeds 3220001-4 R body, 3221001-4 exact; 2,000 R-body and 20,000 exact trials per
frame, p32_dev, machine loaded by gate B and tests): R body / exact type I error C01 0.127 / 0.114, C02 0.117 / 0.113,
D01 0.257 / 0.244, D03 0.144 / 0.143; power 0.83 / 0.84, 0.84 / 0.83, 0.85 / 0.86, 0.78 / 0.78; sd of N under null /
alternative 3.7 / 4.3, 3.5 / 4.3, 8.6 / 6.5, 8.9 / 8.1. The pilot is not a readout.

**Sizing** (`p32_gateS_size.R`, `p32_gateS_sizing.tsv`). Test power at a true difference of 0 is
`1 - 2 * pnorm(1.645 - m / se)`; at least 0.95 needs `m / se >= 1.645 + 1.960 = 3.605`. With n1 R-body and n2 = k n1 exact
trials, `se^2 = s2 (1 + 1/k) / n1`, so `n1 >= s2 (1 + 1/k) (3.605 / m)^2`, where s2 = p(1-p) (p the pilot value closest to
0.5 of the two engines) or the larger pilot per-trial variance of N, and m = 0.01 (power), 0.005 (type I), 1% of the
R-body EN. k = sqrt(cost ratio) from the pilot walls, bounded to [4, 16]. Then x1.2 for pilot uncertainty, rounded up to
R-body runs of 5,000 and exact runs of 50,000.

| Frame | k | n1 needed: power / type I / EN0 / EN1 | Binding | R body | Exact | Predicted test power: power / type I / EN0 / EN1 |
|---|---|---|---|---|---|---|
| C01 | 6.1 | 21,285 / 66,831 / 18,649 / 11,932 | type I | 85,000 (17 runs) | 550,000 (11) | 1.00 / 0.985 / 1.00 / 1.00 |
| C02 | 6.5 | 20,749 / 61,994 / 17,025 / 12,048 | type I | 75,000 (15) | 500,000 (10) | 1.00 / 0.980 / 1.00 / 1.00 |
| D01 | 12.8 | 17,417 / 106,999 / 17,391 / 9,519 | type I | 130,000 (26) | 1,700,000 (34) | 1.00 / 0.980 / 1.00 / 1.00 |
| D03 | 13.9 | 23,973 / 68,480 / 18,273 / 13,323 | type I | 85,000 (17) | 1,200,000 (24) | 1.00 / 0.982 / 1.00 / 1.00 |

Type I error binds everywhere: these bidirectional frames have type I error 0.11 to 0.26 at the captured efficacy
thresholds, and the margin is 0.005. 154 runs, seeds 3230101-3230417 (R body) and 3231101-3231424 (exact), no duplicates
(checked). **Estimated wall time**: 1.88 run-hours at 4 workers per run (R body 1.73, exact 0.15; D01 and D03 R body
dominate at about 25-28 ms per trial), so about 0.5 h with `xargs -P 4` (16 cores). The pilot ran on a loaded machine;
allow up to 1 h.

Classification check (0.80 / 0.10) is part of the pass rule: power in D03 is about 0.78 and in C01, C02, D01 0.83 to 0.85,
and every type I error is above 0.10, all well away from the cut points by many standard errors at these sizes.

## 4. Mismatches and open items

1. **Gate T is not 0 failures.** `test-return-variance.R:242` ("run_scenarios with return_variance=TRUE adds Var_
   columns", base design with `efficacy_stopping_rule_hc = TRUE`, now admitted and run by the exact engine) asserts every
   `Var_` column is >= 0. The exact engine's outcome at seed 71 is 7 efficacy and 93 futility stops in 100 trials, and
   the R body's summary (`R/simulation_driver.R:1221`, not editable) computes `Pr_Reach_Max_N = 1 - 0.07 - 0.93 =
   -1.1e-16`, hence `Var_Pr_Reach_Max_N = -1.1e-18`. The R body and cpp draw 8 / 92 at this seed and give exactly 0. It is
   a floating-point artifact of the R summary exposed by different counts, not an engine error (the same counts through
   the R body would give the same negative). I did not edit this test or the R body: the choice is the coordinator's
   (for example a tolerance such as `>= -1e-12` in that assertion, or a clamp in the summary, which is R-body code).
2. Plan Section 2 matches the R body as read (IL:43-66, 146-166, 238-260; SD:851-937; summary SD:1211-1254). Two points
   the plan leaves open, resolved by the R body: the interim and final success values can differ only when
   `median_pfs_success_threshold_arms` is NULL (the interim reads `null_median_arms`, the final falls back to Inf), so
   C07/C08 use that; and on the vs-ref path V7 already requires a finite efficacy probability threshold, so the "rule
   TRUE with NULL / Inf threshold" frames exist only on the HC path.
3. Item 9 admits the rule FALSE with any threshold value (the R body never reads it), and refuses the rule TRUE with a
   non-numeric threshold even though the R body runs it (never stops). A `TRUE` logical threshold would be "finite" to
   R and compare as 1; it is refused as non-numeric.
4. Gate R (Codex review) and the commit are the coordinator's. Scratch libraries `dev/perf/lib/p32_*` and builds
   `dev/perf/build/p32_*` can be deleted after gate S.

## 5. Summary fix (PI ruling 2026-11-04, commit after a46d34f)

Gate T's one failure (`test-return-variance.R:242`) was a defect in the shared R summary, not the port.
`Pr_Reach_Max_N = 1 - early_eff - early_fut` gives -1.1e-16 at 7 + 93 early stops in 100 trials, so
`Var_Pr_Reach_Max_N` was -1.1e-18. The R body gives the same value at those counts. The phase 1 engine landed on
8 + 92, which cancels exactly. Codex flagged it independently on a46d34f (P2). PI chose the summary fix.

- `R/simulation_driver.R`: `Pr_Reach_Max_N` from the counts, `(n - stop_eff - stop_fut) / n`, and `prop_var` clamps
  p to [0, 1]. The same summary serves every engine, so gate B (engine vs R body) is unaffected.
- Goldens regenerated. Per entry only `Pr_Reach_Max_N` and `Var_Pr_Reach_Max_N` moved, max |diff| 1.11e-16
  (16 of 18 entries; `exact_R03_p31`, `cpp_R03_p31` unchanged). `.Random.seed` identical 18/18. After
  regeneration `p32_golden_identical.R` gives 18/18 table and RNG identical.
- Gate T on build `p32_fix`: 1,793 passed, 0 failed, 37 skipped.
