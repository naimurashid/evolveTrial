# P3-1: general loop, first slice (items 4, 6, 7, 8, 10, 13c)

Date 2026-11-04. Worktree `~/research/evolveTrial-merge`, branch `perf/merge-main`, base c936bde (0.0.0.9003). Plan
`dev/plans/2026-11-04-p3-1-plan.md`. Nothing committed, pushed or installed into `~/R/library`. Final build
`dev/perf/lib/p31_final` (working tree after the Codex fixes), gate build `dev/perf/lib/p31_dev` (before the Codex fixes
and before a comment-only pass over the C++ line references), baseline `dev/perf/lib/p31_base` (git archive of HEAD).
Version 0.0.0.9004. Scripts and logs `dev/perf/phase3/p31_*`.

## Counts

| Item | Result |
|---|---|
| Gate B, all frames (`p31_gateB.log`, build p31_dev) | 29 frames, 4,118 chunks, 13,760 trials, 0 mismatches. 216 chunks end in the R body's chronology error (F < 0 frames); every engine stops with the same message and the same stream |
| Gate B, frames changed after the first run (`p31_gateB_rerun.log`, p31_final) | B09, B09c (new), B11 (changed): 426 chunks, 1,680 trials, 0 mismatches |
| Gate B subset after the Codex P2 fix (`p31_gateB_final.log`, p31_final) | B01, B02, B05, B20, B23: 710 chunks, 2,304 trials, 0 mismatches |
| Per chunk, gate B also checks | instrumented R body = R body; the exact file's replay instantiation = phase 1 (result and stream); no NA in the phase 1 output |
| Layer 1, old designs (`p31_layer1.log`, `p31_layer1_final.log`) | 2,120 chunks, 0 mismatches on p31_final, both before and after the Codex fix |
| Layer 1, new frames | gate B is the layer 1 procedure (masters 3110001 and 3110002 x workers 16, 19, plus one long chunk per scenario) on the P3-1 frames |
| Gate S (`p31_gateS_report.tsv`) | 5 of 5 frames pass after one extension of the sample; first readout 2 of 5 (Section 2) |
| Gate T, full suite (`p31_tests_final.log`, p31_final) | 1,516 passed, 0 failed, 36 skipped, 0 errors, 0 warnings |
| Long tier, `test-engine-equivalence.R` (`p31_tests_long.log`, before the Codex fix) | 515 passed, 0 failed |
| New tests on the old build (`p31_tests_on_base.log`) | 12 failures and 5 errors in `test-engine-equivalence.R` (the new tests detect the change) |
| Codex P3 mutation (exact engine without the j = 0 look, scratch build) | 12 failures, all in the new replay test (`p31_mut_j0_tests.log`) |
| Goldens (`p31_golden_identical.log`) | 12 of 12 reproduced bit for bit by p31_final; the 6 existing entries are stored byte-identical to the pre-P3-1 file |
| G3 worker count, 2 new frames (`p31_g3.log`, `p31_g3_results.csv`) | B16 (unset with thresholds) and B25 (R23, P6 off, dropout): 20 runs, 0 failed |
| P6-on recheck vs 37bb518 (`p31_onpath.log`) | R body 15/15, cpp 15/15, default 12/15. The 3 that differ are the dropout-on frames, now admitted and run by the exact engine (by design, see Section 5) |

Suite reconciliation: 1,305 (EDGEFIX) + 142 (`test-engine-equivalence.R`, 198 to 340) + 6 (`-cpp-vsref`, 228 to 234) + 10
(`-cpp`, 148 to 158) + 52 (`-exact`, 230 to 282) + 1 (`test-chronology-amd8.R`, 39 to 40) = 1,516. One skip added (the
P3-1 long-tier test).

## 1. What changed

- `R/simulate_chunk_cpp.R` predicate: item 4 admits TRUE or FALSE; 6 TRUE or FALSE; 7 finite, or Inf when every arm's cap
  is above 0; 8 NULL or one number; 10 the R body's fallbacks; 13c NULL or one finite number. Labels renamed to say what
  each requires. Packing: explicit `rebalance` and `futility_interim` flags (the threshold is NA and unread when the flag
  is FALSE), final thresholds resolved per arm with the R body's expression (`Inf` / `-Inf` when missing).
- Both C++ engines: the amd8 chronology (cutoff with the clock set to the cutoff, last-slot rule, phantom vs-ref rule,
  j = 0 look only before T_final, post-enrolment looks from the current time strictly before T_final, none under Inf),
  the two invariants with the R body's message text, dropout kept (exact engine: one extra uniform per patient, after the
  event time), no rebalancing, `p = 0` for a non-finite final threshold (phase 1 still draws; the exact engine makes no
  posterior call), no interim futility stop with the draws kept in phase 1. Hard stops for items 4, 6, 7 removed.
- After Codex: both slicers and both rebalancing counts stop when calendar time minus enrolment time is NaN (Section 6).
- `R/simulation_driver.R`: the two stale comments only (lines 584 and 593, see Section 4).
- DESCRIPTION 0.0.0.9004, NEWS section.
- Tests: `test-engine-equivalence.R` (P3-1 frames, routine and long tiers; every-option-unset call; F < 0 error in all
  engines; still-excluded items fall back; 6 golden entries), `test-simulate-chunk-exact.R` (replay = R body on 8 P6-off
  and P6-on frames x 2 scenarios; overflow error), predicate tests updated, `helper-chronology.R` and
  `test-chronology-amd8.R` pinned to the R body (Section 4).

## 2. Gate B detail

Frames are built from the R01, R03 and R23 fixtures (`p31_frames.R`). Branch hits per frame are trials (out of the frame's
trials) in which the instrumented R body (`p31_instrument.R`, a deparsed copy of `simulate_chunk_r` with counters,
checked identical to the R body on every chunk) took the branch. Full table with all counters: `p31_gateB_table.md`.

| Frame | Chunks | Trials | Error chunks | Mismatch | Target branch hits (trials) |
|---|---|---|---|---|---|
| B01 R01 P6 off (F = 6) | 142 | 560 | 0 | 0 | j = 0 look 560, post-enrolment looks 233 |
| B02 R01 off, cutoff 10 | 142 | 560 | 0 | 0 | binding cutoff 250, j = 0 look 310 |
| B03 R01 off, F = 0 | 142 | 560 | 0 | 0 | last-slot rule 238 |
| B04 R01 off, F = -1 | 142 | 560 | 74 | 0 | invariant error in 74 chunks (all engines, same message and stream) |
| B05 R01 off, Inf | 142 | 560 | 0 | 0 | Inf T_final 560, Inf final time 237 |
| B06 R01 P6 on, Inf | 142 | 560 | 0 | 0 | Inf t_star 560, Inf final time 230 |
| B07 / B08 R01 dropout, off / on | 142 each | 560 each | 0 | 0 | dropout censors an event 496 / 493 |
| B09 R01 off, rebalancing NULL | 142 | 560 | 0 | 0 | rebalanced 0; companion B09c (threshold 8) rebalanced 258 |
| B10 R01 off, success threshold NULL | 142 | 560 | 0 | 0 | final Inf fallback 227, interim from null_median_arms 560 |
| B11 R01 off, futility threshold NULL (rerun) | 142 | 560 | 0 | 0 | final -Inf fallback 85 |
| B12 R01 off, fallback chain | 142 | 560 | 0 | 0 | interim futility from success vector 560, Inf 259, -Inf 259 |
| B13 / B14 R01 futility prob NULL, off / on | 142 each | 560 each | 0 | 0 | interim draws with no futility rule 560 / 560 |
| B15 R01 every option unset | 142 | 560 | 0 | 0 | last-slot 560, Inf 560, dropout 471, both fallbacks 560, no futility rule 560 |
| B16 R01 unset except thresholds | 142 | 560 | 0 | 0 | Inf 560, dropout 129, last-slot 6 |
| B17 R03 off | 142 | 560 | 0 | 0 | j = 0 look 560 |
| B18 R03 off, F = 0, cutoff 12 | 142 | 560 | 0 | 0 | cutoff 328, last-slot 46 |
| B19 R03 unset except thresholds | 142 | 560 | 0 | 0 | last-slot 300, Inf 560, dropout 418 |
| B20 R23 off | 142 | 312 | 0 | 0 | phantom vs-ref exit 234, look at a phantom arrival 39, j = 0 look 312 |
| B21 R23 off, cutoff at half accrual | 142 | 312 | 0 | 0 | cutoff 202 |
| B22 R23 off, F = 0 | 142 | 312 | 0 | 0 | phantom look suppressed 86, last-slot 94 |
| B23 / B24 R23 Inf, off / on | 142 each | 312 each | 0 | 0 | Inf T_final 312 / Inf t_star 312 |
| B25 / B26 R23 dropout, off / on | 142 each | 312 each | 0 | 0 | dropout censors an event 312 / 312 |
| B27 R23 off, rebalancing NULL | 142 | 312 | 0 | 0 | rebalanced 0 (B20, same seeds with threshold 50: 179) |
| B28 R23 off, F = -1 | 142 | 312 | 142 | 0 | invariant error in every chunk |
| B29 R23 every option unset (vs-ref) | 142 | 312 | 0 | 0 | phantom suppressed 91, Inf 312, dropout 301 |

Two frames were changed because their target branch was hit 0 times in the first run. B11 (final futility threshold
NULL) hit the `-Inf` fallback 0 times at the captured interim futility probability 0.346, because the looks stop every
trial that would reach the final futility test; it now uses 0.9 (85 hits). B09 cannot show "no rebalancing" on R01,
which never reaches 50 events; B09c (threshold 8, rebalanced in 258 trials) is its companion, and B27 against B20 shows
the same on R23. The first invariant (enrolment, stop or final time above a finite cutoff) cannot fire by construction
and has no frame; the second fires in B04 and B28 and in a test.

No branch had to stay excluded for RNG reasons: every ported branch is bit-identical including `.Random.seed`.

## 3. Gate S detail

Exact engine (default selection) against fresh R-body runs (`use_cpp_engine = FALSE`, 4 FORK workers), both with
`return_variance = TRUE`, seeds 3120xxx (R body) and 3121xxx (exact), none used elsewhere in dev/perf. Test as phase 2 G5
(`p31_gateS_compare.R`): two one-sided tests at 0.05 each, margins 0.01 (power), 0.005 (type I error), 1% of the R-body
mean (EN), plus the 0.80 / 0.10 classification. "Test power" is the probability that both one-sided tests reject at a true
difference of 0, at this run's standard errors.

| Frame | R body / exact trials | Power diff (test power) | Type I diff (test power) | EN0 diff (test power) | EN1 diff (test power) | Result |
|---|---|---|---|---|---|---|
| B02 R01 off, cutoff 10 | 40,000 / 400,000 | -0.0022 (1.00) | +0.0003 (0.99) | -0.008 (1.00) | -0.011 (1.00) | PASS |
| B07 R01 off, dropout | 20,000 / 200,000 | -0.0005 (0.93) | +0.0002 (1.00) | +0.009 (0.71) | +0.000 (1.00) | PASS |
| B15 R01 every option unset | 10,000 / 100,000 | 0 (all inconclusive) | 0 | 0 (N = 32) | 0 | PASS (degenerate) |
| B16 R01 unset except thresholds | 100,000 / 800,000 | +0.0000 (1.00) | -0.0001 (1.00) | +0.002 (1.00) | -0.002 (0.72) | PASS |
| B20 R23 off (vs-ref) | 80,000 / 400,000 | +0.0023 (1.00) | +0.0007 (1.00) | +0.138 (1.00) | +0.065 (1.00) | PASS |

Each frame is run under its null (scenario 1) and alternative (scenario 2). Wall time 2.7 h (09:38 to 12:19), nice'd, on a
machine running about 20 other jobs at nice 0.

The first readout (`p31_gateS_report_first.tsv`; R body 20,000 trials on B02, B16, B20's 40,000, exact 200,000) passed 2
of 5: B02 power (diff -0.006, z_lo 1.4), B16 EN1 (diff +0.049 against a margin of 0.049, test power 0.00) and B20 EN0
(diff +0.18, z_hi 1.6). I then added fresh seeds to those three frames (`p31_gateS_ext.sh`, `runs_ext.txt`) and reran the
same comparison on the pooled runs. This is an extension decided after seeing the result, not a pre-set sample size, so
the pass carries that caveat. B20 EN0 is still +0.14 patients (0.46%, about 2.7 standard errors) with the exact engine
higher, inside the 1% margin. It has the same sign in both readouts, so it is likely a real small difference between the
exact posterior and 800 Monte Carlo draws, of the kind phase 2 recorded on R23.

Gate S was run before the Codex P2 fix. It does not need rerunning: the guard fires only when a calendar time minus an
enrolment time is NaN, which needs both to be infinite. The difference of two finite doubles is never NaN, and gate S
frames have finite arrival times (a rate of 2 per month). The 12 goldens, layer 1 and the gate B subset are unchanged
bit for bit after the fix.

### 3a. Pre-set confirmation (p31c, PI ruling 2026-11-04)

Sizes fixed at the extended readout's final sizes, new seeds 3140xxx (R body) and 3141xxx (exact), installed 0.0.0.9004
(`p31c_gateS.sh`, `p31c_gateS_report.tsv`). Result 4 of 5 frames. B02, B07, B15 and B20 pass on every statistic (B20 EN0
now z 7.1 and 4.5). B16 fails on EN under the alternative only: R body 4.981, exact 4.953, diff -0.028 against a margin
of 0.050, se 0.018, z_lo 1.19, test power 0.72 at this size. The first readout's diff on the same statistic was +0.049,
the opposite sign. Per the PI's rule the 0.0.0.9003 backup was restored into `~/R/library` and P3-1 reopened on this
one frame. PI ruling: rerun B16 alone, sized in advance for test power at least 0.95 (`p31d_gateS_B16.sh`: R body
200,000, exact 1,600,000, se about 0.013, test power about 0.97, seeds 3150xxx and 3151xxx, library built from 325c664).
Reinstall 0.0.0.9004 only if it passes.

### 3b. B16 powered confirmation (p31d)

`p31d_gateS_B16.sh`, sized before running: R body 200,000, exact 1,600,000, seeds 3150xxx and 3151xxx, library built
from 325c664. PASS on every statistic. EN under the alternative: R body 4.932, exact 4.946, diff +0.014 against a margin
of 0.049, se 0.0129 (as sized), z 4.9 and 2.7, test power 0.97. Power 0.0335 vs 0.0337, type I 0.0010 vs 0.0010, EN0
4.123 vs 4.124. The three B16 readouts on this statistic gave +0.049, -0.028 and +0.014, so no systematic difference.
Gate S closed: B02, B07, B15 and B20 from the pre-set confirmation (3a), B16 from this run. 0.0.0.9004 reinstalled into
`~/R/library` 2026-11-04 (backup of 0.0.0.9003 kept).

## 4. Mismatches with the brief or plan

1. The stale comments P3-0 placed at merged lines 579 and 588 are at lines 584 and 593 of c936bde. Both
   fixed in place (still lines 584 and 593). A third stale comment, line 268 ("FALSE = frozen 422eb40"), and the roxygen at lines 210-212 ("The supported set
   requires `looks_stop_at_final = TRUE`", also in `man/run_simulation_pure.Rd`) are now wrong but outside what I could
   edit. Open item.
2. Interim threshold fallback (IL:59-61): `vec[[arm]] %||% vec[1]` errors in `[[` for any non-NULL vector that lacks the
   arm name, so the `vec[1]` fallback is unreachable. Only NULL vectors fall back at the interim. Unnamed or missing-arm
   vectors are admitted only at the final, where the R body falls back to `Inf` / `-Inf`. A NULL interim futility vector
   is admitted only with item 13c NULL (with a futility probability set the R body errors at the first look).
3. Item 7 needed one extra condition. With `max_trial_time = Inf` and a zero cap, the person-time denominator is
   `0 * Inf = NaN` and the R body errors at the first look, while the C++ read it as 0. Such frames stay on the R body.
4. There is no `options(evolveTrial.engine = "r")`. The R body is `options(evolveTrial.use_cpp_engine = FALSE)`, and
   `engine = "cpp"` is bit-identical to it (NEWS says so).
5. P6-on 15/15 against 37bb518 holds for the R body and cpp, but the default engine gives 12/15. The three dropout-on
   frames used to fall back to the R body and are now admitted, so the exact engine runs them, by design.
6. `test-chronology-amd8.R` and `helper-chronology.R` (byte-identical to main until now) trace R-body functions and pin
   R-body values. Once the frames were admitted they ran the exact engine and failed 21 times. They now set
   `use_cpp_engine = FALSE`, and one test also checks that engine "cpp" gives the pinned table.
7. Items 3 and V4 (unnamed caps) end to end make the R body itself error in `make_state()`, so their fallback is checked
   at frame level (as before) and not end to end. Item 12 and V3 are frame-level as well.

## 5. Behavior changes for users

Calls that leave these options at their defaults now run the exact engine, including the paper wrappers on futility-only
designs (`looks_stop_at_final` unset). P6-on frames with dropout on, Inf, no rebalancing, a threshold fallback or no
interim futility rule also move from the R body to the exact engine. Frames admitted in 0.0.0.9003 are unchanged (goldens
6/6, P6-on R body and cpp 15/15).

## 6. Codex review (`p31_codex_review.md`) and how each finding was closed

- **P1:** none.
- **P2, overflow to Inf under `max_trial_time = Inf`.** Reproduced (`p31_codex_p2_repro.R`, `_before.log`): R errors
  "missing value where TRUE/FALSE needed", while cpp and exact returned `Exp_Time = Inf`. Fix: `follow_up_since()` in both
  engines, used by both slicers and both rebalancing counts, stops with the same message when calendar time minus
  enrolment time is NaN. That is where R stops: its slice keeps NaN rows, and `calculate_interval_metrics_fast()` or the
  rebalancing `if()` then fails. After the fix all engines error (`_after.log`). The exact engine has its own streams, so on
  a 5-trial seed it can draw no overflowing trial. The test uses 60 trials. The guard cannot fire on a finite path, and
  goldens 12/12, layer 1 and gate B (B01, B02, B05, B20, B23) are unchanged. New test in `test-simulate-chunk-exact.R`.
- **P3, no deterministic test of the exact engine's P6-off chronology.** New test: the exact file's replay instantiation
  equals the R body itself (result and stream) on R01 off, cutoff 10, F = 0, Inf + dropout + no rebalancing, R01 P6 on +
  Inf + fallbacks, R23 off, R23 cutoff, R23 F = 0 + dropout, both scenarios. Mutation (`p31_mutation_j0.sh`, scratch build
  `dev/perf/lib/p31_mut_j0`, tree untouched): with the exact engine's j = 0 look deleted, the new test fails 12 times and
  nothing else fails.

## 7. Open items

1. Closed by the coordinator: roxygen, Rd and the line-268 comment in `R/simulation_driver.R` updated (documentation and
   comments only, `tools::checkRd` clean).
2. Gate S passed only after a post-hoc sample extension (Section 3). PI ruling 2026-11-04: commit and install now, then
   run a confirmation with the extended sizes fixed in advance on new seeds. If it fails, revert the install and reopen P3-1.
3. B20 EN0 +0.46% (exact above the R body), stable across both readouts. It is within the margin and worth a note for
   P3-2 on R23-type designs.
4. Codex review R of the commits themselves is the coordinator's step.
5. Scratch libraries and builds `dev/perf/lib/p31_{base,dev,final,mut_j0}` and `dev/perf/build/p31_*` can be deleted.
