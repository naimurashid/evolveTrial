# P3-0: merge main (4dc2f1f, amd8) into the engine branch

Date 2026-11-04. Worktree `~/research/evolveTrial-merge`, branch `perf/merge-main`, HEAD dfb4527. State: `git merge
--no-commit --no-ff main`, all conflicts resolved and staged, nothing committed, pushed or installed into `~/R/library`.
`dev/perf/phase3/` (this report, scripts, logs, RDS) is untracked. `~/research/evolveTrial-perf` and
`~/research/evolveTrial` were only read (the build script was copied from the former).

## Counts

| Item | Result |
|---|---|
| Full suite, merged tree (`p30_pkg_tests.log`) | 1234 passed, 0 failed, 35 skipped, 0 errors, 112 s |
| Baseline perf/engine-cpp (37bb518) | 1186 passed, 0 failed, 35 skipped |
| main's own suite, scratch checkout of 4dc2f1f (`p30_main_tests.log`) | 353 passed, 0 failed, 33 skipped. Every file shared with the engine branch has the same counts in both; main's only new file is `test-chronology-amd8.R` (39) |
| Reconciliation | 1186 + 39 (main's chronology tests) + 9 (new P3-0 test) = 1234. Skips unchanged at 35 |
| `test-engine-equivalence.R` alone | 182 passed, 0 failed, 2 skipped (173 before + 9 new expectations, the 2 skips are the long tier) |
| `test-chronology-amd8.R` + `helper-chronology.R` | 39 passed, 0 failed, 0 skipped (alone and in the full suite). Both files byte-identical to main |
| Same chronology file with `options(evolveTrial.looks_stop_at_final = TRUE)` | 28 passed, 11 failed. So main's tests only make sense with P6 off, and under the package default they run with P6 off (no test file leaks the option) |
| P6-off R body, merge vs main, 15 frames (`p30_offpath_vs_main.R`) | 15 of 15 `identical()` (table and `.Random.seed`) |
| P6-on, merge vs 37bb518, same 15 frame shapes (`p30_onpath_vs_37.R`) | R body 15/15, engine "cpp" 15/15, default (exact) 15/15 identical, merged R body vs merged cpp 15/15 |
| Layer 1 harness, 10 designs (`p30_layer1.log`) | 2,120 chunks, 0 mismatches |

Frames for the two behavioral checks: R01, R03, R23 fixtures, each as captured, with a binding cutoff
(0.5 N / accrual rate), F = 0, F = 0 with a binding cutoff, and `max_trial_time = Inf` (P6-off) or
`random_dropout = TRUE` (P6-on). Seed 2026, `return_variance = TRUE`, 60 trials (15 for R23).

## 1. P6-on path is the 37bb518 code

Textual check (`cmp`) of `R/simulation_driver.R`, 37bb518 vs merged:

| 37bb518 lines | merged lines | Content | Result |
|---|---|---|---|
| 1-89 | 1-89 | header, roxygen before the conflict | identical |
| 252-562 | 263-573 | signature tail, checks, setup, accrual loop start | identical |
| 563-588 | 574-599 | P6 cutoff block and P6 phantom vs-ref block | identical |
| 654-697 | 701-744 | enrolment of the arrival | identical |
| 714-733 | 775-794 | P6 post-enrolment `else` branch (`t_star` grid, "end D1 option") | identical |
| 735-748 | 796-809 | final-time block | identical (37bb518 text kept for both paths) |
| 749-1360 | 828-1439 | final analysis, accumulation, engine dispatch | identical |

Differences from 37bb518, all off the P6-on path or no-ops on it: roxygen (90-98 to 90-123), one comment at
249-250, the amd8 block inserted after the P6 blocks and wrapped in `if (!looks_stop_at_final)` (600-633), the
off-branch body at 745-773 (main's code), the invariants wrapped in `if (!looks_stop_at_final)` (810-827), and
main's Windows PSOCK fallback in `run_scenarios()` (1447-1462). One shared line changes: the accrual-phase look
dispatch at 667-669 gains main's `if (!interim_ok_here) {} else` prefix. `interim_ok_here` is set `TRUE` before the
amd8 block and only that block can set it `FALSE`, so on P6-on the prefix is a no-op (the 15/15 identity above and
the 2,120-chunk layer 1 run confirm it).

Two stale comments sit inside the P6 blocks and were left as they are to keep the blocks textually identical: "FALSE
keeps 422eb40 (no time cap on accrual)" (merged 579) and "FALSE keeps the 422eb40 behavior (look first, then break)"
(merged 588). After the merge FALSE is amd8, not 422eb40. Coordinator's call whether to reword them in a later commit.

amd8 invariants (main 704-718) on the P6 path: by construction they hold there (P6 enrols nobody at or after the
cutoff, looks run before T* <= cutoff, final time <= cutoff, and a trial with every arm stopped has its last stop at or
after every enrolment, and vs-ref trials never take the all-stopped branch). So adding them would not change results on
frames where they pass. This was not verified by running them on the P6 path, and they were not added there.

## 2. P6-off path is main

The off branch is main's code verbatim (cutoff closes accrual with `current_time <- max_trial_time`, the
`interim_ok_here` rules (a) and (b), the look at the end of enrolment only if `current_time < t_final_planned`,
post-enrolment looks strictly before `t_final_planned`, invariants). The final-time text is 37bb518's, which
evaluates the same expression as main's `t_final_planned`. The 15/15 bitwise identity with a build of 4dc2f1f
(including binding cutoffs, F = 0 and vs-ref) and the 39/39 chronology tests confirm it.

## 3. random_dropout and amd8

main has no `random_dropout`. It always draws `runif(1, 0, censor_max_time_sim)` and uses it as the dropout time,
which is the 422eb40 censoring model. On the engine branch `random_dropout = TRUE` (default) makes the same draw and
keeps it; `FALSE` makes the same draw and then sets the time to `Inf` (RNG stream unchanged). So main's semantics
assumed the 422eb40 censoring model, and the merged default reproduces it exactly (the P6-off identity runs drop
`random_dropout` from the fixtures, so the merge uses its default TRUE). No conflict. Neither default changed.
`random_dropout = FALSE` with P6 off is a new combination (amd8 chronology, no dropout) that main never offered.
It runs the R body.

## 4. Supported set

`R/simulate_chunk_cpp.R` is unchanged by the merge. Items 4 (`looks_stop_at_final` TRUE) and 6 (`random_dropout`
FALSE) still gate admission, so no amd8 frame reaches either fast engine. New test at the end of
`tests/testthat/test-engine-equivalence.R`: "P6-off frame with a binding max_trial_time falls back and equals the R
body". For R01 (cutoff 10 months) and R23 (cutoff at half the accrual time), option unset vs
`options(evolveTrial.use_cpp_engine = FALSE)`: identical table and RNG, no notice, the cutoff binds (expected N
below each arm's maximum), and an explicit `engine = "exact"` errors "outside the supported set". 9 expectations.

## 5. Version, NEWS, man

DESCRIPTION `Version: 0.0.0.9003`. main's auto-merged DESCRIPTION edits came along (RoxygenNote replaced by
`Config/roxygen2/version: 8.0.0`, `SystemRequirements`). NEWS: new 0.0.0.9003 section (merge, amd8 on P6-off, P6-on
unchanged, exact default unchanged, random_dropout default), then "0.0.0.9002 (engine branch: exact engine as
default)" and "0.0.0.9002 (main: amd8 chronology, BATON-C freeze Amendment 8)", both bodies as written. The engine
section still says "a call with every option at its default still runs the R body, unchanged". That was true when
written and the 0.0.0.9003 section says what changed. `man/run_simulation_pure.Rd` edited by hand to match the
roxygen (main's `min_follow_up_at_final` and `max_trial_time` text plus one scoping sentence; the `looks_stop_at_final`
FALSE clause now describes amd8 instead of 422eb40). No roxygen run. `tools::checkRd` clean.

Found, not fixed (pre-existing at 37bb518, not in this brief): `tools::codoc` reports `random_dropout` in the code
but not in `run_simulation_pure.Rd`.

## Build and commands

`dev/perf/phase3/p30_install.sh merge` (copy of the phase 2 `ds_install.sh`, pointed at this worktree) builds
`dev/perf/lib/p30_merge`. Suite: `Rscript dev/perf/harness/run_pkg_tests.R dev/perf/lib/p30_merge <out.rds> [filter]`.
Scratch libs `dev/perf/lib/p30_main` (git archive of main) and `dev/perf/lib/p30_37` (git archive of 37bb518), build
trees under `dev/perf/build/p30_*`, all deletable. Layer 1: `Rscript dev/perf/harness/run_layer1.R --lib
dev/perf/lib/p30_merge --no-cache`.

## Layer 1

`Rscript dev/perf/harness/run_layer1.R --lib dev/perf/lib/p30_merge --no-cache` (fresh R-body references, no cache):
2,120 chunks compared, 0 mismatches (0 with result differing, 0 with stream differing), 304 s. Same count as the G4
readout. Output `dev/perf/harness/out/layer1_20261104_070807.csv`, log `dev/perf/phase3/p30_layer1.log`.
