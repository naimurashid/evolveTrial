# P2-3 report: local gate harness for G5 and G6

Date 2026-11-04. Worktree `~/research/evolveTrial-perf`, branch `perf/engine-cpp`, HEAD 0854a58. Nothing committed, no
package source touched, no other repo written, no Longleaf. `dev/perf/lib/dev_p2` did not exist, so p2 mode itself is
**untested** (see Section 4).

## 1. Files

In `dev/perf/phase2/`:
- `p2_cert.R`: derived from `layer2/layer2_cert.R`, modes `ref`, `dev`, `p2`.
- `p2_gate.sh`: runs the 75 reference rows in p2 mode.
- `p2_compare.R`: G5 comparator, writes `g5_report.tsv` and `g5_report.md`, exit 0 PASS, 1 FAIL, 2 input error.
- `p2_selftest.sh`: the comparator self-tests (below).
- `p2_timing.sh` and `p2_t7.R`: G6 timing (adapted from `step7_timing/scripts/run.sh` and `t7.R`).

Under `dev/perf/evidence/phase2_gate/`: `selftest/` (self-test reports and inputs). `rows/`, `run.log` and `timing/` are
created by the gate and timing scripts when they run.

## 2. Changes from layer2_cert.R (all tagged "P2:" in the file)

- P2-a: usage `dev|ref|p2`; modes `ref` and `dev` unchanged (checked, Section 3).
- P2-b: p2 library default `lib/dev_p2` (env `P2_LIB`), `options(evolveTrial.engine = "exact")`, asserted after load.
- P2-c: identity. Worktree HEAD SHA, whether R/, src/, DESCRIPTION differ from HEAD, and env `P2_ENGINE_SHA` are recorded
  (`worktree_sha`, `worktree_dirty`, `p2_engine_sha_env`), not required. `L2_ENGINE_SHA` is not read in p2. Required:
  `P6RecertEvolveTrialSHA` = 7565a17, no `VerifyM4PatchMD5`, both `.simulate_chunk_cpp` and `.simulate_chunk_exact` in the
  namespace. Runner/staged md5 checks unchanged.
- P2-d: per-process probe files (workers are FORKs) for both entries, one line per call: pid, trials = `length` of the first
  argument, cfg `random_dropout`, `min_follow_up_at_final`, `looks_stop_at_final`, `mode`. Wrappers use `function(...)`, so
  they do not depend on the exact formals; the formals are recorded (`exact_formals`).
- P2-e: p2 engine check = zero interim_check probe lines, zero P6 frame events, zero trials through `.simulate_chunk_cpp`,
  exact trials summed over all processes = (number of scenarios) x R, cfg FALSE / 6 / TRUE; the run_scenarios-level
  terms are as in dev. The echo string carries `exact[lines pids trials expect ...]`.
- P2-f: the rds adds `EN_null_var`, `EN_alt_var`, `power_var`, `type1_var`, `var_attr`. These come from the `variance`
  attribute that `run_scenarios(return_variance = TRUE)` leaves on the wrapper result, and they are the **variance of the
  mean**. Per-trial variance = value x R. Checked: R01 seed 900001, value x R is 50.8, 49.4, 34.9, 31.5, 30.6 at R = 100,
  200, 500, 1000, 4000 (heavy-tailed, so noisy at small R, converging), so the scale is per-trial/R.
- Layer 2 header comments kept below a P2 header block.

## 3. Tests run

p2_cert.R, R = 200, 19 workers, seed 900001, R01 and R23, modes ref and dev: all eight capture checks TRUE, engine_ok
TRUE. dev reproduced ref exactly at R = 200 (R01 power 0.88 / type1 0.025 / EN 11.885 / 29.400; R23 0.76 / 0.055 /
59.25 / 88.875). R01 at R = 10,000 in ref mode reproduced reference row 900001 (8245, 150, 11.247, 28.291) in 24 s.

p2_compare.R self-tests (`p2_selftest.sh`, rerunnable, exit 0 only if all hold; each also checks the verdict line):
1. reference vs itself: PASS, 10 of 10.
2. phase 1 layer 2 dev outputs (`evidence/layer2_l2a/compare.csv`) vs reference: PASS, 10 of 10.
3. power +0.02 on R23 (counts raised by 0.02 R in each seed): FAIL, 9 of 10, R23 power_ok FALSE.
4. type I +0.01 on R01: FAIL, 9 of 10, R01 type1_ok FALSE.
5. extra: EN_null +2% on R02: FAIL, 9 of 10.
6. rds input path: five real ref-mode rows of R01 (R = 10,000, seeds 900001 to 900005), relabelled mode p2, with --subset:
   SUBSET-PASS (1 of 1; report in `selftest/rds_path_R01/`). The unrelabelled ref-mode rds fail on `rows_ok`, as intended.
   One seed alone fails EN_null at zero difference (z = 1.4 < 1.645), see 5.1.

Variance for tests 1 to 5 (CSV inputs carry none) comes from `selftest/var.csv`: one ref-mode R = 2000 run per design
(first seed), per-trial variance = variance attribute x R.

G5 output pooling: counts summed over seeds, EN means averaged over seeds (equal R), EN variance by the law of total
variance (mean within-seed per-trial variance plus variance of the seed means, divisor S).

## 4. Not tested

- p2 mode of p2_cert.R, `p2_gate.sh` end to end (only `--dry-run` with `P2_LIB` pointing at dev_l2: lists 75 rows), and
  `p2_timing.sh` (syntax checked only; it refuses to start without dev_p2 or above load 1.0). p2 mode without a library stops
  at the identity check. Run `p2_cert.R R01 900001 200 p2 out.rds` as soon as dev_p2 exists, before the 75-row run.
- Assumptions about the P2-2 interface that the first p2 run will confirm: the trial index vector is the first argument of
  `.simulate_chunk_exact`, cfg is the second positional argument or named `cfg` and has fields `random_dropout`,
  `min_follow_up_at_final`, `looks_stop_at_final`, `mode`; the driver looks the entry up in the namespace at call time (as
  phase 1 does), so `assignInNamespace` wrapping reaches the FORK workers; trials arrive in blocks that sum to R.

## 5. Ambiguities and mismatches in the plan (reported, not resolved)

1. EN test power. With a 1% margin, at zero true difference the pooled z for the lower and upper test is about 3.2 (R01,
   5 seeds) and 1.4 for a single seed. A true EN difference of about half the margin gives z near the critical 1.645, so
   a pass certifies only differences below roughly 0.5 to 0.7 of the margin. Designs with 5 seeds and a large per-trial
   variance (R23, R26) have the weakest margin. Not changed.
2. Variance source: "the phase 2 run's variance applied to both samples" is implemented as the pooled per-trial variance
   above. The plan does not say how to pool over seeds.
3. Classification ties: "same side of 0.80 / 0.10" is implemented as power >= 0.80 and type I <= 0.10 (the constraint
   directions). The report lists designs whose pooled reference value is within 3 SE of a threshold; none is known yet.
4. G6 wording: "R23 under 8 s per 1,000 trials at one worker on one worker" is duplicated. The step 7 scripts time one
   `run_scenarios` call with 1,000 trials per scenario and two scenarios, so "per 1,000 trials" is read as elapsed seconds per
   call (column `per_1000_s` = elapsed). The 75.1 s baseline in the plan is not in `timing_step7.csv` (R23 ref at 1 worker
   there: 235.8 s, dev_l2 25.4 s, dev_l2 16 workers 3.4 s), so which quantity 75.1 s measures should be confirmed before
   the 8 s target is read against these numbers. Efficiency is computed as T(1) / (16 T(16)) from medians.
5. Plan Section 3 names the reference rows as the comparison for both engines, but those rows store no variance and no
   per-seed trial vectors, so a variance-aware comparison of the reference side was impossible, as the plan anticipated.
6. The step 7 log started at load 1.25, above this script's default limit of 1.0 (P2_MAX_LOAD); timing needs a quiet machine,
   so run it when no gate job is running.
