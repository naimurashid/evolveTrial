# P2-2 report: phase 2 chunk engine, counter-based generator and driver

Date 2026-11-04. Worktree `~/research/evolveTrial-perf`, branch `perf/engine-cpp`, HEAD 106b79c, nothing committed or
pushed, installed only to `dev/perf/lib/dev_p2` (plus scratch libraries under `/tmp/claude-1000/p2build`). Plan
`dev/plans/2026-11-04-phase2-plan.md` Section 1 items 2 and 3, gates G3 and G4 (Section 3).

**Result.** G3 PASS (100 of 100 runs, 10 designs). G4 PASS: package tests 831 passed, 0 failed, 33 skipped (783 + 48 new);
phase 1 layer 1 on `dev_p2`: 2,120 chunks, 0 mismatches. Driver-level check of the default mode against `dev_l2`: identical tables and
`.Random.seed` (12 of 12). Replay identity (the copied trial logic against the phase 1 engine): identical on all 20
scenario frames of the 10 designs. Sanity R01, R23 at R = 2,000: max |z| 1.61. `p2_cert.R` p2 engine check: R01 and R23
`engine_ok TRUE`. Timing at one worker: R01 8.82 s (phase 1, `dev_l2`) vs 0.079 s (phase 2) per 1,000 trials, R23 9.97 s vs 0.109 s (about 110x and 90x). Profile: look slicing is 6.2% of R01 trial time and 75.5% of R23
trial time, so the plan's P2-4 trigger (30% of R01 or R23) is met by R23.

## 1. Design

### 1.1 Engine (`src/simulate_chunk_exact.cpp`)

The phase 1 trial loop (`src/simulate_chunk.cpp:595-777`) is copied line for line into one template
`run_trial<Rng, Post>`, together with copies of `Cfg`/`unpack` (same packed list, same checks, messages prefixed
`simulate_chunk_exact`), `TrialState`, `slice_arm`, `gates_pass_vsref`, `trt_can_enroll` and both interim looks.
`src/simulate_chunk.cpp` is not touched. Only two things vary, as policies:

- `Rng`: arrival gap, arm choice, event time, dropout time.
- `Post`: single-arm interim futility probability, final efficacy and futility probabilities, and the PH step plus
  `Pr(logHR >= m)` and `Pr(logHR < m)`.

Two instantiations are compiled.

| Entry | Rng | Post | Use |
|---|---|---|---|
| `.simulate_chunk_exact(sim_indices, cfg, scenario_seed)`, `rng = false` | `CounterRng` | `ExactPost` | phase 2 engine |
| `.simulate_chunk_exact_replay(sim_indices, cfg)` | `RStreamRng` (R's stream, phase 1 order, dropout uniform drawn and discarded) | `McPost` (phase 1 Monte Carlo: `draw_posterior_hazard_samples_cpp` + `calculate_median_survival_matrix_cpp`; a copy of `ph_sample_core` without `keep_medians`) | test hook |

The replay instantiation must return exactly what `.simulate_chunk_cpp` returns on the same `.Random.seed` and leave the
same `.Random.seed`. It does, on both scenarios of all 10 designs (R = 200; R01 and R23 also in the package tests at
R = 120). Since the exact engine runs the same template code, this checks that the copied trial logic (arrival loop,
phantom-arrival rules, caps, gates, look schedule, rebalancing, slicing, final analysis, error paths that a run reaches)
is the phase 1 logic. What the replay cannot check is the two exact policies themselves; they are short and listed
here.

`ExactPost` (P2-1 functions, `src/posterior_exact.h`):
- HC interim futility: `median_probs(E, PT, alpha, beta, lens, K, futility_thr_interim).below` (phase 1 `mean(med < thr)`).
- HC final: efficacy `median_probs(..., success_thr).above`, futility `median_probs(..., futility_thr).below`, each a
  separate exact evaluation on the same slice and the current (rebalanced) `lens`; futility only when efficacy fails, as
  in phase 1. Status 2 (invalid input) gives NaN, which the unchanged phase 1 code turns into an error at the interim and
  into "not passed" at the final, as phase 1 does with an NA tally. Status 1 (rate replaced by 1e-6) emits one R
  warning per call (unreachable in the supported set, P2-1 item 4).
- PH: `ph_beta_mode_var_cpp` exactly as phase 1 calls it, then `loghr_probs(bmean, bsd, m)`: interim futility `.ge` at
  `log_margin_interim`, final efficacy `.lt` at `log_eff_final`, final futility `.ge` at `log_margin_final`. Nothing is
  drawn; the draw counts in cfg are unpacked and checked but not used.

`CounterRng` (header comment of the file):
- `seed64 = (uint32 s1 << 32) | uint32 s2` from the two-integer scenario seed.
- Per trial: `key = mix64(seed64 ^ mix64(trial_index + 0x9E3779B97F4A7C15))` (`mix64` = splitmix64 finalizer, a
  bijection, so distinct indices give distinct keys under one seed); xoshiro256++ state = four successive splitmix64
  outputs from `key` (the seeding its authors recommend).
- `U = ((x >> 12) + 0.5) * 2^-52`, in `[2^-53, 1 - 2^-53]`. (The common `(x >> 11) * 2^-53` form plus a half would
  round the top value to exactly 1; the 52-bit form cannot.)
- Gap `-log(U) / rate` (as `-log(U) * (1 / rate)`, the scale phase 1 passes to `R::rexp`); event time
  `scale * (-log U)^(1 / shape)` (R's `rweibull` formula); arm by phase 1's `ProbSampleNoReplace` port with `rT = U`, so
  the same probabilities, and no draw when one arm is eligible; dropout draws nothing and stores `Inf` (phase 1 stores
  `Inf` after discarding its uniform). Order per arrival: gap, look, arm, event time.

Return value: `sim_indices` and per-trial outcomes in that order (`final_n`, `final_events` n x A double;
`stop_eff`, `stop_fut`, `final_eff`, `final_fut`, `final_inc` n x A integer; `final_time`). `.exact_reduce(blocks, cfg,
n_sims)` requires each index 1..n exactly once over the blocks, accumulates in trial-index order with phase 1's plain
double `+` (`simulate_chunk.cpp:780-789`), stores raw values of trial i at position i, and returns the phase 1 chunk list
(same 12 names, same types, arm names on the sums).

Other exports, test hooks only: `.exact_rng_uniforms(seed, trial, n)`, `.exact_rng_raw(state, n)`.

### 1.2 Driver (`R/simulate_chunk_exact.R`, `R/simulation_driver.R`)

- `options(evolveTrial.engine)`: `"cpp"` (default) or `"exact"`; any other value stops. `"exact"` takes effect only
  when `evolveTrial.use_cpp_engine` is TRUE; FALSE still selects the R body (tested).
- In `run_simulation_pure`, after the `chunk_probe` seam: in exact mode only, `.chunk_cpp_assert_supported(frame, cfg)`
  (stops outside the phase 1 set, before any draw), then the scenario seed `sample.int(.Machine$integer.max, 2L,
  replace = TRUE)`, one call per scenario dispatch, so `set.seed()` before `run_scenarios` reproduces the run. Nothing
  else touches R's generator in exact mode (the C++ entry is `rng = false`).
- Dispatch: the existing three sites. Sequential (few simulations, or `parallel_replicates = FALSE`): one call over
  `1:n`. Parallel: the cluster is set up exactly as before (`makeCluster` FORK/PSOCK or the external cluster), then
  `parallel::clusterApplyLB(cl, blocks, .exact_block_worker, cfg, seed)` instead of the `sample.int` seed list and
  `parLapply`. Blocks are contiguous index ranges of `getOption("evolveTrial.exact_block_size")` trials, default
  `min(250, max(1, ceiling(n / (8 * workers))))`. `.exact_block_worker` is a namespace function (the cluster ships a
  namespace reference, not the frame) that looks `.simulate_chunk_exact` up at call time, so the harness wrappers reach
  FORK workers.
- `aggregate_results` receives one partial list (from `.exact_reduce`), so it, `run_scenarios` and BATON see the phase 1
  names and meanings.
- Default mode: the code path is the old one with an untaken branch; the only additions it executes are two
  `getOption` reads.

## 2. Files changed

New: `src/simulate_chunk_exact.cpp`, `R/simulate_chunk_exact.R`, `tests/testthat/test-simulate-chunk-exact.R` (48
expectations, about 11 s), `dev/perf/phase2/g3_check.R`, `g3_results.csv`, `p2_2_sanity.R`, `p2_2_sanity.csv`,
`p2_2_sanity.log`, `p2_2_timing.R`, `p2_2_timing.log`, `p2_2_profile.R`, `g4_pkg_tests.log`, `g4_pkg_tests.rds`, `g4_layer1.log`, this
report. Modified: `R/simulation_driver.R` (+22 / -2, the three dispatch sites and the exact preamble), and
`R/RcppExports.R`, `src/RcppExports.cpp` (regenerated by `Rcpp::compileAttributes()`, additions only). Unchanged:
`src/simulate_chunk.cpp`, `R/simulate_chunk_cpp.R`, `src/posterior_exact.*`, `src/Makevars`, every harness file.

`dev/perf/lib/dev_p2` is installed from a copy of the worktree (git-tracked plus the new untracked files, `dev/`
excluded) with one DESCRIPTION line appended, `P6RecertEvolveTrialSHA: 7565a17c…` (HEAD descends from 7565a17 and
58a4850, checked), see 4.1. `EngineCppSHA` is not written, since the build is not a commit.

## 3. Gates and checks

### G3, determinism (`dev/perf/phase2/g3_check.R`, `g3_results.csv`)

Each design's capture, `run_scenarios(..., seed = cert_seed, return_variance = TRUE)` at R = 2,000, engine exact. Runs:
workers 1, 4, 16, 19 (`parallel_replicates = TRUE`) x block sizes 13 and 128, plus one sequential run; each compared
with `identical()` to the first on the returned table and on `.Random.seed`. A wrapper counts trials through
`.simulate_chunk_exact` (all 4,000 per run, every run) and makes `.simulate_chunk_cpp` stop if reached. Then seed + 1 at
19 workers must differ.

| | runs | pass |
|---|---|---|
| same seed (9 per design) | 90 | 90 identical tables and `.Random.seed` |
| seed + 1 (1 per design) | 10 | 10 different tables |

**G3 PASS**, 38 s for all 100 runs.

### G4, no regression

- Package tests (`dev/perf/harness/run_pkg_tests.R dev/perf/lib/dev_p2`): **831 passed, 0 failed, 33 skipped, 0 errors,
  0 warnings** (783 before; the 48 new are `test-simulate-chunk-exact.R`). Log `g4_pkg_tests.log`.
- Layer 1 (`dev/perf/harness/run_layer1.R --lib dev/perf/lib/dev_p2`, defaults): **2,120 chunks compared, 0 mismatches** (0 with the result differing, 0 with the stream differing), all 10 designs, masters 900001 to 900003, 16 and 19 workers, R = 50 plus the 200-trial long chunk; `dev/perf/harness/out/layer1_20261104_004741.csv`. It ran on the first `dev_p2` install (the worktree build before the DESCRIPTION line of 4.1 was added; the source files are identical, checked with `cmp`). Log `g4_layer1.log`.
- Additional, driver level (layer 1 tests the chunk, not the driver): `run_scenarios` on R01, R23, R04, R26 with
  `dev_l2` (phase 1) and `dev_p2` (default option, and option set explicitly to `"cpp"`), each at R = 2,000 / 4 workers
  (C++), R = 300 sequential (C++), and R = 120 / 4 workers with `use_cpp_engine = FALSE` (R body): tables and
  `.Random.seed` identical in 12 of 12 comparisons, both ways.

New tests (`test-simulate-chunk-exact.R`): xoshiro256++ outputs for state (1, 2, 3, 4) against an independent Python
implementation (first output `rotl(5, 23) + 1 = 41943041`); per-trial uniforms against the same implementation of the
seeding; uniforms in (0, 1), KS on 1e5, lag-1 correlation and KS of the first uniform over 4,000 consecutive trials;
replay identity with `.simulate_chunk_cpp` (R01, R23, both scenarios); split, block-order and single-trial invariance,
seed dependence, R's stream untouched, reducer errors on a missing or repeated index, output names equal to phase 1's;
refusal outside the supported set; driver: one `sample.int(…, 2)` per scenario, identical at 1 / 2 / 3 workers and
block sizes 9 / 64, seed dependence, bad option value; `use_cpp_engine = FALSE` gives the R body under either option.

### Sanity (smoke test, not G5; `p2_2_sanity.R`, `p2_2_sanity.csv`)

R = 2,000, cert seed 900001, sequential, both engines from `dev_p2`. z = (exact - cpp) / sqrt(var_exact + var_cpp).

| Design | Scenario | Arm | Metric | Phase 1 | Phase 2 | z |
|---|---|---|---|---|---|---|
| R01 | 1 (null) | Experimental | type I | 0.0115 | 0.0135 | 0.57 |
| R01 | 1 | Experimental | E[N] | 11.340 | 11.402 | 0.35 |
| R01 | 2 (alt) | Experimental | power | 0.810 | 0.823 | 1.06 |
| R01 | 2 | Experimental | E[N] | 27.955 | 28.248 | 1.10 |
| R23 | 1 (null) | Triplet | type I | 0.080 | 0.082 | 0.23 |
| R23 | 1 | Doublet / Triplet | E[N] | 29.20 / 29.89 | 29.10 / 29.91 | -0.25 / 0.04 |
| R23 | 2 (alt) | Triplet | power | 0.8075 | 0.8230 | 1.26 |
| R23 | 2 | Doublet / Triplet | E[N] | 43.93 / 46.99 | 44.42 / 47.42 | 1.52 / 1.37 |

All |z| <= 1.61 (PET futility included in the CSV). One pattern for G5 to watch, not a finding at this R: under the
alternative both designs show slightly higher power, fewer futility stops and larger E[N] with the exact engine (z 1.1
to 1.6, same sign in both designs). The R23 reference rows have power near 0.80, so the G5 classification check matters
there.

### Harness (`p2_cert.R`, R = 200, p2 mode)

`P2_LIB=<P6-field build> Rscript dev/perf/phase2/p2_cert.R R01 900001 200 p2 <out>` and the same for R23 (its first
reference seed is also 900001): both `engine_ok TRUE`, echo `exact[lines=200 pids=38 trials=400 expect=400 rd=FALSE F=6
lsf=TRUE mode=hc|vsref]`, `cpp[trials=0]`, all eight capture checks TRUE. R01 power 0.825, type I 0.015, E0N 11.50, E1N
28.28; R23 power 0.890, type I 0.050, E0N 56.18, E1N 95.19 (R = 200, noisy). Both checks were run first against a scratch library and then again against the final `dev_p2` (default `P2_LIB`), same echo. After the final install G3 (100 of 100) and the package tests (831 / 0 / 33) were rerun on it, same results. The harness was not edited.

### Timing (rough, not G6; `p2_2_timing.R`)

One process, `run_scenarios` on both scenarios, median of three after a warm-up. The machine was not quiet (load average
about 8 from other work; layer 1 was running), so these are rough; G6 should be run on a quiet machine with the step 7
scripts. Log `p2_2_timing.log`.

| Design | Phase 1 `dev_l2`, R = 1,000 | Phase 2 `dev_p2`, R = 1,000 | Phase 2, R = 10,000 | Ratio |
|---|---|---|---|---|
| R01 | 8.82 s | 0.079 s | 0.077 s | about 110x |
| R23 | 9.97 s | 0.109 s | 0.103 s | about 90x |

Seconds per 1,000 trials. The plan's G6 target (R23 under 8 s per 1,000 trials at one worker) is met by about 70x on
this measure; see P2-3 Section 5 item 4 for how the step 7 scripts define the quantity.

## 4. Mismatches between the plan, the brief and the code (reported, not improvised around)

1. **The harness needs a library field that a worktree build does not have.** `p2_cert.R` p2 mode requires
   `P6RecertEvolveTrialSHA = 7565a17` in the library's DESCRIPTION (P2-c), which only `layer2/build_dev_tree.sh` writes,
   and that script builds from a commit. P2-2 may not commit, and P2-1 had installed `dev_p2` without the field, so the
   first p2 run stopped ("is not a build_dev_tree.sh-style build on 7565a17"). I installed `dev_p2` from a copy of the
   worktree with that one line appended (lineage asserted as the script does) and no `EngineCppSHA`. The coordinator
   should rebuild `dev_p2` with `build_dev_tree.sh` from the P2-2 commit before the G5 run.
2. **mclapply.** The brief names `mclapply` in `R/simulation_driver.R`; the driver uses `parallel::parLapply` on a
   `makeCluster` (FORK on unix) or caller-supplied cluster. I kept that cluster code and used `clusterApplyLB` for
   dynamic dispatch, so external clusters (BATON) and PSOCK still work. `mclapply(mc.preschedule = FALSE)` would fork
   one process per block.
3. **"One draw" of R's generator.** The scenario seed is one `sample.int(.Machine$integer.max, 2L, replace = TRUE)` call
   (two integers, 62 bits), as the brief allows ("two ints"). It consumes two uniforms of R's stream (more only under
   rejection sampling), not one.
4. **Arm draw with one eligible arm.** Phase 1 draws a uniform even then; the exact engine does not (same
   probabilities). Single-arm designs therefore draw two numbers per patient (gap, event time).
5. **Cluster start-up dominates the exact engine.** At 0.08 to 0.11 s per 1,000 trials, a 19-worker FORK cluster built
   per scenario call (phase 1 behavior, kept) costs more than the trials: G3's 19-worker runs of 4,000 trials took 0.15
   to 1.7 s, mostly start-up. G6's 16-worker efficiency target will read poorly for that reason, not because of the
   engine. Not changed (the plan asks for the phase 1 driver structure); for the coordinator.
6. **P2-4 trigger.** The plan's rule (look slicing at 30% or more of R01 or R23 time) is met by R23 (75.5%), but R23
   already runs at about 0.11 to 0.14 s per 1,000 trials against the 8 s target, so P2-4's absolute gain is small.
   Whether to run it is the coordinator's call.
7. **Test count in the plan.** G4 in the plan says 723 passed; the brief says 783 (after P2-1). Both are history; the
   current count is 831.
8. **R version.** `dev_p2` is built with R 4.6.2 (P2-1 item 3); `dev_l2` with R 4.6.0. The driver-level identity
   above (dev_l2 vs dev_p2, default mode) held anyway.
9. **Warning text.** The status-1 warning is "Posterior beta parameter is non-positive. Setting to 1e-6.", once per
   call, without the interval number phase 1 prints. Unreachable in the supported set.
10. **Duplicate or NA indices.** Phase 1 accepts repeated or NA indices in a chunk; `.simulate_chunk_exact` stops on an
    NA index and `.exact_reduce` requires 1..n exactly once. The driver never produces either.

## 5. Profile (`p2_2_profile.R`, scratch timer build)

Scope timers (`std::chrono::steady_clock`) around `slice_arm`, the exact single-arm evaluation, the PH Laplace step, the
per-arrival rebalancing scan and the whole trial, in a scratch copy of the tree (the timers are not in the package).
`perf` is not usable (`perf_event_paranoid` = 4). One process, R = 2,000 per scenario, both scenarios.

| | R01 | R23 |
|---|---|---|
| wall per 1,000 trials (timer build) | 0.102 s | 0.135 s |
| look slicing (`slice_arm`) | **6.2%** (40,332 calls) | **75.5%** (329,770 calls) |
| exact single-arm probability | 90.6% (23,095 calls) | 0 |
| PH Laplace (`ph_beta_mode_var_cpp`, returns an R list) | 0 | 10.2% (79,935 calls) |
| rebalancing scan | 1.1% | 5.3% |

Timer overhead is about 40 to 50 ns per scope, under 5% of the totals. R23 slices both arms at every look whether or not
the gates pass (330k slices for 80k posterior evaluations), which is where the time goes.

## 6. Reproduce

```
cd ~/research/evolveTrial-perf
Rscript -e 'Rcpp::compileAttributes(".")'
# dev_p2 as built here (see 4.1; the coordinator should use build_dev_tree.sh from the commit instead):
#   copy git-tracked + new files (excluding dev/) to a tree, append the P6RecertEvolveTrialSHA line, R CMD INSTALL -l dev/perf/lib/dev_p2 <tree>
Rscript dev/perf/phase2/g3_check.R dev/perf/lib/dev_p2 2000                    # G3, about 40 s
Rscript dev/perf/harness/run_pkg_tests.R dev/perf/lib/dev_p2 <out.rds>          # G4 tests, about 60 s
Rscript dev/perf/harness/run_layer1.R --lib dev/perf/lib/dev_p2                 # G4 layer 1
Rscript dev/perf/phase2/p2_2_sanity.R dev/perf/lib/dev_p2 2000                  # sanity
Rscript dev/perf/phase2/p2_cert.R R01 900001 200 p2 <out.rds>                   # harness, also R23
Rscript dev/perf/phase2/p2_2_timing.R <lib> <cpp|exact> 1000 R01,R23            # timing
Rscript dev/perf/phase2/p2_2_profile.R <timer build lib> 2000                   # profile
```

The timer build: copy the tree, add `PT _pt(k)` scope timers at the five places in Section 5 and an `.exact_profile`
export, `compileAttributes`, install to a scratch library (done here under `/tmp/claude-1000/p2build/prof`).
