# evolveTrial simulator profile (engine 58a4850, P6 + dropout off), 2026-11-03

Scratch dir: `/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf/` (called `$S` below).
No repository was edited or committed. evolveTrial was installed only into `$S/Rlib`.

## Headline

| Design | Path | s per 1,000 trials, 1 worker (R = 2,400 per scenario, production FORK path) | null / alt split | 16-worker efficiency (both scenarios) |
|---|---|---|---|---|
| R01 cohortA_H0_fut (Scenario A, single_arm) | `interim_check_hc`, C++ posterior | **13.4** | 4.9 / 21.9 | 0.64 (speedup 10.3) |
| R23 doublet_triplet_H0_fut (Scenario B, multi_arm vs-ref PH) | `interim_check_vs_ref`, **R-level** posterior | **75.1** | 47.1 / 103.1 | 0.56 (speedup 8.9) |

At the certification fidelity (R = 10,000 per scenario, 20,000 trials per wrapper call), one evaluation costs about 268 core-s (R01) and 1,500 core-s (R23) at 1-worker speed.
In-process sequential runs (`parallel_replicates = FALSE`, R = 300, both scenarios) gave 13.4 (R01) and 85.5 (R23) s per 1,000. The R23 gap from 75.1 is the trial mix at R = 300. A 1-worker FORK run at R = 300 gave 77 s per 1,000, so most of the gap comes from the sample, not the execution mode.

Main finding. In the multi-arm vs-ref path, the posterior sampler that runs is the pure-R `sample_vs_ref_medians_ph` from `R/posterior_helpers.R:105`. DESCRIPTION has no Collate field, so files collate alphabetically. That means `posterior_helpers.R` loads after `posterior_cpp_dispatchers.R` and overwrites the C++ dispatcher of the same name (`posterior_cpp_dispatchers.R:241`, which calls `sample_vs_ref_medians_ph_cpp`). I confirmed this in the installed build: `getSrcFilename(evolveTrial:::sample_vs_ref_medians_ph)` returns `posterior_helpers.R` line 105, and `.use_cpp_posterior()` returns TRUE.
That one function is 67% of R23 wall time. The 800-draw R loop calls `rgamma` and the R `calculate_median_survival_piecewise` twice per draw.
Exploratory check, in memory only in a scratch session: I swapped in the dispatcher body with `assignInNamespace`. R23 then ran at **34.7 s per 1,000 (2.5x faster)**, but the operating characteristics at the same seed were **not identical** (R = 300: type I 0.0867 -> 0.090, power 0.790 -> 0.750, PET_fut also moved). So the C++ sampler uses the RNG stream differently or computes the mode differently. Adopting it changes the engine's output and would need an equivalence (distributional) check before use.

## Setup (exact commands)

```bash
S=/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf
git -C ~/research/evolveTrial worktree add $S/et58 58a4850
MAKEFLAGS=-j8 R_LIBS=$S/Rlib R CMD INSTALL --with-keep.source --no-test-load -l $S/Rlib $S/et58   # -O2 (system default flags)
R_LIBS=$S/Rlib Rscript -e 'library(evolveTrial); find.package("evolveTrial")'                   # -> $S/Rlib/evolveTrial, 0.0.0.9001
R_LIBS=$S/Rlib Rscript $S/capture.R R01 $S/args_R01.rds      # build F0 wrapper, capture run_scenarios call
R_LIBS=$S/Rlib Rscript $S/capture.R R23 $S/args_R23.rds
R_LIBS=$S/Rlib Rscript $S/profile.R R01 300        # Rprof(interval 0.01, line.profiling = TRUE), sequential
R_LIBS=$S/Rlib Rscript $S/profile.R R23 300
R_LIBS=$S/Rlib Rscript $S/profile.R R01 300 noprof # timing without Rprof (also R23)
Rscript $S/bucket.R $S/rprof_R01_R300.out $S/buckets_R01_R300.csv   # stack-based buckets (also R23)
$S/run_scale.sh > $S/scale.log                      # 32 runs: {R01,R23} x scen {1,2} x w {1,4,8,16} x {internal,external}, R = 2,400
R_LIBS=$S/Rlib Rscript $S/counts.R R01 100          # per-trial call counts (also R23)
R_LIBS=$S/Rlib Rscript $S/counts.R R23 300 swap     # exploratory C++ dispatcher swap
```

Machine: Intel Core Ultra 7 265, 20 cores (8 performance + 12 efficiency cores, no SMT), R 4.6.0.

### How the call was built (capture.R)

- Designs come from `atbo-longleaf/results/retune_t1/cert_design.csv`, row `cert_seed == 900001`, `theta_json`. They are R01 (task T04, protocol xphil) and R23 (task T11, protocol bws). The task rows come from `pilot/retune/staged/t1_tasks.tsv`.
- F0 code at `/home/naimrashid/Downloads/atbo-f0` (1863cecb) is sourced in that protocol's order (the `SRC` list copied from `pilot/retune/R/retune_run.R`). The scenario row is parsed with `parse_single_scenario()` from `scenarios/cohortB_singlearm_production.csv` or `scenarios/cohortB_multiarm_production.csv`, and the wrapper is `spec$wrapper`. `SLURM_CPUS_PER_TASK = 1` gives `cached_num_workers = 1`.
- Engine settings match what the retune certify phase installs: `options(evolveTrial.looks_stop_at_final = TRUE)`, plus a `run_scenarios` wrapper that asserts the incoming `min_follow_up_at_final == 0` and `random_dropout` absent, then sets `min_follow_up_at_final = 6` and `random_dropout = FALSE`.
- The wrapper was called as `W0(as.list(theta), fidelity = "high", seed = 900001)`, with `fidelity_levels["high"]` set to 300 in a child environment (the TEST_R pattern). The captured `base_args`, `scens` and `run_scenarios(..., parallel = FALSE, seed = 900001, return_variance = TRUE)` were then replayed directly with `num_simulations`, `parallel_replicates`, `num_workers` (and `cluster`) changed.
- Full argument lists: `$S/args_R01.txt`, `$S/args_R23.txt` (dput). The key values follow.

R01 theta: `total_n 32, eff 0.893541, fut 0.345921, ev 3, beat 1.285360, margin 0.727348, prior_strength 0.216732, n_intervals 8`. Single-arm `Experimental`, 5,000 posterior draws, beat 1.285, max_trial_time 72, null and alt scenarios.
R23 theta: `total_n 103, eff 0.950199, fut 0.372153, ev 11, beat 1.005807, margin 0.272355, hr_margin 0.033790, prior_strength 0.305119, n_intervals 15`. Arms Doublet (reference) and Triplet, `compare_arms_option = TRUE`, `use_ph_model_vs_ref = TRUE`, 800 posterior draws, max 52/51 per arm, accrual 2/month, 12 intervals on the 24-month grid, `rebalance_after_events = 50`. Scenarios are medians (6, 6) for null and (6, 9) for alt.
Both: `parallel_replicates = TRUE, num_workers = 1L, cluster_type = "auto"` as built by the F0 wrapper (no `cluster` object), `min_follow_up_at_final = 6, random_dropout = FALSE`.

Sanity check at R = 300 (wrapper): R01 power 0.830, type I 0.033. R23 power 0.800, type I 0.083. These agree with the certified values (0.821 / 0.014 and 0.804 / 0.082) within R = 300 Monte Carlo error.

## 1. Single-worker profile, R = 300 per scenario (600 trials)

| Design | wall (no Rprof) | wall (Rprof) | Rprof samples |
|---|---|---|---|
| R01 | 8.03 s | 8.20 s | 807 |
| R23 | 51.3 s | 53.2 s | 5,243 |

Raw outputs: `$S/rprof_<D>_R300.out`, `$S/summaryRprof_<D>_R300.txt` (by.self and by.total top 40, by-line self top 40/60), `$S/summaryRprof_<D>_R300.rds` (function and line summaries), `$S/buckets_<D>_R300.csv`.

### Top self time (function level)

R01: `draw_posterior_hazard_samples_cpp` 61.1%, `[<-.data.frame` 3.5%, `calculate_median_survival_matrix_cpp` 2.5%, `simulate_chunk` 2.4%, `pmin` 2.2%, `data.frame` 1.4%, `vapply` 1.4%, `extract_events_fast` 1.2%, `interim_check_hc` 1.0%, `[<-` 1.0%, `[.data.frame` 1.0%. Everything else is under 1% each.
R23: `calculate_median_survival_piecewise` (R) 31.0%, `rgamma` 17.6%, `sample_vs_ref_medians_ph` (R, self) 14.8%, `pmin` 2.2%, `length` 2.2%, `[<-.data.frame` 1.9%, `vapply` 1.8%, `simulate_chunk` 1.7%, `calculate_interval_metrics_fast` 1.5%, `[.data.frame` 1.4%, `%in%` 1.2%, `/` 1.0%, `data.frame` 1.0%.

### Top total time

R01: `interim_check_hc` 76.0%, `calculate_current_probs_hc` 60.4%, `draw_posterior_hazard_samples` 61.3%, `slice_arm_data_at_time` 14.5%, `data.frame` 8.1%, `extract_events_fast` 6.4%, `[<-` 6.2%, `get_active_registry` 4.0%.
R23: `interim_check_vs_ref` 83.8%, `sample_vs_ref_medians_ph` 66.7%, `calculate_median_survival_piecewise` 32.2%, `rgamma` 18.9%, `slice_arm_data_at_time` 17.4%, `data.frame` 6.7%, `extract_events_fast` 5.8%, `get_active_registry` 5.3%, `calculate_interval_metrics_fast` 4.6%.

### Hot lines (line profiling, package installed with keep.source)

- R23 `posterior_helpers.R#130`: `med_ctrl[i] <- calculate_median_survival_piecewise(...)` and the `rgamma` line around it (#129-132) take 23% self and about 40% total. `posterior_helpers.R#19-44` is the body of the per-interval R loop in `calculate_median_survival_piecewise`, about 18%.
- `state_management.R#137` (the `data.frame(observed_time, event_status)` construction in `slice_arm_data_at_time`): 7.1% (R01) and 6.4% (R23). `state_management.R#111` (`get_active_registry` row subset of the pre-allocated registry data.frame): 3.6% and 4.8%. `#144` (`stats::median` of follow-up): about 1.4-2%. `#54` (person-time loop, `pmin`): 1.1-2.2%.
- `simulation_driver.R#672-675` (four `state$registries[[arm]][idx, col] <- value` data.frame cell assignments per enrolled patient): 7.0% (R01) and 4.2% (R23).
- `simulation_driver.R#503-514` (`extract_events_fast`, called on every enrollment until rebalancing is done): 2.5-3%.

### Self time by bucket (each Rprof sample assigned from its full call stack)

| Bucket | R01 | R23 |
|---|---|---|
| Posterior sampling, C++ call (`draw_posterior_hazard_samples_cpp`, median matrix) | **63.6%** | 0% |
| Posterior sampling, R level: median-survival loop | 0% | **32.2%** |
| Posterior sampling, R level: RNG (`rgamma`, `rnorm` in the sampler) | 0% | **18.0%** |
| Posterior sampling, R level: sampler loop other / wrapper | 1.2% | **17.1%** |
| Look logic: data slice (`slice_arm_data_at_time`, `get_active_registry`, `calculate_interval_metrics_fast`) | **17.2%** | **20.4%** |
| Look logic: `extract_events_fast` (rebalance event count) | 3.7% | 2.8% |
| Look logic: gates and decision code | 2.5% | 1.9% |
| Per-trial loop: registry/state data.frame bookkeeping | 7.3% | 4.2% |
| Per-trial loop: other R code (`is_eligible`, `sample` set-up, etc.) | 4.1% | 2.9% |
| Per-trial loop: patient RNG (`rexp`, `rweibull`, `runif`, `sample`) | 0.4% | 0.5% |
| Outside `simulate_chunk` (setup, aggregation) | ~0% | ~0% |

Generic R overhead (innermost self time, cutting across buckets), R01 / R23: `[<-.data.frame` 3.5 / 1.9%, `[.data.frame` 1.0 / 1.4%, `[[.data.frame` 0.6 / 0.6%, `data.frame` + `as.data.frame` + `deparse`/`.deparseOpts`/`make.names` (the deparse cost comes from `data.frame()` naming its columns) 3.4 / 2.5%, `pmin`/`pmax` 2.4 / 2.5%, `vapply` 1.4 / 1.8%, `%in%` 0.7 / 1.2%, `match.arg` 0.2 / 0.3%, `sys.call` 0.7 / 0.7%. `rbind` 0 / 0 and `do.call` 0 (self) in both. All data.frame machinery together is about 9% (R01) and 7% (R23).

Per-trial call counts (`counts.R`, R = 100 per scenario):
- R01: 9.1 `interim_check`, 10.1 `slice_arm_data_at_time`, 5.9 posterior sampler calls of 5,000 draws each.
- R23: 41.5 `interim_check`, 81.8 `slice_arm_data_at_time` (two arms per look), 19.8 vs-ref posterior calls of 800 draws each. That is about 16,000 R-level loop iterations per trial, each with an `rgamma` call and two R median computations.

## 2. Parallel scaling, R = 2,400 per scenario

Each scenario is a separate process, so the R23 1-worker runs stay under 5 minutes (113 s and 248 s).
- **internal**: the production path. `parallel_replicates = TRUE`, `num_workers = w`, no cluster object, so `run_simulation_pure` calls `parallel::makeCluster(w, "FORK")` and stops it inside **every** call. That is two clusters per F0 wrapper evaluation (null and alt). CPU = master `proc.time` self + child (FORK children are reaped, so `user.child + sys.child` is populated).
- **external**: `create_simulation_cluster(w, "FORK")` is called once and passed as `cluster`. Worker CPU is read by `clusterCall(proc.time)` before and after.

Cluster set-up: 0.009 s (1 worker) to 0.05 s (16 workers), and stop 0.004-0.028 s. Creating a fresh FORK cluster per call is therefore not a cost. The code comment at `simulation_driver.R:961` ("cluster spawn/teardown ~5-10 sec") describes PSOCK, not FORK. Master CPU is under 0.1 s in every run.

| Design, scen | mode | w=1 wall | w=4 | w=8 | w=16 | eff 4 / 8 / 16 | CPU per rep, ms (w = 1 -> 16) |
|---|---|---|---|---|---|---|---|
| R01 null | internal | 11.69 s | 3.16 | 1.77 | 1.30 | 0.93 / 0.83 / 0.56 | 4.87 -> 7.02 (1.44x) |
| R01 null | external | 11.80 | 3.21 | 1.76 | 1.28 | 0.92 / 0.84 / 0.58 | 4.90 -> 6.87 (1.40x) |
| R01 alt | internal | 52.64 | 13.12 | 7.65 | 4.97 | 1.00 / 0.86 / 0.66 | 21.9 -> 29.4 (1.34x) |
| R01 alt | external | 52.45 | 13.78 | 7.46 | 4.92 | 0.95 / 0.88 / 0.67 | 21.7 -> 29.2 (1.34x) |
| R23 null | internal | 112.96 | 28.56 | 17.34 | 12.68 | 0.99 / 0.81 / 0.56 | 46.6 -> 71.1 (1.52x) |
| R23 null | external | 120.60 | 30.14 | 18.16 | 12.88 | 1.00 / 0.83 / 0.59 | 49.5 -> 73.8 (1.49x) |
| R23 alt | internal | 247.52 | 68.02 | 37.86 | 27.68 | 0.91 / 0.82 / 0.56 | 102.8 -> 165.2 (1.61x) |
| R23 alt | external | 258.65 | 72.07 | 40.86 | 30.16 | 0.90 / 0.79 / 0.54 | 107.4 -> 182.0 (1.69x) |

Both scenarios together (internal): R01 64.3 s -> 6.27 s at 16 workers (speedup 10.3, efficiency 0.64). R23 360.5 s -> 40.4 s (speedup 8.9, efficiency 0.56).
Full table: `$S/scaling_table.csv`; per-run rows `$S/scale/*.rds`; log `$S/scale.log`.

Reading the scaling. Wall time follows per-replicate CPU almost exactly. CPU per replicate rises 34-69% from 1 to 16 workers, while master CPU and cluster overhead are negligible. So the loss of efficiency is per-core throughput, not parallel bookkeeping. On this CPU, workers beyond 8 land on efficiency cores, and all-core turbo drops. Efficiency at 4 workers (0.90-1.00) shows the code itself scales. These numbers will differ on Longleaf's homogeneous nodes. Chunking is one equal chunk per worker (`splitIndices`), so trial-length variance also leaves some workers idle at the end. That probably matters most for the null scenarios, where trials stop at varied times.

## 3. Seconds per 1,000 simulated trials at 1 worker

- R01: **13.4 s** (null 4.9, alt 21.9). Sequential in-process R = 300 gives 13.4.
- R23: **75.1 s** (null 47.1, alt 103.1). Sequential in-process R = 300 gives 85.5. A 1-worker FORK run at R = 300 gives 77.1, so the difference is mostly the R = 300 trial sample.
- R23 with the shadowed C++ dispatcher swapped in (exploratory, results not identical): 34.7 s (sequential, R = 300).

## 4. Growth patterns

- No `rbind` inside the trial loop (0 samples on any stack). The registry is pre-allocated (`max_total_patients_per_arm` rows per arm) and filled by index.
- Per-patient cell assignment `state$registries[[arm]][idx, col] <- v` goes through `[<-.data.frame` four times per patient, and each call copies the column. The cost is linear in N per patient, so roughly O(N^2) per trial, but N is 32-103 here (4-7% of time).
- Every look calls `get_active_registry()` (a row subset that makes a new data.frame), builds another `data.frame` in `slice_arm_data_at_time`, and recomputes interval metrics over **all** patients enrolled so far. The number of looks grows with trial duration (R23: 41.5 looks and 82 slices per trial), so slicing is O(N x looks) per trial, which is the 17-20% "data slice" bucket. Incremental per-arm summaries or plain vectors in place of data.frames would remove most of it.
- `extract_events_fast` runs over the whole registry on every enrollment until 50 events have accrued (O(N^2) per trial, about 3%).
- The posterior sampler cost per look does not depend on N. It is (draws x intervals), and in R23 the draw loop is in R.

## Where this brief conflicts with the code

1. **Engine of record.** `results/retune_t1/t1_provenance.csv` records evolveTrial `7565a17` + `m4_nodropout.patch` (md5 83de2c8b), and `retune_run.R` refuses any library whose DESCRIPTION lacks `VerifyM4PatchMD5` / `P6RecertEvolveTrialSHA`. The brief's 58a4850 is 7565a17 plus the `random_dropout` formal (the only diff is in `R/simulation_driver.R`), and a plain install of it carries neither build field. So I did not use `p6r_install()` or `verify_m4_install()`, which would fail closed or add per-call probes. I reproduced their effect with a small `run_scenarios` wrapper instead. Timings exclude the harness's interim_check instrumentation.
2. **"Between-arm".** Scenario B `doublet_triplet_*` runs with `trial_type = multi_arm` (`compare_arms_option = TRUE`, `interim_check_vs_ref`, PH vs-reference posterior) through `run_simulation_pure`. It is not the hybrid SA-to-BA path (`compute_hybrid_oc_cpp` / `hybrid_sim_rcpp.R`). If "between-arm" meant the hybrid engine, that path was not measured.
3. **Where design rows live.** The theta values are in `results/retune_t1/cert_design.csv` and `t1_summary_designs.csv` (`theta_json`). `t1_provenance.csv` holds build provenance only. The P6 `p6recert.R` driver targets the 32 E2 P-B designs (designs.rds), not the T1 designs, so the call construction was taken from `pilot/retune/R/retune_run.R` (certify phase).
4. **"R = 300 simulations".** Each wrapper or `run_scenarios` call simulates two scenarios, so R = 300 is 600 trials per profile run. The scaling runs used R = 2,400 per scenario, one scenario per process.
5. **`create_simulation_cluster`** exists, but the production F0 wrapper never passes a cluster. A fresh FORK cluster is created per `run_simulation_pure` call (two per evaluation), at negligible cost.

Cleanup: `git -C ~/research/evolveTrial worktree remove $S/et58` was run at the end. The `$S/../et58_ro` worktree was there before this task and was left alone.
