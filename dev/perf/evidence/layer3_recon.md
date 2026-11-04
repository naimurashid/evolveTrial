# Layer 3 reconnaissance: replaying E2c search evaluations (plan step 6)

Read-only research, 2026-11-03. Nothing was edited, installed or submitted. The Longleaf ssh control master was open,
so four read-only `ls`/`sha256sum`/`grep` checks were run there. No other Longleaf commands were run.

Abbreviations: AL = `~/Downloads/atbo-longleaf` (HEAD 1416e3d6), F0 = `~/Downloads/atbo-f0` (HEAD 1863cecb, the
E2C_F0_SHA), WT = `~/research/evolveTrial-perf`, PL = AL`/pilot`.

## 1. Headline answers

| Question | Answer |
|---|---|
| Engine used by E2c equals the port's reference? | **Yes.** E2C_ET_SHA `7565a17c625573d2a697d0d4a85efb2fa6c58ad9` + M4 patch md5 `83de2c8b2de1902630eddf00c334b508` (PL/e2c/E2C_FREEZE.md, Section 2). The Longleaf E2c library `/proj/rashidlab/nur2/e2_ll/R_libs_verify_m4/evolveTrial/DESCRIPTION` carries `P6RecertEvolveTrialSHA: 7565a17…` and `VerifyM4PatchMD5: 83de2c8b…` (Built R 4.5.0, 2026-10-21). These are the same two fields that `layer2_cert.R:120-122` requires of `lib/ref`, the plan's engine of record (plan Section 0). BATON 0.4.0 (same DESCRIPTION check), F0 1863cecb. Bit identity against the stored E2c values is the right expectation. |
| Where the ledgers are | On Longleaf at `/proj/rashidlab/nur2/e2_ll/e2c/results/PB-<method>-<seed>/sim_ledger.csv` (plus `results_certrecovery/`). Local copies at `~/pilot_r1/e2c/results_live/` (20 runs, with `f0_out/` and `run_state.rds`), checked against `~/pilot_r1/e2c/SHA256SUMS` (all 8 files for the candidates OK) and against Longleaf (sim_ledger sha256 prefixes equal for all three candidate runs: 0e0ddee4…, fb342584…, cc617992…). A third copy is in AL`/reproducibility/build/reproducibility_materials/data/e2c_raw/`. |
| Replayable in isolation? | **Yes.** Every search evaluation's seed is recorded in the ledger. The wrapper calls `set.seed(seed)` immediately before the simulator (F0 `R/warmstart_wrappers.R:885`), and `run_scenarios(..., seed = seed)` reseeds again (F0 `scripts/sim_fun_evolveTrial_unified.R:380-386`). No RNG state carries over from earlier evaluations. BATON's counter only chooses the seed, and the seed is stored. |
| Problems covered | E2c has **one problem only**: P-B, `doublet_triplet_admissible_fut`, which is **between-arm** (multi_arm wrapper, 9 coordinates). There is no single-arm problem in E2c. The candidates below therefore cover different runs, methods and K values. |

## 2. How an E2c search evaluation calls evolveTrial

### 2a. Job and environment

- Search jobs: `--cpus-per-task=16 --mem=48G --hint=nomultithread --constraint=amd`, general partition
  (PL/e2c/longleaf/submit_e2c.sh:176, `E2C_BIND` in PL/e2c/longleaf/e2c_env.sh). Certification jobs: 19 CPUs (:179).
- Job environment (PL/longleaf/job.sbatch:72-75): `env -i`-style block with `R_ENVIRON_USER=/dev/null
  R_PROFILE_USER=/dev/null`, `R_LIBS=$LL_FROZEN:$LL_COMP` (LL_FROZEN = `/proj/rashidlab/nur2/e2_ll/R_libs_verify_m4`),
  **`BATON_DISABLE_KSNAP=1`**, `SLURM_CPUS_PER_TASK=16`, `R_PARALLEL_PORT`, `LANG=C.UTF-8`, `TZ=America/New_York`.
  With `PILOT_E2C=1` it adds PILOT_E2C_FREEZE / _P6R / _M4S / _ENGINE_DIR (:77-83). Modules gcc/12.2.0 r/4.5.0 (freeze
  Section 1).
- BATON calls only: the wrapper is called inside `with_env(c(BATON_F0_CONF_ENV, BO_BUDGET_S1/2/3 = 700/300/200,
  SLURM_CPUS_PER_TASK = 16))` (PL/R/opt_baton.R:34-37, 236-238, 269). BATON_F0_CONF_ENV = BO_LOW/MED/HIGH_REPS
  3000/5000/10000, BO_USE_BUFFER TRUE, BO_CONSTRAINT_BUFFER 0.02, BO_VALIDATION_TIER_SIZES 5,5,
  BO_VALIDATION_MIN_STOP 10. **No BO_N_INIT and no BO_SEED** (layer 2 sets both, layer2_cert.R:108-112).
- Comparator calls (GA/SA/PSO/cBO): `problem$simulate(th, R, seed)` with no BO_* variables (PL/R/problems.R:427-431).
- What the wrapper path actually reads from the environment: `SLURM_CPUS_PER_TASK` / `MC_CORES` when the wrapper is
  built (F0 `R/warmstart_wrappers.R:57,67`, cached as `cached_num_workers`), `BO_SCALING_PROFILE` at parse time
  (F0 `R/scenario_parser.R:606`, unset in E2c, so "LITE"), and **`BATON_DISABLE_KSNAP` at call time** (F0
  `scripts/sim_fun_evolveTrial_unified.R:155-157`). The BO_* variables are read only by `run_three_stage_warmstart`
  (F0 `R/warmstart_core.R:1834-1977`), not by the wrapper. So the BATON and comparator paths run the same simulator
  call. The replay should still set the BATON block for BATON rows (cheap, and exact to the record), and must leave
  MC_CORES and BO_SCALING_PROFILE unset.
- Engine switches (PL/e2c/R/e2c_engine.R:127-222, called from `e2c_activate()` :259-301 before any problem is
  built): `verify_m4_install("off")` and then `p6r_install("p6", "probe")`, the same Tier 1 files and md5 values that
  layer 2 checks (p6r 52d2b43c…, m4s 405f7d60…). So `random_dropout = FALSE`, `min_follow_up_at_final = 6`, and option
  `evolveTrial.looks_stop_at_final = TRUE`, as in layer 2. E2c also installs an **outermost** echo wrapper on
  `run_scenarios` (it evaluates `base_args` once and passes the value, and it checks P6 frame events per call). The
  stored echo for the chosen jobs reads `p6[calls=1 F=6 opt=TRUE workers=32 … build=7565a17 lib=expected]
  dropout[calls=1 switch=off passed=FALSE …]` (32 = 16 FORK children x 2 scenarios), and SUMMARY.txt reads
  263/263 engine-ok (baton-1001) and 266/266 (baton-1003).

### 2b. Wrapper construction (differs from layer 2)

E2c builds the P-B wrapper this way (PL/R/problems.R:358-371, 372-432):
1. `source("R/scenario_parser.R")` in the F0 working directory (problems.R:68-74). That file itself sources
   `R/warmstart_wrappers.R` (which sources `scripts/sim_fun_evolveTrial_unified.R` and `R/scenario_configs.R`),
   `R/scenario_configs.R` and `R/globals_loader.R` (F0 scenario_parser.R:10-12).
2. With `SLURM_CPUS_PER_TASK=16` set: `parse_warmstart_scenarios("scenarios/cohortB_multiarm_production.csv",
   verbose = FALSE)`, with the default `filter_enabled = TRUE`, and pick `scenario_id == "doublet_triplet_admissible_fut"`.
3. Assertions: `cached_num_workers == 16`, `trial_type == "multi_arm"`, fidelity low/high = 3000/10000, weights
   0.5/1/0, objective `weighted_EN_N`.
4. The BATON process additionally sources `R/warmstart_core.R` and loads tidyverse and BATON
   (opt_baton.R:212-222). I checked that warmstart_core.R defines no function that the parser, wrapper or simulator
   files also define, so the wrapper path is the same in both processes.

Layer 2 replays Tier 1 R24, which is **this same scenario** (t1_tasks.tsv rows T13-T17, `doublet_triplet_admissible_fut`,
cohortB_multiarm_production.csv line 34). It used protocol `bws`: it sources four files in its own order and calls
`parse_single_scenario` on the filtered row (layer2_cert.R:140-150). This is the I2 issue again, so layer 3 must use
the E2c parse path above and not layer 2's.

### 2c. Seeds

- BATON: root `baton_rng_seed(s) = 1e9 + s * 1e5` (opt_baton.R:44-48). BATON 0.4.0 `bo_calibrate` passes
  `root + eval_counter` (opt_baton.R:14-15). Checked on all five BATON ledgers: every search row has
  `seed == 1e9 + search_seed*1e5 + eval_index` (243, 198, 246, 249 and 186 rows, zero exceptions). Re-screen rows
  (BATON validation, R 10,000) use root + 15001… instead.
- Comparators: `sim_seed("search", s, i) = s * 1e5 + i` (PL/R/evaluator.R:67-71), so e.g. GA seed 1001, eval 274 gives 100100274.
- A BATON evaluation is recorded by the decorator after the call (opt_baton.R:69-108, `evaluate(passthrough=…,
  rng_isolate = FALSE)`). The seed column is the seed that BATON passed to the wrapper.

### 2d. R, workers, chunking

R = 3,000 (`fidelity = "low"`) for every E2c search row (the five BATON runs all stopped early in stage 1 on
"validation_success", so no stage-3 evaluations at high R exist). Workers = 16, giving 16 chunks of 185x2, 187x2 and
188x12 per scenario (rng_driver_inventory.md Section 2). Each evaluation makes one `run_scenarios` call with two
scenarios.

## 3. What the ledger stores, and at what precision

`sim_ledger.csv` is written by `utils::write.table(sep=",", qmethod="double")` (PL/R/evaluator.R:150-157; row built
at :266-283). Columns: row_id, run_id, method, problem, search_seed, eval_index, purpose, R, seed, theta_total_n,
theta_eff, theta_fut, theta_ev, theta_beat, theta_margin, theta_hr_margin, theta_prior_strength, theta_n_intervals,
objective, power, type1, power_se, type1_se, objective_se, n_rep, error_flag, invalid_design, error_msg, n_warnings,
wall_s, cpus, cpus_from_job, worker_count, alloc_core_s, port, host, git_sha_H, git_sha_F0, timestamp, cum_spend,
rng_isolate.

- write.table prints doubles to **15 significant digits**. So power and type1 recover exact integer counts (x 3000;
  integral on every search row checked), but objective and the *_se columns are rounded at 15 digits.
- **The theta columns are also rounded at 15 digits.** Compared with the theta that was actually simulated, the
  maximum relative difference is 4.3e-15 to 4.8e-15 per BATON run (for example eff 0.98155195893057501 in the ledger
  against 0.98155195893057512 simulated). The exact theta is in:
  - BATON: `f0_out/doublet_triplet_admissible_fut_warmstart_result.rds`, `$history_combined`. The `fidelity == "low"`
    rows are in ledger order with `eval_id == eval_index` (checked for all 5 runs). The list column `theta[[i]]` is
    identical to the per-parameter columns. The same table also holds `power`, `type1`, `objective`, **`EN_null`, `EN_alt`**
    at full double precision.
  - Comparators: `run_state.rds`, `$history` (columns eval_index, the 9 parameters, objective, power, type1,
    search_feasible, ledger_row) at full precision. It has no EN_null or EN_alt.
- The ledger has **no EN_null or EN_alt**. For comparator rows the comparison is limited to counts and objective
  (objective is a deterministic function of total_n and EN_null under weights 0.5/1/0).

## 4. Conflicts with the plan and with layer 2 (reported, not worked around)

1. **Plan step 6 says "from the stored ledgers … compared with the ledger values".** The ledger theta is rounded at
   15 significant digits, so replaying from the ledger would simulate a theta that differs in the last bits from the
   one E2c simulated. A continuous threshold off by one ulp can flip a decision. The replay must take theta from the
   rds files in Section 3, and the comparison should use the rds values (full double) plus the ledger counts.
2. **The ledger has no EN_null or EN_alt**, which layer 2 compares. Proposal: compare power_count and type1_count
   (integers) and objective, EN_null and EN_alt with `identical()` on the full doubles from history_combined. For
   comparator rows, counts and objective only. This is stricter than layer 2's 4-decimal rule, and it is achievable
   because the stored values are full precision.
3. **Parse path.** Layer 2 for this scenario used `bws` (parse_single_scenario). E2c used parse_warmstart_scenarios
   (Section 2b). Layer 3 must copy E2c's path.
4. **Environment block.** Layer 2 sets BO_N_INIT=90 and BO_SEED=<task seed> (layer2_cert.R:108-112). E2c set neither
   (BATON rows: BATON_F0_CONF_ENV plus budgets; comparator rows: no BO_*). Neither variable reaches the wrapper path
   (Section 2a), but layer 3 should use E2c's block exactly.
5. **CPU binding.** E2c searches ran with `--constraint=amd` (E2C_BIND). layer2.sbatch:1-9 sets `--hint=nomultithread`
   but not `--constraint=amd`. Layer 3 should add it to match the stored rows' hardware class (hosts c1514xx).
6. **E2c's own engine check cannot be reused with the dev library.** The echo in e2c_engine.R requires P6 frame
   events with F 6 from the R chunk body, and the C++ body produces none (layer2_cert.R header L2-e). Use layer 2's
   L2-d/L2-e probe and check (zero interim_check probe lines, zero P6 frame events, C++ trials = 2 x 3000, cfg
   FALSE / 6 / TRUE).

None of these touches the engine identity question. The engine matches.

## 5. Recipe: `layer3_search.R`, derived from layer2_cert.R with the changes tagged

Usage: `Rscript layer3_search.R <run_id> <eval_index> <dev|ref> <out.rds>`, one evaluation per process.

Keep from layer2_cert.R unchanged: L2-b paths via env vars (L3_AL, L3_F0, L3_LIB, L3_BATON_LIB), L2-c library
identity (ref: 7565a17 + 83de2c8b; dev: P6RecertEvolveTrialSHA 7565a17 + EngineCppSHA), L2-d C++ probe wrapper, L2-e
engine check, L2-f FORCE_R negative control, the RNGkind assertion, the "no packages loaded before the env block"
assertion, and the F0 git SHA check.

Change:
- **L3-a input row.** Read `sim_ledger.csv` and the rds for `<run_id>` from L3_RESULTS (default
  `~/pilot_r1/e2c/results_live`; on Longleaf `/proj/rashidlab/nur2/e2_ll/e2c/results`). Select `purpose == "search"
  & eval_index == i`, and assert R == 3000, worker_count == 16, error_flag FALSE, invalid_design FALSE, git_sha_F0 ==
  1863cecb. Seed = ledger `seed`, then assert the derivation (BATON: `1e9 + search_seed*1e5 + i`; others:
  `search_seed*1e5 + i`). Theta = `history_combined$theta[[k]]` (BATON, `fidelity=="low" & eval_id==i`) or the
  run_state history row (comparators). Assert that theta agrees with the ledger theta to 1e-13 relative, and record the
  sha256 of both input files.
- **L3-b env block** (set before any package load). `BATON_DISABLE_KSNAP=1`, `SLURM_CPUS_PER_TASK=16`, MC_CORES and
  BO_SCALING_PROFILE unset, no stray BO_*. For BATON rows also BATON_F0_CONF_ENV plus BO_BUDGET_S1/2/3=700/300/200,
  and no BO_N_INIT and no BO_SEED. WORKERS = 16, asserted against SLURM_CPUS_PER_TASK.
- **L3-c wrapper build** as in Section 2b (source scenario_parser.R in F0, parse_warmstart_scenarios on
  cohortB_multiarm_production.csv, pick doublet_triplet_admissible_fut). For BATON rows also source warmstart_core.R
  and load BATON 0.4.0. Assert `cached_num_workers == 16`, multi_arm, weights 0.5/1/0, fidelity 3000/10000.
- **L3-d switches.** Install the capture hook and then `verify_m4_install("off")` and `p6r_install("p6","probe")` as
  layer 2 does (innermost hook). Optionally build the wrapper after the switches, as E2c did (no semantic difference).
- **L3-e call.** In the F0 working directory: `m <- W(theta, fidelity = "low", seed = SEED)`.
- **L3-f checks.** One run_scenarios call with 2 scenarios. Cutpoints = K (unsnapped, so K = 15 really runs 15
  intervals). random_dropout FALSE, min_follow_up_at_final 6, looks_stop_at_final TRUE, `num_workers == 16`,
  `num_simulations == 3000`, counts integral.
- **L3-g compare** (separate `layer3_compare.R`). `power_count`, `type1_count` equal to the ledger counts, and
  `identical()` on objective, EN_null, EN_alt (BATON) or objective (comparator) against the rds values. Pass means
  every row passes.
- sbatch: copy layer2.sbatch with `--cpus-per-task=16`, `--constraint=amd`, `--hint=nomultithread`, r/4.5.0,
  `--time=00:15:00`, and an array of 3 tasks.

Local dry run first. I6 showed that R 4.6.0 reproduced R 4.5.0 rows, so a local `ref` run of candidate 2 at 16
workers (about 40 s) is a cheap check of the recipe before Longleaf.

## 6. Candidate evaluations (all P-B, between-arm, R 3,000, 16 workers)

Values are from the local copies (sha256-verified). Counts come from the ledger. The doubles come from the rds at full
precision (%.17g).

**C1. PB-baton-1001, eval 1** (BATON, first evaluation of the process, K = 8, a snap value, so it serves as the
baseline). Seed 1100100001, host c151407, wall 95.8 s, ledger row 1.
theta: total_n=169; eff=0.98155195893057512; fut=0.69588270453425749; ev=10; beat=5.007874832100546;
margin=0.13693782703355989; hr_margin=0.13919840686398352; prior_strength=0.29445178564637897; n_intervals=8.
Expected: power_count 2369, type1_count 72, objective 107.83444444444444, EN_null 80.913333333333327,
EN_alt 135.74266666666665.

**C2. PB-baton-1003, eval 95** (BATON, K = 15, outside {4,6,8,12}, so it tests the unsnapped path; search-feasible).
Seed 1100300095, host c151412, wall 40.7 s.
theta: total_n=96; eff=0.96252958471348915; fut=0.37682610586226611; ev=12; beat=11.410913408958198;
margin=0.14067301857802603; hr_margin=0.011587663647343611; prior_strength=0.29970528989394662; n_intervals=15.
Expected: power_count 2476, type1_count 260, objective 57.69155555555556, EN_null 71.472000000000008,
EN_alt 90.003999999999991.

**C3. PB-ga-1001, eval 274** (GA comparator path `problem$simulate`, no BO_*, K = 15). Seed 100100274, host c151415,
wall 40.4 s. Theta from `run_state.rds$history`:
total_n=121; eff=0.97500158712464924; fut=0.4703726912881312; ev=8; beat=14.319480040939588;
margin=0.44552224744691205; hr_margin=0.15855964914027593; prior_strength=0.45450102810635407; n_intervals=15.
Expected: power_count 2571, type1_count 211, objective 75.913277777777779 (no EN values stored).

Alternates: baton-1001 eval 58 (K 13, seed 1100100058, 50.2 s), baton-1003 eval 122 (K 15, 46.3 s), ga-1001 eval
270 (K 13, 46.8 s). Every E2c search row checked has error_flag FALSE, invalid_design FALSE and n_warnings 0.

## 7. Cost

Stored reference-engine walls for E2c search evaluations at 16 workers: median 34-54 s per BATON run, 90th
percentile about 95 s, maximum 257 s. GA-1001: median 63 s, maximum 398 s. The candidates take 95.8, 40.7 and
40.4 s under the reference engine. Layer 2's C++ engine ran R24 (this scenario) at R 10,000 and 19 workers in 4.5-16 s,
against about 56 s for the reference engine. Layer 2 billed 7.12 core-h for 75 tasks, about 0.095 core-h per task
including start-up. Estimate for 3 dev tasks at 16 CPUs: about 0.3 core-h, and at most about 0.5 core-h if start-up
dominates. That is inside the plan's "under 1 core-hour". Adding a Longleaf `ref` control for the 3 rows would cost
about 0.8 core-h more (177 s x 16 plus overhead) and is not needed, because the stored rows are the reference.

## 8. Open items for the coordinator

- Decide whether the comparison is `identical()` on the full doubles (proposed, since the stored values are full
  precision) or layer 2's 4-decimal rule.
- Confirm that the plan wording "ledger values" may be read as ledger counts plus the run's rds values (Conflict 1).
- C1 is the first call in its process (its 95.8 s includes the first FORK cluster). This has no bearing on results,
  because every call starts from `set.seed`.
