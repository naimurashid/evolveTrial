# Review of the layer 3 harness (dev/perf/layer3), 2026-11-03

Scope: can the layer 3 gate report PASS on vacuous, wrong or stale input, and does `layer3_search.R` change the E2c
call path beyond its L3-a to L3-g hunks. Read-only review: grep, line ranges, `diff layer2_cert.R layer3_search.R`,
`git status --porcelain`, `git diff --quiet edbbb83 HEAD -- dev/perf/layer3` (exit 0, harness unchanged since the
staged commit), and one local `readRDS` of `out/local/dev/*.rds`. Nothing was run on Longleaf, so the Longleaf-side
state (staged md5 files, `$L2/lib` binary) is UNCONFIRMED. Call-path claims about E2c were checked against
`~/Downloads/atbo-longleaf/pilot/R/problems.R:358-432`, `pilot/R/opt_baton.R:1-110,205-275`, the F0 worktree at
1863cecb, and BATON 0.4.0 `invoke_simulator`.

Counts: Critical 0, High 0, Medium 0, Low 8, plus 4 notes. No path was found by which the gate reports PASS on
vacuous, wrong-engine, mismatched, missing or reduced-precision input through the scripts as written.

## Primary question: paths to a false PASS

| path | verdict | where it is stopped |
|---|---|---|
| R engine silently used instead of C++ | not reachable | The dispatcher (evolveTrial 1c91871 `R/simulation_driver.R:955-957`) runs C++ only when the option is TRUE and the config is supported, else the R body. Any R-body trial leaves `cp_trials < nscen * Rf` (search:200, 212) and the run exits 3 before saving. `simulate_chunk_cpp` never falls back (driver:947), so no chunk is counted by the C++ probe and then recomputed in R. `L3_FORCE_R` is dropped by `env -i` (sbatch:139-143), and compare rejects `force_r` rows (compare:48). Local negative control confirms the probe count reaches 0 (record Section 5). |
| wrong engine commit | not reachable via the scripts | sbatch:125 and stage:174 grep `EngineCppSHA`; search:101-103 requires `EngineCppSHA == L3_ENGINE_SHA`, the P6 field, no M4 field and the C++ entry point; compare:48 and :57 require every row's `engine_sha` to equal `L3_ENGINE_SHA` and be unique. Layer 2 finding M2 is fixed. Residual: the identity is a DESCRIPTION label, see L3. |
| ledger rows mismatched to candidates | not reachable | Compare joins on the output's own `run_id`/`eval_index` (compare:20-21), requires exactly one `purpose == "search"` row (compare:28), membership in the three candidates (compare:45) and no duplicate key (compare:46). The search asserts R, workers, error flags, F0 sha and the seed derivation from the same row (search:64-69), and theta from the rds agrees with the ledger theta to 1e-13 (search:124). |
| missing rows counted as pass | not reachable | `missing == 0` and `nrow(rows) == 3` are both required (compare:51, 57). A task that fails before `saveRDS` leaves no rds; a task whose checks fail after saving exits 4 with `checks` FALSE and compare rejects it (compare:34). |
| doubles at reduced precision | not reachable | `identical()` on objective, EN_null, EN_alt against the rds doubles (compare:30), counts as integers (compare:28-29). Reference doubles are re-read by compare from the rds (compare:22-26), not taken from the output's `expected` field. See L1 for the one vacuous case of `identical()`. |
| compare reading the wrong directory | not reachable | Gate reads `$L3/out/$TAG` (gate:95, 97), the directory the tasks write (sbatch:126); submit refuses a used tag (submit:192) and sbatch refuses an existing output (sbatch:127). Compare and search both read `$L3_RESULTS`, whose sha256 were checked at staging against the desktop copies (stage:166, 176), and compare re-hashes both files against each output's recorded hashes (compare:31-32). |
| `afterany` compare passing when tasks failed | not reachable | A failed or timed-out task produces no rds (missing > 0) or a partial rds (`readRDS` error, compare exits non-zero). `--no-requeue` (sbatch:111). Gate exits with the compare status (gate:103). |
| stale outputs | not reachable | Fresh tag directory per submission, per-task `OUTF` refusal, `recursive = TRUE` listing only within the tag directory. |

## Low

### L1. compare: `identical()` accepts a NULL or zero-length pair (layer3_compare.R:30, 38-40)
- `dbl_ok` is `identical(x$objective, e$objective)` and likewise for EN_null and EN_alt, with no length-1 or finiteness
  check. If the wrapper result lacks a field (`m[["EN_null"]]` is NULL, search:254) and the stored table lacks the
  column (`h$EN_null` is NULL, compare:25), `identical(NULL, NULL)` is TRUE and the row reports "E0[N] identical".
- Failure scenario: a future candidate drawn from a run whose history has no EN_null/EN_alt (as the GA rds already
  has none), or an objective column renamed. The counts would still be checked, so the gate would not be fully
  vacuous, but the "doubles identical" claim would be.
- Not reachable for the three candidates: the local outputs hold length-1 numerics (objective 107.834..., 57.6915...,
  75.9132...; EN_null 80.913..., 71.472..., 76.804...), and the record (Section 1) re-read C2's stored E0[N] and E1[N].
- Fix direction: require `length(v) == 1L && is.finite(v)` on both sides before `identical()`.

### L2. compare: a zero-row ledger match or missing rds stops with an R error instead of the FAIL line (layer3_compare.R:20-27, 33-42)
- With no matching ledger row, `led_p` is `numeric(0)` and `data.frame(..., power_count_ledger = integer(0))` errors
  ("differing number of rows"); a missing ledger or rds file errors in `read.csv`/`readRDS`; a missing `counts` or
  `input_sha256` element errors on `[[`. R exits 1, so the gate fails closed, but without the summary line. This is
  layer 2's L1 in a new place; the original L1 (NA in `pass`) is fixed by the `isTRUE` and `%in%` coercions
  (compare:33-35, 47-48).

### L3. Library identity is a label, not a binary check (layer3.sbatch:125, layer3_stage.sh:174, layer3_search.R:101-104)
- Layer 3 reuses `$L2/lib` and checks only the `EngineCppSHA` DESCRIPTION line and the presence of
  `.simulate_chunk_cpp`. No hash of `evolveTrial/libs/evolveTrial.so` or of the installed R files is compared with what
  layer 2 passed with.
- Failure scenario: a manual `R CMD INSTALL` of another tree into `$L2/lib` carrying the same DESCRIPTION fields.
  Not reachable via the scripts: `layer2_stage.sh` is the only writer of `$L2/lib` and it refuses an existing install
  (layer 2 H1, now fixed at layer2_stage.sh:43, does not recur since layer 3 installs nothing).
- Fix direction: record `md5sum` of the installed `.so` at layer 2 staging and assert it in layer3_stage.sh and sbatch.

### L4. `L3_ENV` override is honoured from the submitting shell (layer3.sbatch:120, layer3_gate.sh:94, layer3_submit.sh:189)
- All three source `${L3_ENV:-.../harness/layer3_env.sh}`, and sbatch exports the submitter's environment by default.
  A stray `L3_ENV` would make the tasks and the gate agree on a different env file (engine sha, library, results
  directory), and compare would pass against that engine. Not reached for l3a: the evidence compare.csv shows
  `engine_sha` 1c918711... and outputs under `engine_cpp_l3/out/l3a`. Fix direction: drop the override, or have the
  gate assert that the sourced file's md5 equals the line in `$L3/harness.md5`.

### L5. Compare does not assert provenance fields (layer3_compare.R:33-42)
- `host`, `r_version` and `lib` are written to the CSV but not checked. A local dev rds copied by hand into
  `$L3/out/<tag>` would pass (it carries the same engine sha and identical numbers). Not reachable via the scripts
  (fresh tag, per-task refusal). Fix direction: require `lib == L3_LIB` and an `R version 4.5.0` prefix in the
  Longleaf gate.

### L6. run_local.sh masks failures and forwards `L3_FORCE_R` (run_local.sh:3-7)
- No `set -euo pipefail`; `Rscript ... | tee | cut` returns `cut`'s status, so a failed candidate does not stop the
  loop or the script. `${L3_FORCE_R:+...}` forwards a value exported in the caller's shell. Variables are unquoted.
  Not a false-PASS path: compare counts outputs and rejects `force_r` rows. Local-only driver.

### L7. Untagged call-path deviations in layer3_search.R (second question)
The diff against layer2_cert.R shows every numeric-path change under an L3 tag, but the following differences from
E2c's actual call path are not stated as deviations. None is numerically reachable, and bit identity on 3 of 3
confirms that for these candidates.
- a. Parse environment (search:78-82 then :133). Layer 3 parses the scenario with the BATON block set on BATON rows.
  E2c handed BATON `problem$wrapper`, which `f0_parse_multiarm` built with no BO_* set (problems.R:358-369,
  opt_baton.R:250-255); only the call ran inside `with_env(BATON block)`. In F0, BO_USE_BUFFER and
  BO_CONSTRAINT_BUFFER are read only in `scenario_configs.R:493, 502, 1318, 1322` for constraint bounds, which the
  wrapper does not use; the wrapper path reads only SLURM_CPUS_PER_TASK, MC_CORES (warmstart_wrappers.R:57, 67) and
  BATON_DISABLE_KSNAP (sim_fun:157). Not reachable.
- b. Load and source order (search:87-94, 131-133). Layer 3 attaches tidyverse and BATON, then sources
  scenario_parser.R, then warmstart_core.R, then parses. E2c sourced scenario_parser.R and parsed first, and loaded
  tidyverse, BATON (from R_LIBS, not `lib.loc`) and warmstart_core.R later (opt_baton.R:211-221). F0 code lives in
  globalenv, which masks attached packages either way, and the recon checked that warmstart_core.R redefines nothing
  on the wrapper path. Not reachable.
- c. Scenario selection by list name `sc[[SID]]` (search:134) rather than E2c's "exactly one `scenario_id` match"
  (problems.R:375-376). A duplicated name would pick the first silently, but search:136 checks `scenario_id`. Not
  reachable for this CSV (md5-pinned in f0_files.md5 line 42).
- d. Removed from layer 2's W without a tag: the atomic-theta coercion (layer2_cert.R `COERCE`) and
  `attributes(W) <- attributes(W0)`. Layer 3 always passes a list (search:115, 120), as GA's `problem$simulate`
  does (`as.list(theta)`, problems.R:427-431) and as BATON 0.4.0 does (`invoke_simulator` forwards theta unchanged;
  initial points are `as.list` of history rows). Nothing reads `attr(W, ...)` on this path. Not reachable.
- e. E2c's outermost `run_scenarios` echo wrapper (recon Section 2a) is not installed. It forces `base_args` once and
  draws no random numbers. Documented in recon Section 4 item 6, not in the script header.
- Fix direction: add one header line listing a to e as known, non-numeric deviations.

### L8. Temporary directory leak on "no free port" (layer3.sbatch:129-137)
Layer 2's L3 recurs: `TMPD` is created before the port search and the EXIT trap is set after it. Harmless.

## Layer 2 findings, carried over or fixed

| layer 2 | status in layer 3 |
|---|---|
| H1 (install without deps library) | not applicable: layer 3 installs nothing; fixed in layer2_stage.sh:43. Residual binary-identity gap is L3. |
| M1 (commit guard passes on untracked files) | fixed: stage:156 uses `git status --porcelain`, which lists untracked files. Probe: empty output today, all nine files tracked. |
| M2 (compare does not pin engine sha) | fixed: compare:14, 48, 57. |
| L1 (NA mode gives R error) | fixed for NA fields; recurs in another form as L2. |
| L2 (`read` on empty substitution) | not reachable: candidates come from a bash array with an index range check (sbatch:123-124). |
| L3 (TMPD leak) | recurs, L8. |
| L4 (`--wrap` under /bin/sh) | fixed: the gate is a `#!/bin/bash` script (gate:1). |
| L5 (misleading ssh message) | fixed in effect: stage:158 probes reachability first. |
| L6 (`rsync` without `--delete`) | not reachable: stage:159 refuses an existing harness or results directory; results go by `cp`. |
| L7 (untagged WT line) | resolved: search:35 now sits under the "paths (L2-b)" heading, which the header declares unchanged. New untagged items are L7 above. |

## Notes

- N1. The engine echo records `cfg$mode` (`mode=vsref` in out/local/vdev/PB-baton-1003_95.log), but the dev check
  (search:211-213) does not assert it. It adds nothing for numerical identity, since the counts and doubles are
  compared.
- N2. `out/local/` holds `vdev/` and `vref/` (dated 21:05, after `dev/` and `ref/`), which record Section 5 does not
  mention. Compare is never pointed at `out/local` as a whole. If it were, the duplicates would make it FAIL.
- N3. BUILD_REPORT.md line 3 says "nothing committed", but the harness is committed (0bba80f, edbbb83). Stale text.
- N4. Brief vs code: the brief describes compare as checking counts and doubles. It also checks input sha256
  (compare:31-32), mode, `force_r`, the nine capture checks and duplicates. Record decision 2 matches the code. No
  conflict that changes the verdict.

## Verdict

The gate's PASS for l3a (array 8072537, compare 8072538) is sound as far as a desktop review can establish. Every
false-PASS path in the brief is closed by a check in the search, the compare or the submission scripts. The
remaining items harden the gate for later reuse: L1 (NULL-safe `identical()`), L3 (binary hash of the reused
library) and L4 (env-file override). They do not reopen the l3a result. `layer3_search.R` changes nothing on the E2c
numeric call path beyond its tagged hunks. The five untagged differences in L7 are not numerically reachable and
should be listed in the header.
