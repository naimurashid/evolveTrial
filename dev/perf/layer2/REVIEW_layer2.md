# Review of the layer 2 harness (dev/perf/layer2), 2026-11-03

Scope: correctness defects only. Read-only review: diff, grep, line ranges, one `git diff --quiet` probe. Nothing was
run on Longleaf, so every Longleaf-side finding is UNCONFIRMED unless it says otherwise. Not redone: the local 10-row
pass, the L2_FORCE_R failure test, and the compare negative tests (design record Section 1).

Counts: Critical 0, High 1, Medium 2, Low 7, plus 3 notes on plan conformance.

## High

### H1. layer2_stage.sh:183 installs the engine without the deps library on the library path
- What: the remote heredoc runs `R CMD INSTALL --no-multiarch -l "$L2_LIB" "$L2_TREE"` after `l2_modules`. `R_LIBS`
  is set only on the later `Rscript -e` line (185), and `R_LIBS_USER` keeps its login default. evolveTrial at 1c91871
  has `Imports: data.table, ggplot2, magrittr, progress, Rcpp, ...` and `LinkingTo: Rcpp, RcppArmadillo`.
- Failure scenario: (a) if the user's default R 4.5 library on Longleaf lacks these packages, INSTALL stops with
  "dependencies ... are not available" and the stage aborts (fails closed, but layer 2 cannot be staged as written);
  (b) if the default user library has them at other versions, the engine compiles against Rcpp/RcppArmadillo headers
  that differ from those in `R_libs_verify_m4`, which the run then loads at runtime. Case (b) produces a build that
  the load check (lines 185-193) accepts, and whose provenance is not the certify library.
- Fix direction (not applied): run the INSTALL with `R_LIBS="$L2_DEPS_LIB" R_LIBS_USER=/dev/null` (or the same
  pair as the load check), and record `packageVersion("Rcpp")`/`"RcppArmadillo"` in the load check output.
- Verdict: UNCONFIRMED (needs Longleaf). The code path is certain; which of (a) or (b) happens depends on
  `~/R/x86_64-pc-linux-gnu-library/4.5` on Longleaf.

## Medium

### M1. layer2_stage.sh:156 "commit first" guard passes on untracked files
- `git diff --quiet HEAD -- dev/perf/layer2 ...` ignores untracked files. Right now `git ls-files dev/perf/layer2`
  lists 0 files and the command exits 0.
- Failure scenario: stage is run before the commit (the current state). layer2_cert.R and the reference rows are still
  pinned by md5 (lines 158-159), but layer2_compare.R, layer2.sbatch, layer2_submit.sh, layer2_env.sh and
  build_dev_tree.sh go to Longleaf with no commit that records them. The gate script that decides PASS is then not
  reproducible from the repository.
- Fix direction: also require `git ls-files --error-unmatch` on each HF file plus build_dev_tree.sh, or
  `git status --porcelain -- <paths>` empty.
- Verdict: CONFIRMED (exit 0 with zero tracked files, 2026-11-03).

### M2. layer2_compare.R does not pin the engine commit
- Lines 38 and 45 require `engine_sha` to be unique and non-NA, not equal to `L2_ENGINE_SHA`. layer2_cert.R enforces
  equality with the environment variable at run time, and layer2.sbatch passes the env.sh value, so on the intended
  path this holds. Failure scenario: outputs from a run staged at another commit (env.sh edited, new tag) are compared
  with no statement of which engine passed; the gate line prints the sha but does not check it.
- Fix direction: pass the expected sha (argument or env) and require equality.
- Verdict: UNCONFIRMED as a live failure (needs a mis-staged run); the missing assertion is certain.

## Low

- L1. layer2_compare.R:35,43,45. A row with `mode` NA (an RDS without the field) makes `m$pass` NA; `any(!m$pass)` and
  `if (ok)` then error, and R exits 1. Fails closed, but with an R error instead of the FAIL line. `isTRUE`-style
  coercion on `mode == "dev"` would give a clean FAIL. Not a pass-on-wrong-input path.
- L2. layer2.sbatch:122-123. Under `set -e`, `read` on an empty process substitution (awk `exit 3` on a bad header, or
  no row i) returns 1 and the script exits before the "no reference row" message. Fails closed, silent.
- L3. layer2.sbatch:127-135. TMPD is created before the port search, and the EXIT trap is set only after a port is
  found, so "no free port" leaks `/tmp/l2cpp_*`. Harmless.
- L4. layer2_submit.sh:213. `--wrap` scripts run under `/bin/sh`; layer2_env.sh uses `${L2_ENGINE_SHA:0:12}` (bash
  substring). On Longleaf (RHEL) `/bin/sh` is bash, so this works; on a dash system the compare job would fail on
  source. Quoting of the wrap string itself is correct: `\$O`, `\$2`, `\$3`, `\":\"` and `\\n` arrive as intended.
- L5. layer2_stage.sh:160. Any ssh failure (key, network) reports "refusing: $L2_LIB/evolveTrial exists". Misleading
  message only.
- L6. layer2_stage.sh:172. `rsync -a` without `--delete` into an existing `$L2/src/<tree>` would keep stale files; the
  md5 check covers listed files only. Not reachable on a first stage (line 160 refuses an existing install, though not
  an existing tree directory).
- L7. layer2_cert.R:56. `normalizePath(..., mustWork = FALSE)` for WT is an untagged line inside the L2-b hunk. It is
  needed on Longleaf (no `~/research` there) and does not touch the call path. Tag it for the header's "nothing else"
  claim.

## Question-by-question

### 1. capture_cert.R vs layer2_cert.R
Every hunk of the diff is L2-tagged except L7 above. The changes that touch what produces the numbers: LIB (the chosen
library, L2-b), CSEED (taken from the argument and checked to be one of the design's cert seeds, then the same
`cr` row lookup; the task seed `SEED` for BO_SEED is unchanged), and the `.simulate_chunk_cpp` wrapper. The ENV block,
the F0 source order, W0, base_args, the capture hook order (innermost, before M4 and P6), `W(th, fidelity = "high",
seed = CSEED)` and `counts()` are byte-identical to capture_cert.R. The extra work in W (`ncp`, `ncap`, the cpp echo)
reads files and lists and draws no random numbers. No defect.

### 2. Dev engine check (L2-e)
- Partial R fallback: any R-body trial is not counted by the wrapper, so `cp_trials < nscen * Rf` and the check fails.
  This positive count is the guard that carries the check. The two negative terms (`length(pr) == 0`,
  `NROW(ev) == 0`) are vacuous if the probe directory is missing or `p6r_read_events` errors (NULL gives NROW 0), but
  that cannot let an R trial through because the count would still be short. No defect.
- Wrong library: refused unless the DESCRIPTION carries P6 7565a17, `EngineCppSHA == L2_ENGINE_SHA`, no M4 patch
  field, and the namespace path equals LIB (checked twice, before and after BATON loads). A tree installed from a
  modified source with the same fields would pass; the stage md5 check on the tree is the only defence (and see H1 for
  the headers it compiles against).
- Wrong worker count: `SLURM_CPUS_PER_TASK` must equal 19 (sbatch:118 and cert), and the capture check asserts
  `num_workers == 19` in base_args. The cluster size is derived from that argument (simulation_driver.R:971-978).
  `cp_pids` is echoed but not asserted; asserting it would add nothing for numerical identity.
- Stale probe files: CPPD is under `OUT = tempfile()` inside a per-job `TMPDIR`, and the check uses the delta since the
  call began. Per-pid files avoid interleaved appends. No defect.
- Valid input wrongly refused: none found. `format(6)` is "6" for 6 and 6L; FALSE/TRUE formats are fixed. Zero-length
  chunks add 0.
- RNG and results: the wrapper calls `Sys.getpid`, `sprintf`, `format`, `cat` and then the original closure with the
  same arguments; none touches `.Random.seed`. The dispatcher (simulation_driver.R:949-953) resolves
  `.simulate_chunk_cpp` in the namespace at call time, so FORK workers created after installation use the wrapper. The
  wrapper cannot change results.

### 3. layer2_compare.R
Wrong seeds (not in ref), duplicates, ref-mode rows, force_r rows, failed checks, R != 10000, and missing rows all fail,
and the full gate needs `nrow(m) == 75`, so no vacuous pass was found (NA gives an error exit, L1). On precision: the
reference CSV is write.csv output at 15 significant digits (weighted_EN_N shows 13 decimals), so its stored precision is
not 4 decimals. EN_null and EN_alt have at most 4 decimals in all 75 rows because they are means of integers over
10,000 trials. A difference in the true value is therefore at least 1e-4 and the 4-decimal comparison catches it, so it
is adequate. Comparing at full precision (or also comparing `weighted_EN_N`, which carries the extra digits) would be
the stricter reading at no cost. Suggestion, not a defect.

### 4. Shell
- layer2_stage.sh heredoc: unquoted `EOF`, so `$L2_HARNESS` expands on the desktop (same value as on Longleaf, fine);
  every remote variable and `$(...)` is escaped, and the R code's `d\$Field` arrives as `d$Field`. No expansion bug.
  The defect is H1, not quoting.
- layer2_submit.sh `--wrap`: correct (L4 only).
- layer2.sbatch awk: the CSV has no quotes (0 `"` characters) and 20 fields on every line, column 7 is `seed`, and the
  header test guards both column 1 and column 7. Correct.
- Nothing found that would make all 75 tasks fail, apart from H1 preventing the stage (and therefore the submit,
  whose library check at submit.sh:207 would refuse).

### 5. Plan conformance (engine-cpp-plan.md Step 5, inventory I8)
- N1. Estimate: the plan says 3 to 5 core-hours; the design record estimates 7.8 (13.4 worst case). The record
  explains the gap. The plan requires the coordinator to state the estimate to the PI before submitting, so this is a
  required disclosure, not a defect.
- N2. The plan says "install lib/dev on Longleaf"; the harness installs a `git archive` tree with two DESCRIPTION
  fields appended. The record (Section 1, item 1) documents the reason (the frozen P6 protocol refuses a library
  without the P6 field). Benign, documented deviation; the `P6RecertEvolveTrialSHA: 7565a17` field on a non-7565a17
  build is a lineage claim that the PI should know is being made.
- N3. I8 (19 workers, asserted): satisfied by sbatch:118, the cert's SLURM check, and the `num_workers_19` capture
  check.
