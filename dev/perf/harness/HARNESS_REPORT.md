# Step 2b harness report

Worktree `~/research/evolveTrial-perf`, branch `perf/engine-cpp` (nothing committed). R 4.6.0, x86-64. Libraries used: `dev/perf/lib/ref` (agent A's build, DESCRIPTION has `VerifyM4PatchMD5: 83de2c8b...`) and `dev/perf/lib/dev` (this branch). Every install named `-l dev/perf/lib/dev`.

## Files

Package source: `R/simulation_driver.R` (seam, +7 -1 lines), `src/Makevars`, `src/r_semantics.h`, `src/r_semantics_exports.cpp`, `R/RcppExports.R` (+36), `src/RcppExports.cpp` (+106), `tests/testthat/test-r-semantics.R`.
Harness: `dev/perf/harness/` has `chunk_harness.R`, `layer0_seam_check.R`, `layer0_compare.R`, `run_layer1.R`, `run_pkg_tests.R`, `discrimination_check.R`, `tests/test-harness.R`, `out/` (logs, csv, layer 0 rds), `cache/` (gitignored). R 4.6.0 sources read: `dev/perf/build/rsrc/R-4.6.0/` (tarball plus extracted files, gitignored).

## 2b.1 Baseline (unmodified HEAD c2a88d5 in dev lib)

- Install 35 s. Package tests: 314 expectations passed, 0 failed, 33 skipped (6 files: hybrid-ba-margin, parallel-execution, path-parity, rcpp-equivalence, rcpp-r-equivalence, trial-modes), 0 errors, 0 warnings, 35 s.
- After step 2b (seam, Makevars, helpers, new test file): 347 passed, 0 failed, 33 skipped, 0 errors, 28 s. The 33 added passes are `test-r-semantics.R`. Nothing got worse.

## 2b.2 Seam

Rename `simulate_chunk` to `simulate_chunk_r` (line 468), then after `aggregate_results` the five lines of the plan (now lines 946-950). No other `simulate_chunk` reference existed before line 946. No random draw happens before the probe: `chunk_frame` leaves `.Random.seed` identical (test), and the parent stream after `run_scenarios` equals the stream the chunk plan reproduces (H3).

## 2b.3 Layer 0 (R = 200, 4 workers, seed 900001, `return_variance = TRUE`, `looks_stop_at_final` option TRUE)

| Design | ref vs unmodified HEAD (before seam) | ref vs dev, seam (after seam) | ref vs dev, final lib (with Makevars and helpers) |
|---|---|---|---|
| R01 (2 rows x 23 cols) | table TRUE, `.Random.seed` TRUE | TRUE, TRUE | TRUE, TRUE |
| R23 (4 rows x 23 cols) | TRUE, TRUE | TRUE, TRUE | TRUE, TRUE |

Result: `identical()` TRUE on the table and on `.Random.seed` for both designs in every comparison (12 of 12). The first column also shows that the reference library equals the unmodified branch. Layer 0 must be rerun after any later R-code change (`layer0_seam_check.R` then `layer0_compare.R`). Args used: `dev/perf/evidence/args_R01.rds`, `args_R23.rds` (the early K-snapped captures).

## 2b.4 and 2b.5 Harness and self-tests

`Rscript -e 'testthat::test_file("dev/perf/harness/tests/test-harness.R")'` against dev lib: 26 expectations passed, 0 failed, 0 skipped.

| Test | Expectations | Result |
|---|---|---|
| H6 option guard (`chunk_frame` errors when option unset) | 1 | pass |
| Frame: RNG untouched, `looks_stop_at_final` TRUE, `simulate_chunk_r` present, R recorded | 4 | pass |
| H1 reference vs itself, R01 and R23, 20 trials | 4 | pass |
| H2 composition, trial by trial vs chunk, 30 trials, R01 and R23 | 4 | pass (stop condition not triggered) |
| H3 plan emulation, R01, R = 200, 4 workers | 3 | pass: 40 of 40 summary columns identical (2 scenarios x 20 columns incl. all `Var_*`), parent `.Random.seed` equal. Also run for R23: 40 of 40, stream equal |
| H4 stream fault at trial 7 of 20 | 3 | pass: k = 7, kind "rng" |
| H5 value fault at trial 12 of 20 | 4 | pass: k = 12, kind "result", per-trial sum delta differs by exactly 1 |
| Reference cache (miss, hit, key changes with seed) | 4 | pass |

Layer 1 self-check (`run_layer1.R --designs R01,R23 --cand simulate_chunk_r --no-cache`, R = 50, workers 16 and 19, masters 900001-900003, plus a 200-trial chunk per design and scenario): 424 chunks compared (212 per design: 2 scenarios x 3 masters x (16 + 19) chunks = 210, plus 2 long chunks), 0 mismatches, 187 s. With the cache, R01 alone: 212 chunks, 0 mismatches, 27 s cold and 14 s warm (212 cache files). A faulty candidate (extra `runif(1)` before sim index 3) run through `run_layer1`: 17 chunks, 2 mismatches (chunk 1 and the long chunk), both bisected to k = 3, kind "rng". The CLI writes `out/layer1_<timestamp>.csv` and exits 1 on a mismatch. The full 10-design run was not made (agent A's args set was still being checked; all ten files `args_R01..R28.rds` are in `dev/perf/evidence/args/`, and the runner picks them up by default).

## 2b.6 and 2b.7 Floating-point helpers and unit tests

`test-r-semantics.R`: 33 expectations passed, 0 failed (7 s). Cases are `identical()` comparisons.

| Helper (wrapper) | Copied from (R 4.6.0) | Cases |
|---|---|---|
| `sum_ld` (`rsem_sum`) | `src/main/summary.c` `rsum()`, `do_summary()` (`zcum.r += tmp` from 0) | 2,002 inputs (2,000 seeded, two 1e5 vectors) + 8 NA/Inf/overflow vectors |
| `mean_ld` (`rsem_mean`) | `summary.c` `real_mean()` | 2,002 + 9 special (NA, NaN, Inf, overflow, empty) |
| `median_r` (`rsem_median`) | `stats:::median.default`, `src/library/stats/R/median.R` | 2,002 + 8 special |
| `quantile7`, `quantile7_sorted` (`rsem_quantile7`) | `stats:::quantile.default`, `quantile.R`, type 7 | 2,000 (p random or from 10 fixed values) + 300 full `seq` probability vectors K = 3..16 (about 2,900 quantiles) + empty input |
| `find_interval` (`rsem_find_interval`) | `src/appl/interv.c` `findInterval2()` via `src/main/util.c` `do_findinterval()` | 2,000 cases x about 40 points, all 8 flag combinations, 0 to 17 breaks, ties in breaks, points on and around breaks, outside, NaN; plus the engine's flag set |
| `mean_logical` (`rsem_mean_logical`) | `summary.c` `logical_mean()` | all (count, n) for n = 1..2,000 (2,003,000 cases incl. n + 1 counts), n = 0, plus every count for n in 800, 3511, 5000, 5998, 7107, 8846, 10000 (41,269 cases) |
| `seq_probs` (`rsem_seq_probs`) | `base:::seq.default`, `src/library/base/R/seq.R` | K = 1..40, with `K + 1` double and `K + 1L` integer |
| `rebalance_cutpoints` (`rsem_rebalance_cutpoints`) | `R/simulation_driver.R:1331-1398` built from the above | 2,000 random (length 50-300, Weibull 1.3, medians 3.8-9, some beyond 120, ties, zeros, values at the cap, cap 24 or 120, K 3..16) + 9 degenerate |

Plus a direct long double test (`c(1, rep(1e-16, 10))` separates long double from double accumulation, and `rsem_ld_digits() == 64`). `discrimination_check.R` shows the tests can fail: on the seeded inputs a naive double accumulator differs from `sum()` on 1,540 of 2,000 inputs and `sum/n` from `mean()` on 1,574 of 2,000; `i/K` differs from `seq()` for 31 of 40 K.

`src/Makevars` is `PKG_CXXFLAGS = -ffp-contract=off` (confirmed on every compile line). `compileAttributes()` was first run on the unmodified tree and changed nothing, so the Rcpp 1.1.2 regeneration is clean. After adding the helpers the diff of `R/RcppExports.R` and `src/RcppExports.cpp` is additions only (142 insertions, 0 deletions), and all eight new wrappers use `rng = false` (no `RNGScope` in the diff). The helpers are internal (not in NAMESPACE): tests reach them with `getFromNamespace`.

## Conflicts between the plan and the code

1. **Baseline command.** `test_dir` with a relative `lib.loc` fails (testthat changes directory, so lazy-load files are not found: 88 errors on the first attempt). `run_pkg_tests.R` uses `normalizePath(lib)`. The same applies to the coordinator's verification command.
2. **`R CMD INSTALL --preclean` deletes the tracked file `src/symbols.rds`.** I ran `git checkout src/symbols.rds` after each install. Check `git status` for ` D src/symbols.rds` before any commit.
3. **Bisection `kind`.** The plan lists "rng" | "result" | "both" and expects "rng" in H4. An extra draw before trial 7 also changes trial 7's own values (verified: events 4 vs 22, final N 15 vs 32), so a three-way split would report "both". I used: "rng" whenever the stream state differs at the first diverging trial, "result" when only values differ, and added `rng_differs` and `result_differs` booleans so nothing is lost. "both" is never returned.
4. **`find_interval`** calls R's own `findInterval2()` (declared in `R_ext/Utils.h`, present in R 4.5.0 and 4.6.0) instead of a copy of `interv.c`. The plan asks for a copy. The routine is R's by construction, and no GPL source is pasted. The ilo hint does not change the result for sorted breaks; NaN gives NA as in `do_findinterval`. Engine flags (`state_management.R:27-31`): `left.open = FALSE`, `rightmost.closed = FALSE`, `all.inside` default FALSE.
5. **`mean_logical`.** The plan and inventory 4.36 say `count / n` differs from `mean()` in about 1 of 8,000 cases. For n up to 2,000 it never differs (0 of 2,003,000); it first differs at n = 3,511 (count 2,613) and occasionally above (about 330 of 2e6 random pairs with n up to 1e4). The plan's "n up to 2,000" test therefore cannot tell the two apart, so I added the larger-n test above. The engine's posterior draws (800 to 10,000) are in the range where it matters.
6. **`seq_probs`** is `0 + i * (del / n1)` with `del = 1`, `n1 = K`, last point exactly 1 (from `seq.default`), the same as `i * (1.0 / K)` for K >= 2; K = 1 and K = 0 follow the shorter branches.
7. **`pmax`/`pmin` in `rebalance_cutpoints`** follow `do_pmin` (replace when `tmp > ra` or `tmp < ra`), which keeps the sign of a `-0`. `length(event_times) < max(2, K)` counts NA entries before `na.rm` removes them, as R does.
8. **Cache key.** The plan keys on the deparsed body, the serialized args, indices and seed. I added `looks_stop_at_final`, `EVOLVETRIAL_USE_CPP`, the R version and the md5 of the package source files the R body reaches (all `R/*.R` except `simulation_driver.R` and `RcppExports.R`; `src/*.cpp|h` except `RcppExports`, `r_semantics*`, `simulate_chunk*`), so a change to the posterior or state code cannot return a stale reference. Delete `dev/perf/harness/cache/` to force a rerun.
9. **Port collisions.** While agent A runs, `makeForkCluster` sometimes fails with "port 11434 cannot be opened" (random port clash). The layer 0 and H3 runs were repeated until they went through; no result differed between attempts.
10. Layer 0 and the tests use the early `args_R01.rds` and `args_R23.rds` (K-snapped, R23 has 12 intervals), as the plan allows; the recaptured args from step 2a have not been run through layer 1.
