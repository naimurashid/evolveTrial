# Step 7, readout: timing table and gate results for the bit-identical C++ engine (phase 1)

Engine under test 1c91871 on branch `perf/engine-cpp`, library `dev/perf/lib/dev_l2` (EngineCppSHA field), against the
reference library `dev/perf/lib/ref`. Local only, not pushed.

## 1. Timing table (desktop, 2026-11-03 21:27 to 21:34 EST)

Production path (`run_scenarios`, `parallel_replicates = TRUE`, `num_workers` w, FORK cluster made inside the call),
both scenarios of the captured call, 1,000 trials, seed 900001, `looks_stop_at_final = TRUE`. Each ref run was followed
directly by its dev run. Script and log in `dev/perf/evidence/step7_timing/` (`scripts/t7.R`, `scripts/run.sh`,
`run.log`, `timing_step7.csv`).

| design | workers | ref s | dev s | speed-up | output digest ref = dev | 1-minute load at start (ref, dev) |
|---|---|---|---|---|---|---|
| R01 (single-arm) | 1 | 26.3 | 17.6 | 1.50 | yes | 1.3, 2.9 |
| R01 | 16 | 5.2 | 3.3 | 1.60 | yes | 4.4, 6.3 |
| R23 (between-arm) | 1 | 235.8 | 25.4 | 9.29 | yes | 7.6, 13.7 |
| R23 | 16 | 32.2 | 3.4 | 9.38 | yes | 13.6, 21.3 |

The digest covers the `run_scenarios` result and `.Random.seed` after the call. The 1-worker and 16-worker digests
differ from each other, as they should, because each worker count draws from different streams. Within each pair ref
and dev are identical.

A baton-dosefinding job from another session (`n1t_search.R`, about 16 R processes) started on the desktop at about
21:28, so every row after the first ran on a loaded 20-core machine. The R01 1-worker row is the only clean one. It
matches step 3 (26.4 and 16.4 s, 1.59 times). The R23 1-worker ref time of 235.8 s is above step 4's best of 163.3 s
and inside its loaded range (163 to 230 s), so the R23 ratios are indicative. They agree with step 4's 8.7 times.

What the table shows: single-arm designs gain about 1.5 times, because about 64% of the time is gamma draws that bit
identity keeps (step 3). Between-arm designs gain about 9 times, and about 92% of the remaining C++ time is the 800 by K
gamma draws of the PH posterior (step 4). The design target of under 8 s per 1,000 trials for R23 is not met at 1
worker (25.4 s here, 18.75 s in step 4 for both scenarios, about 9.4 s per scenario). The ratio holds at 16 workers.

## 2. Gate results, all layers

| layer | what it checks | result | record |
|---|---|---|---|
| 0 | R01 and R23 end to end, ref vs dev, option FALSE and TRUE: table and `.Random.seed` | identical, 4 of 4 | step 4 commit 1c91871 |
| 1 | chunk-level C++ against the R body, all 10 designs | 2,120 chunks, 0 mismatches, with and without the cache | step 4 commit 1c91871 |
| package tests | evolveTrial test suite on dev | 756 expectations, 723 passed, 0 failed, 33 skipped | step 4 commit 1c91871 |
| 2 | certify capture path on Longleaf, 75 reference rows at R 10,000 | PASS 75 of 75, 7.12 core-h | step 5 record Section 6, `evidence/layer2_l2a/` |
| 3 | three E2c search evaluations on Longleaf at R 3,000 from recorded seeds | PASS 3 of 3, 0.15 core-h | step 6 record Section 6, `evidence/layer3_l3a/` |
| 3 local | the same three on the desktop, ref and dev, plus the FORCE_R negative control | PASS 3 of 3 each, control caught | step 6 record Section 5 |

Longleaf spend for steps 5 and 6 was 7.27 core-hours. Every gate in the plan passed. The C++ engine is bit identical
to the engine of record (7565a17 plus m4 patch 83de2c8b) on every path tested, from single chunks to whole E2c search
evaluations.

## 3. What stays open (carried from handoffs 03e to 03g)

1. P6 build field states lineage (step 5 record Section 1).
2. Review items L2, L3 and L7.
3. Codex review of the branch, after 2026-11-04 06:36.
4. Step 6 record Section 7 (review) is unfilled.
5. Questions 2 and 5 of handoff 2026-11-03e.
6. Phase 2 (statistically equivalent engine, draw count as fidelity) is the route to the R23 target. It needs the PI's
   decision and is outside phase 1.
7. Branch not pushed, engine not installed outside `dev/perf/lib`.
