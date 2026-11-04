# G5 statistical equivalence report

Verdict: **FAIL** (9 of 10 designs pass)

Input: CSV /home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/phase2_gate/selftest/pert_type1.csv with per-trial variances from /home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/phase2_gate/selftest/var.csv (row checks on mode and engine NOT applicable)
Engine identity: n/a (csv input)
Reference: /home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/acceptance_set_reference_rows.csv
Margins: power 0.010, type I 0.005, EN 1% of the pooled reference mean; each one-sided test at level 0.05 (z > 1.6449).

Each z pair is (z_lo, z_hi) of the two one-sided tests; a test rejects when its z exceeds the critical value.

| design | seeds | rows ok | power ref | power p2 | diff | z (lo, hi) | ok | type I ref | type I p2 | diff | z (lo, hi) | ok | EN0 ref | EN0 p2 | z (lo, hi) | ok | EN1 ref | EN1 p2 | z (lo, hi) | ok | class | design |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R01 | 5 | yes | 0.8214 | 0.8214 | +0.0000 | 4.1, 4.1 | yes | 0.0141 | 0.0241 | +0.0100 | 17.4, -5.8 | NO | 11.2223 | 11.2223 | 3.3, 3.3 | yes | 28.2242 | 28.2242 | 5.5, 5.5 | yes | yes | FAIL |
| R02 | 5 | yes | 0.8609 | 0.8609 | +0.0000 | 4.6, 4.6 | yes | 0.0154 | 0.0154 | +0.0000 | 6.4, 6.4 | yes | 13.2804 | 13.2804 | 4.2, 4.2 | yes | 27.6693 | 27.6693 | 7.6, 7.6 | yes | yes | PASS |
| R03 | 5 | yes | 0.9605 | 0.9605 | +0.0000 | 8.1, 8.1 | yes | 0.0123 | 0.0123 | +0.0000 | 7.2, 7.2 | yes | 21.4467 | 21.4467 | 8.2, 8.2 | yes | 29.8521 | 29.8521 | 44.6, 44.6 | yes | yes | PASS |
| R04 | 5 | yes | 0.8251 | 0.8251 | +0.0000 | 4.2, 4.2 | yes | 0.0057 | 0.0057 | +0.0000 | 10.5, 10.5 | yes | 13.6777 | 13.6777 | 4.8, 4.8 | yes | 28.0653 | 28.0653 | 6.7, 6.7 | yes | yes | PASS |
| R23 | 10 | yes | 0.8045 | 0.8045 | +0.0000 | 5.6, 5.6 | yes | 0.0824 | 0.0824 | +0.0000 | 4.1, 4.1 | yes | 59.2282 | 59.2282 | 7.6, 7.6 | yes | 91.1411 | 91.1411 | 14.1, 14.1 | yes | yes | PASS |
| R24 | 10 | yes | 0.8101 | 0.8101 | +0.0000 | 5.7, 5.7 | yes | 0.0863 | 0.0863 | +0.0000 | 4.0, 4.0 | yes | 63.8694 | 63.8694 | 11.4, 11.4 | yes | 83.7993 | 83.7993 | 22.2, 22.2 | yes | yes | PASS |
| R25 | 10 | yes | 0.8167 | 0.8167 | +0.0000 | 5.8, 5.8 | yes | 0.0937 | 0.0937 | +0.0000 | 3.8, 3.8 | yes | 70.8033 | 70.8033 | 16.0, 16.0 | yes | 81.5041 | 81.5041 | 31.0, 31.0 | yes | yes | PASS |
| R26 | 5 | yes | 0.8126 | 0.8126 | +0.0000 | 4.1, 4.1 | yes | 0.0853 | 0.0853 | +0.0000 | 2.8, 2.8 | yes | 60.7887 | 60.7887 | 6.4, 6.4 | yes | 86.0825 | 86.0825 | 13.6, 13.6 | yes | yes | PASS |
| R27 | 10 | yes | 0.8145 | 0.8145 | +0.0000 | 5.8, 5.8 | yes | 0.0935 | 0.0935 | +0.0000 | 3.8, 3.8 | yes | 74.2236 | 74.2236 | 20.5, 20.5 | yes | 82.1388 | 82.1388 | 31.8, 31.8 | yes | yes | PASS |
| R28 | 10 | yes | 0.8117 | 0.8117 | +0.0000 | 5.7, 5.7 | yes | 0.0949 | 0.0949 | +0.0000 | 3.8, 3.8 | yes | 66.8567 | 66.8567 | 14.0, 14.0 | yes | 79.6104 | 79.6104 | 27.1, 27.1 | yes | yes | PASS |

Failing designs: 
- R01: complete=TRUE rows_ok=TRUE power_ok=TRUE type1_ok=FALSE ENnull_ok=TRUE ENalt_ok=TRUE power_class_ok=TRUE type1_class_ok=TRUE

Classification sensitivity (reference pooled value within 3 standard errors of its threshold): 
- R23: power is 2.5 SE from 0.80, type I is -14.4 SE from 0.10
