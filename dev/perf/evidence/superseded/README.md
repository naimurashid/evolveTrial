# Superseded captures

`args_R01` and `args_R23` were captured by `scripts/capture.R` on 2026-11-03 without the certify job's environment
block. `BATON_DISABLE_KSNAP` was unset, so R23 was snapped from 15 intervals to 12 (plan inventory I1). They were
used only for the profiling in `profile_report.md`, which measures time, not results. The gate inputs are the
recaptures in `../args/`, which reproduce their certification rows exactly (`../capture_report.md`). Do not use these
two files for any identity test.
