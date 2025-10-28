# Repository Guidelines

## Project Structure & Module Organization
Core simulation and decision logic lives in `R/` (`state_management.R`, `interim_logic.R`, `simulation_driver.R`). Statistical helpers such as posterior samplers and predictive calculations sit alongside them. Tests are under `tests/testthat`, documentation topics in `man/` with longer-form walkthroughs in `vignettes/`. Artifacts for pkgdown live in `_pkgdown.yml`, while `renv/` pins package dependencies for reproducible runs. Keep auxiliary scripts (e.g., scenario grids) in the project root or `inst/` so they do not pollute the package namespace.

## Build, Test, and Development Commands
Use `devtools::load_all()` during iterative work; it refreshes the package namespace without reinstalling. Run `devtools::test()` for targeted validation and `devtools::check()` before publishing to ensure R CMD check passes cleanly. For local script runs, `Rscript path/to/file.R` keeps environments isolated. After major updates, `devtools::install()` followed by `library(evolveTrial)` mirrors the user experience. `pkgdown::build_site()` regenerates browsable documentation once substantive API changes land.

## Coding Style & Naming Conventions
Follow tidyverse-style R code: two-space indentation, `{` on the same line, and `snake_case` identifiers (`min_person_time_frac_per_arm`). Reuse package utilities such as `%||%` and `coalesce_num()` instead of re-implementing fallbacks. Keep diagnostic `message()` calls concise and prefixed (e.g., `[vsREF gate]`). When adding arguments, thread them through `run_simulation_pure()` formals so downstream helpers can pick them up via `modifyList()`.

## Testing Guidelines
New decision logic or gating tweaks must include regression tests under `tests/testthat/`. Prioritize scenarios that stress multi-arm comparisons, proportional gate scaling, and predictive probability fallbacks. Use deterministic seeds (e.g., `set.seed(4242)`) and low simulation counts for unit tests; save large grids for ad-hoc scripts. Document required outputs (PET, alpha, expected N) in expectations so future refactors flag behavior shifts early.

## Commit & Pull Request Guidelines
Write imperative, present-tense commit titles under ~72 characters (see `git log`: “Add one-time interval rebalancing support”). Group related edits—code, docs, and tests—into a single commit when feasible. PRs should summarize design intent, note any simulation diagnostics (attach the key rows from `res[...]`), and call out configuration changes affecting calibration. Link to upstream issues or design discussions so future agents understand the rationale. Before requesting review, ensure `devtools::check()` is clean and `results.txt` or other intermediate files are ignored or removed.
