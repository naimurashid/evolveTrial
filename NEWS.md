# evolveTrial 0.0.0.9000

## Package Maintenance

* Fixed R CMD check warnings and notes for CRAN compliance:
  - Removed non-ASCII characters (emojis, Greek letters, smart quotes) from source code
  - Replaced `library()` calls with `requireNamespace()` + `::` for proper package dependencies
  - Reorganized DESCRIPTION: moved optional packages (gt, dplyr, openxlsx, etc.) to Suggests
  - Removed unused imports (ggpubr, survival, survminer)
  - Added VignetteBuilder field for proper vignette support

* Documentation improvements:
  - Added complete `@param` documentation for internal helper functions
  - Fixed documentation mismatches in `explore_early_stopping_from_cal()`
  - Added `@keywords internal` tags for non-exported functions

* Code quality enhancements:
  - Added `stats::rnorm` to namespace imports
  - Added missing global variable declarations
  - Updated .Rbuildignore to exclude development files
  - Improved test coverage with actual utility function tests

## Bug Fixes

* None in this development version

## New Features

* None in this development version (package in active development)
