# P3-2: golden fixture update (copy of p31_golden.R). Recomputes all 18 golden entries (the 12 existing and the 6 P3-2 entries of
# tests/testthat/test-engine-equivalence.R) with the given build, checks that the 12 existing entries are identical()
# (table and RNG, no tolerance) to the stored ones, and writes the fixture with the 12 stored objects copied verbatim plus
# the 6 new ones. Usage: Rscript dev/perf/phase3/p32_golden.R <lib> [--write]
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
library(evolveTrial, lib.loc = normalizePath(a[1]))
library(testthat)
src <- readLines("tests/testthat/test-engine-equivalence.R")
# take the helpers and golden_specs/golden_run from the test file, nothing else
grab <- function(start, end_pat) { i <- grep(start, src, fixed = TRUE)[1]; j <- i - 1L + grep(end_pat, src[i:length(src)])[1]; src[i:j] }
eval(parse(text = c(grab("fx <- function(D)", "^fx"), grab("golden_path <- function()", "^golden_path"),
                    grab("run_design <- function(", "^}"), grab("golden_specs <- list(", "^)$"),
                    grab("golden_run <- function(g)", "base_mod = if"))))
fx <- function(D) readRDS(file.path("tests/testthat/fixtures", paste0("args_", D, ".rds")))
gpath <- "tests/testthat/fixtures/engine_golden.rds"
old <- readRDS(gpath)
cat("stored entries:", paste(names(old), collapse = ", "), "\n")
now <- lapply(golden_specs, golden_run)
for (k in names(old)) cat(sprintf("%-14s tab identical %s, rng identical %s\n", k, identical(now[[k]]$tab, old[[k]]$tab),
                                  identical(now[[k]]$rng, old[[k]]$rng)))
new <- c(old, now[setdiff(names(golden_specs), names(old))])[names(golden_specs)]
for (k in setdiff(names(golden_specs), names(old))) cat(sprintf("%-14s new: Exp_N %s, power/T1 %s, NA %s\n", k,
  paste(round(now[[k]]$tab$Exp_N, 3), collapse = " "), paste(round(now[[k]]$tab$Type_I_Error_or_Power, 4), collapse = " "),
  anyNA(now[[k]]$tab)))
if ("--write" %in% a) {
  stopifnot(all(vapply(names(old), function(k) identical(now[[k]]$tab, old[[k]]$tab) && identical(now[[k]]$rng, old[[k]]$rng), TRUE)))
  saveRDS(new, gpath, version = 2)
  chk <- readRDS(gpath)
  cat("written; existing entries identical to the previous file:",
      all(vapply(names(old), function(k) identical(chk[[k]], old[[k]]), TRUE)), "\n")
}
