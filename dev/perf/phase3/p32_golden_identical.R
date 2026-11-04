# P3-2: all 18 golden entries reproduced bit for bit (identical(), no tolerance) by a build. Usage: Rscript ... <lib>
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
library(evolveTrial, lib.loc = normalizePath(a[1]))
src <- readLines("tests/testthat/test-engine-equivalence.R")
grab <- function(start, end_pat) { i <- grep(start, src, fixed = TRUE)[1]; j <- i - 1L + grep(end_pat, src[i:length(src)])[1]; src[i:j] }
eval(parse(text = c(grab("run_design <- function(", "^}"), grab("golden_specs <- list(", "^)$"),
                    grab("golden_run <- function(g)", "base_mod = if"))))
fx <- function(D) readRDS(file.path("tests/testthat/fixtures", paste0("args_", D, ".rds")))
gold <- readRDS("tests/testthat/fixtures/engine_golden.rds")
before <- readRDS("dev/perf/phase3/p32_engine_golden_before.rds")
cat("stored file: existing 12 entries identical to the pre-P3-2 file:",
    all(vapply(names(before), function(k) identical(gold[[k]], before[[k]]), TRUE)), "\n")
for (k in names(golden_specs)) {
  now <- golden_run(golden_specs[[k]])
  cat(sprintf("%-14s tab identical %s, rng identical %s\n", k, identical(now$tab, gold[[k]]$tab), identical(now$rng, gold[[k]]$rng)))
}
