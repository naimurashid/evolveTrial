# counts.R -- per-trial call counts (interim looks, posterior sampler calls, registry size) and, for R23, an
# exploratory in-memory swap of the shadowed C++ dispatcher for sample_vs_ref_medians_ph (scratch session only).
#   Rscript counts.R <design> <R> [swap]
a <- commandArgs(TRUE); DES <- a[1]; R <- as.integer(a[2]); SWAP <- identical(a[3], "swap")
S <- "/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf"
suppressPackageStartupMessages(library(evolveTrial))
stopifnot(startsWith(find.package("evolveTrial"), file.path(S, "Rlib")))
options(evolveTrial.looks_stop_at_final = TRUE)
x <- readRDS(file.path(S, sprintf("args_%s.rds", DES)))
ba <- x$base_args; ba$num_simulations <- R; ba$parallel_replicates <- FALSE; ba$num_workers <- 1L
ns <- asNamespace("evolveTrial")
if (SWAP) {
  # the dispatcher body from posterior_cpp_dispatchers.R:241, which posterior_helpers.R:105 overrides at collation
  src <- readLines(file.path(S, "et58/R/posterior_cpp_dispatchers.R"))
  i0 <- grep("^sample_vs_ref_medians_ph <- function", src); i1 <- grep("^sample_vs_ref_medians_ph_r <- function", src) - 1
  f <- eval(parse(text = src[i0:i1])[[1]]); environment(f) <- ns
  utils::assignInNamespace("sample_vs_ref_medians_ph", f, ns = "evolveTrial")
  pt0 <- proc.time(); res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE); pt <- proc.time() - pt0
  cat(sprintf("SWAP %s R=%d x %d: elapsed %.2f s; %.2f s per 1000 trials\n", DES, R, length(x$scens), pt[["elapsed"]], 1000 * pt[["elapsed"]] / (R * length(x$scens))))
  print(as.data.frame(res)[, c("scenario", "Arm_Name", "Type_I_Error_or_Power", "PET_Futility")])
  ref <- readRDS(file.path(S, sprintf("timing_seq_%s_R%d.rds", DES, R)))$res
  cat("identical OC to the shipped (R-level) sampler at the same seed:",
      isTRUE(all.equal(ref$Type_I_Error_or_Power, as.data.frame(res)$Type_I_Error_or_Power, tolerance = 0)), "\n")
  quit()
}
cnt <- new.env(); cnt$look <- 0L; cnt$post <- 0L; cnt$slice <- 0L; cnt$ee <- 0L
wrap <- function(name, key) { o <- get(name, envir = ns); force(o)
  utils::assignInNamespace(name, function(...) { cnt[[key]] <- cnt[[key]] + 1L; o(...) }, ns = "evolveTrial") }
wrap("interim_check", "look"); wrap("slice_arm_data_at_time", "slice")
if (isTRUE(ba$compare_arms_option)) wrap("sample_vs_ref_medians_ph", "post") else wrap("draw_posterior_hazard_samples", "post")
res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE)
nt <- R * length(x$scens)
cat(sprintf("%s per trial (R=%d x %d scen): interim_check %.1f, slice_arm_data_at_time %.1f, posterior sampler calls %.1f (x %d draws each)\n",
            DES, R, length(x$scens), cnt$look / nt, cnt$slice / nt, cnt$post / nt, ba$num_posterior_draws))
cat("max_total_patients_per_arm:", ba$max_total_patients_per_arm, " registry preallocated rows per arm = that value\n")
