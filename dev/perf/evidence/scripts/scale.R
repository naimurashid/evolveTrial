# scale.R -- parallel scaling of the captured call.
#   Rscript scale.R <design> <scen 1|2> <R> <workers> <mode internal|external>
# internal = production path (F0 wrapper): parallel_replicates TRUE, num_workers w, cluster_type "auto", no cluster
#            object, so run_simulation_pure makes and stops a FORK cluster inside every call.
# external = one FORK cluster from create_simulation_cluster(w, "FORK") passed as `cluster`; per-worker CPU read
#            with clusterCall(proc.time) before and after.
a <- commandArgs(TRUE); DES <- a[1]; SC <- as.integer(a[2]); R <- as.integer(a[3]); W <- as.integer(a[4]); MODE <- a[5]
S <- "/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf"
suppressPackageStartupMessages(library(evolveTrial))
stopifnot(startsWith(find.package("evolveTrial"), file.path(S, "Rlib")), W <= 16)
options(evolveTrial.looks_stop_at_final = TRUE)
x <- readRDS(file.path(S, sprintf("args_%s.rds", DES)))
ba <- x$base_args; ba$num_simulations <- R; ba$parallel_replicates <- TRUE; ba$num_workers <- W
stopifnot(identical(ba$min_follow_up_at_final, 6), identical(ba$random_dropout, FALSE))
scens <- x$scens[SC]
setup_s <- NA_real_; stop_s <- NA_real_; wcpu <- NA_real_
if (MODE == "external") {
  t0 <- proc.time()[["elapsed"]]; cl <- create_simulation_cluster(W, "FORK"); setup_s <- proc.time()[["elapsed"]] - t0
  ba$cluster <- cl
  cp <- function() sum(unlist(parallel::clusterCall(cl, function() { p <- proc.time(); p[["user.self"]] + p[["sys.self"]] })))
  c0 <- cp()
}
pt0 <- proc.time()
res <- run_scenarios(ba, scens, parallel = FALSE, seed = 900001, return_variance = TRUE)
pt <- proc.time() - pt0
if (MODE == "external") {
  wcpu <- cp() - c0
  t0 <- proc.time()[["elapsed"]]; parallel::stopCluster(cl); stop_s <- proc.time()[["elapsed"]] - t0
}
res <- as.data.frame(res)
row <- data.frame(design = DES, scen = SC, R = R, workers = W, mode = MODE, elapsed = pt[["elapsed"]],
                  master_cpu = pt[["user.self"]] + pt[["sys.self"]], child_cpu = pt[["user.child"]] + pt[["sys.child"]],
                  worker_cpu_clustercall = wcpu, cluster_setup_s = setup_s, cluster_stop_s = stop_s,
                  power_or_t1 = res$Type_I_Error_or_Power[nrow(res)], s_per_1000 = 1000 * pt[["elapsed"]] / R)
print(row)
saveRDS(row, file.path(S, "scale", sprintf("scale_%s_s%d_R%d_w%d_%s.rds", DES, SC, R, W, MODE)))
