# Layer 0: the seam does not change the R path. (plan 2b.3, breakage B1)
#   Rscript dev/perf/harness/layer0_seam_check.R <lib> <args.rds> <out.rds>
# Runs run_scenarios on a captured design at R = 200, 4 workers, seed 900001, return_variance = TRUE, and saves the
# returned table together with .Random.seed. Run it against dev/perf/lib/ref and dev/perf/lib/dev, then compare the
# two output files with identical() (see layer0_compare.R).
a <- commandArgs(TRUE)
if (length(a) != 3L) stop("usage: Rscript layer0_seam_check.R <lib> <args.rds> <out.rds>")
lib <- normalizePath(a[1]); args_rds <- a[2]; out_rds <- a[3]
stopifnot(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
suppressPackageStartupMessages(library(evolveTrial, lib.loc = lib))
stopifnot(dirname(getNamespaceInfo("evolveTrial", "path")) == lib)
options(evolveTrial.looks_stop_at_final = TRUE)
d <- readRDS(args_rds)
base_args <- d$base_args
base_args$num_simulations <- 200L
base_args$num_workers <- 4L
t0 <- Sys.time()
res <- evolveTrial::run_scenarios(base_args, d$scens, parallel = FALSE, seed = 900001, return_variance = TRUE)
saveRDS(list(res = res, rng = .Random.seed, lib = lib, design = d$design,
             wall = as.numeric(difftime(Sys.time(), t0, units = "secs")),
             rver = as.character(getRversion())), out_rds)
cat(sprintf("layer0 %s lib=%s rows=%d wall=%.1fs\n", d$design, basename(lib), nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
