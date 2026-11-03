# Layer 0: the seam does not change the R path. (plan 2b.3, breakage B1)
#   Rscript dev/perf/harness/layer0_seam_check.R <lib> <args.rds> <out.rds> [use_cpp_engine]
# Runs run_scenarios on a captured design at R = 200, 4 workers, seed 900001, return_variance = TRUE, and saves the
# returned table together with .Random.seed. Run it against dev/perf/lib/ref and dev/perf/lib/dev, then compare the
# two output files with identical() (see layer0_compare.R).
# Optional 4th argument (step 3): TRUE or FALSE, sets options(evolveTrial.use_cpp_engine = ...) before the run. FALSE checks
# the seam (the R body runs); TRUE (the package default) runs the C++ chunk engine where its configuration is supported.
# The reference library has no dispatcher and ignores the option.
a <- commandArgs(TRUE)
if (!length(a) %in% c(3L, 4L)) stop("usage: Rscript layer0_seam_check.R <lib> <args.rds> <out.rds> [use_cpp_engine]")
lib <- normalizePath(a[1]); args_rds <- a[2]; out_rds <- a[3]
use_cpp <- if (length(a) == 4L) as.logical(a[4]) else NA
if (length(a) == 4L && is.na(use_cpp)) stop("use_cpp_engine must be TRUE or FALSE")
if (!is.na(use_cpp)) options(evolveTrial.use_cpp_engine = use_cpp)
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
             rver = as.character(getRversion()), use_cpp_engine = getOption("evolveTrial.use_cpp_engine")), out_rds)
cat(sprintf("layer0 %s lib=%s use_cpp_engine=%s rows=%d wall=%.1fs\n", d$design, basename(lib),
            format(getOption("evolveTrial.use_cpp_engine", "unset")), nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
