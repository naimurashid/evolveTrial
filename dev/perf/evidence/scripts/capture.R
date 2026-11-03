# capture.R -- build the F0 wrapper for one T1 cert design, apply the P6 + dropout-off engine settings the
# retune certify phase applies (p6r_install("p6") + verify_m4_install("off")), and capture the exact
# run_scenarios(base_args, scens, seed) call the wrapper makes. Read-only on all repos.
#   Rscript capture.R <design R01|R23> <out.rds>
a <- commandArgs(TRUE); DES <- a[1]; OUTF <- a[2]
S <- "/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf"
F0 <- "/home/naimrashid/Downloads/atbo-f0"
AL <- "/home/naimrashid/Downloads/atbo-longleaf"
Sys.setenv(SLURM_CPUS_PER_TASK = "1")
suppressPackageStartupMessages({ library(evolveTrial); library(tidyverse); library(BATON) })
stopifnot(startsWith(find.package("evolveTrial"), file.path(S, "Rlib")))
cd <- read.csv(file.path(AL, "results/retune_t1/cert_design.csv"), stringsAsFactors = FALSE)
cr <- cd[cd$design == DES & cd$cert_seed == 900001, ]; stopifnot(nrow(cr) == 1)
TK <- read.delim(file.path(AL, "pilot/retune/staged/t1_tasks.tsv"), stringsAsFactors = FALSE, na.strings = "")
tr <- TK[TK$task == cr$task, ]; stopifnot(nrow(tr) == 1)
SRC <- list(rws = c("R/scenario_parser.R", "R/warmstart_core.R", "R/scenario_configs.R", "R/warmstart_wrappers.R", "scripts/sim_fun_evolveTrial_unified.R"),
            bws = c("R/scenario_parser.R", "R/warmstart_wrappers.R", "R/warmstart_core.R", "scripts/sim_fun_evolveTrial_unified.R"),
            xphil = c("R/warmstart_wrappers.R", "R/warmstart_core.R", "R/scenario_parser.R", "scripts/sim_fun_evolveTrial_unified.R"))
old <- setwd(F0)
suppressMessages(suppressWarnings(for (f in SRC[[tr$protocol]]) source(f, local = FALSE)))
d <- readr::read_csv(tr$csv, show_col_types = FALSE)
spec <- suppressMessages(parse_single_scenario(dplyr::filter(d, scenario_id == tr$scenario_id), verbose = FALSE))
setwd(old)
W0 <- spec$wrapper
theta <- jsonlite::fromJSON(cr$theta_json)
# engine settings (P6 + dropout off), as the retune certify phase installs them
options(evolveTrial.looks_stop_at_final = TRUE)
ns <- asNamespace("evolveTrial"); orig_rs <- get("run_scenarios", envir = ns)
cap <- new.env()
utils::assignInNamespace("run_scenarios", function(base_args, scens, ...) {
  stopifnot(identical(base_args$min_follow_up_at_final, 0), is.null(base_args$random_dropout))
  base_args$min_follow_up_at_final <- 6; base_args$random_dropout <- FALSE
  cap$base_args <- base_args; cap$scens <- scens; cap$dots <- list(...)
  orig_rs(base_args, scens, ...)
}, ns = "evolveTrial")
e <- new.env(parent = environment(W0)); fl <- environment(W0)$fidelity_levels; fl[["high"]] <- 300
assign("fidelity_levels", fl, envir = e); environment(W0) <- e
setwd(F0); t0 <- Sys.time(); res <- W0(as.list(theta), fidelity = "high", seed = 900001); wall <- as.numeric(difftime(Sys.time(), t0, units = "secs")); setwd(old)
utils::assignInNamespace("run_scenarios", orig_rs, ns = "evolveTrial")
saveRDS(list(design = DES, task = cr$task, scenario_id = tr$scenario_id, trial_type = tr$trial_type, protocol = tr$protocol,
             csv = tr$csv, theta = theta, base_args = cap$base_args, scens = cap$scens, dots = cap$dots,
             wrapper_res = res, wrapper_wall = wall), OUTF)
cat(DES, tr$scenario_id, tr$trial_type, "wrapper wall (R=300, workers from SLURM_CPUS_PER_TASK=1):", round(wall, 2), "s\n")
print(res[intersect(names(res), c("power", "type1", "EN_null", "EN_alt", "n_rep"))])
