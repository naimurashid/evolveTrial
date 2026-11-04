# layer3_search.R -- layer 3 of the engine C++ port (plan step 6, dev/plans/2026-11-03-step6-layer3.md): one stored E2c
# search evaluation replayed through E2c's call path with a chosen engine library.
#
#   Rscript layer3_search.R <run_id> <eval_index> <dev|ref> <out.rds>
#
# Derived from dev/perf/layer2/layer2_cert.R. Changes from it carry an "L3-" tag (recon Section 5):
#   L3-a  input row: ledger row and rds values of <run_id> at <eval_index> (seed, counts from the ledger, full-precision
#         theta and expected doubles from the run's own rds). Seed derivation and theta agreement are asserted.
#   L3-b  environment block of E2c (no BO_N_INIT, no BO_SEED, BATON block on BATON rows only), 16 workers.
#   L3-c  wrapper built as E2c built it: parse_warmstart_scenarios on cohortB_multiarm_production.csv.
#   L3-d  engine switches as layer 2 (capture hook innermost, then M4 off, then P6 probe).
#   L3-e  the call: W(theta, fidelity = "low", seed = SEED), R = 3,000.
#   L3-f  checks: 16 workers, 3,000 simulations, two scenarios, cutpoints = K, counts integral.
#   L3-g  the output holds everything layer3_compare.R needs (the comparison itself is in that file).
# Kept from layer 2 unchanged: L2-b path variables (here L3_*), L2-c library identity, L2-d C++ probe, L2-e engine check,
# L2-f FORCE_R negative control (variable L3_FORCE_R), the RNGkind assertion, the no-packages-before-env-block assertion.
# Read-only on every input; writes only <out.rds> and a temporary directory.
`%||%` <- function(x, y) if (is.null(x)) y else x
a <- commandArgs(TRUE)
if (length(a) != 4L) stop("usage: Rscript layer3_search.R <run_id> <eval_index> <dev|ref> <out.rds>")
RUN <- a[1]; EVI <- as.integer(a[2]); MODE <- a[3]; OUTF <- a[4]
stopifnot(grepl("^PB-(baton|ga)-[0-9]+$", RUN), !is.na(EVI), MODE %in% c("dev", "ref"))
if (file.exists(OUTF)) stop("output exists: ", OUTF)
FORCE_R <- identical(Sys.getenv("L3_FORCE_R"), "1")                                                 # L2-f
if (FORCE_R && MODE != "dev") stop("L3_FORCE_R applies to dev only")
ENGINE_SHA <- Sys.getenv("L3_ENGINE_SHA")
if (MODE == "dev" && !grepl("^[0-9a-f]{40}$", ENGINE_SHA)) stop("dev mode needs L3_ENGINE_SHA (40 hex)")
IS_BATON <- grepl("^PB-baton-", RUN)
WORKERS <- 16L                                                                                      # L3-b
BATON_VERSION <- "0.4.0"
SID <- "doublet_triplet_admissible_fut"; CSVREL <- "scenarios/cohortB_multiarm_production.csv"

# ---- paths (L2-b) ------------------------------------------------------------------------------------------------
envv <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
WT <- normalizePath("~/research/evolveTrial-perf", mustWork = FALSE)
AL <- normalizePath(envv("L3_AL", "~/Downloads/atbo-longleaf"))
F0 <- normalizePath(envv("L3_F0", "~/Downloads/atbo-f0"))
LIB <- normalizePath(envv("L3_LIB", file.path(WT, "dev/perf/lib", if (MODE == "dev") "dev_l2" else "ref")))
BL <- normalizePath(envv("L3_BATON_LIB", "~/R_libs_frozen/jasa_r1_v040"))
RES <- normalizePath(envv("L3_RESULTS", "~/pilot_r1/e2c/results_live"))                             # L3-a
P6R <- normalizePath(file.path(AL, "pilot/p6recert/R/p6r_protocol.R"))
M4S <- normalizePath(file.path(AL, "pilot/verify/R/verify_m4_switch.R"))
RR <- file.path(AL, "pilot/retune")
OUT <- tempfile(paste0("l3_", RUN, "_"), tmpdir = envv("TMPDIR", tempdir())); dir.create(OUT, recursive = TRUE); OUT <- normalizePath(OUT)

# ---- inputs: engine switch files, F0 files, git check (as layer 2) -------------------------------------------------
md5f <- function(f) unname(tools::md5sum(f))
if (md5f(P6R) != "52d2b43cf69feccffeca53a9924d4b80" || md5f(M4S) != "405f7d60e1a4be532821e3ba603227bc")
  stop("engine switch files differ from the certify manifests")
chk_list <- function(listing, base) { x <- read.table(listing, stringsAsFactors = FALSE, col.names = c("md5", "path"))
  bad <- x$path[md5f(file.path(base, x$path)) != x$md5 | !file.exists(file.path(base, x$path))]
  if (length(bad)) stop("md5 mismatch: ", paste(bad, collapse = ", ")) }
chk_list(file.path(RR, "staged", "f0_files.md5"), F0)
F0_SHA <- "1863cecb928eca56d4e52fe202542630568a5f35"
git_sha <- tryCatch(suppressWarnings(system2("git", c("-C", F0, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
if (length(git_sha) && !identical(git_sha, F0_SHA)) stop("F0 worktree at ", git_sha, ", not 1863cecb")

# ---- L3-a: the stored row ----------------------------------------------------------------------------------------
LEDF <- file.path(RES, RUN, "sim_ledger.csv")
RDSF <- file.path(RES, RUN, if (IS_BATON) "f0_out/doublet_triplet_admissible_fut_warmstart_result.rds" else "run_state.rds")
sha256 <- function(f) sub(" .*", "", system2("sha256sum", f, stdout = TRUE))
led <- utils::read.csv(LEDF, stringsAsFactors = FALSE)
lr <- led[led$purpose == "search" & led$eval_index == EVI, ]
if (nrow(lr) != 1L) stop("ledger: not exactly one search row for eval ", EVI)
if (lr$R != 3000L || lr$worker_count != 16L || !identical(lr$error_flag, FALSE) || !identical(lr$invalid_design, FALSE) ||
    !startsWith(lr$git_sha_F0, "1863cecb")) stop("ledger row fails the L3-a assertions")
SSEED <- lr$search_seed
SEED <- as.integer(lr$seed)
if (SEED != if (IS_BATON) as.integer(1e9 + SSEED * 1e5 + EVI) else as.integer(SSEED * 1e5 + EVI)) stop("seed derivation fails")
PN <- c("total_n", "eff", "fut", "ev", "beat", "margin", "hr_margin", "prior_strength", "n_intervals")

# ---- L3-b: environment block, set before any package load ---------------------------------------------------------
stray <- grep("^BO_", names(Sys.getenv()), value = TRUE)
if (length(stray)) stop("BO_* variables already set (fail closed): ", paste(stray, collapse = ", "))
if (nzchar(Sys.getenv("MC_CORES")) || nzchar(Sys.getenv("BO_SCALING_PROFILE"))) stop("MC_CORES or BO_SCALING_PROFILE set")
SLW <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (nzchar(SLW) && !identical(as.integer(SLW), WORKERS)) stop("SLURM_CPUS_PER_TASK ", SLW, " != workers ", WORKERS)
ENV <- c(BATON_DISABLE_KSNAP = "1", SLURM_CPUS_PER_TASK = as.character(WORKERS))
if (IS_BATON) ENV <- c(ENV, BO_LOW_REPS = "3000", BO_MED_REPS = "5000", BO_HIGH_REPS = "10000", BO_USE_BUFFER = "TRUE",
                       BO_CONSTRAINT_BUFFER = "0.02", BO_VALIDATION_TIER_SIZES = "5,5", BO_VALIDATION_MIN_STOP = "10",
                       BO_BUDGET_S1 = "700", BO_BUDGET_S2 = "300", BO_BUDGET_S3 = "200")
do.call(Sys.setenv, as.list(ENV))
if (length(setdiff(loadedNamespaces(), c("compiler", "graphics", "utils", "grDevices", "stats", "datasets", "methods", "base", "tools"))))
  stop("packages loaded before the environment block: ", paste(loadedNamespaces(), collapse = " "))

# ---- libraries (L2-c) ---------------------------------------------------------------------------------------------
.libPaths(unique(c(LIB, .libPaths())))
loadNamespace("evolveTrial", lib.loc = LIB)
if (dirname(normalizePath(getNamespaceInfo("evolveTrial", "path"))) != LIB) stop("evolveTrial not from ", LIB)
if (IS_BATON) {
  suppressPackageStartupMessages({ library(tidyverse); library(BATON, lib.loc = BL) })
  if (!identical(as.character(utils::packageVersion("BATON")), BATON_VERSION)) stop("BATON is ", utils::packageVersion("BATON"), ", not ", BATON_VERSION)
  if (dirname(normalizePath(getNamespaceInfo("BATON", "path"))) != normalizePath(BL)) stop("BATON not from ", BL)
}
if (dirname(normalizePath(getNamespaceInfo("evolveTrial", "path"))) != LIB) stop("evolveTrial reloaded from another library")
ED <- utils::packageDescription("evolveTrial", lib.loc = LIB)
if (MODE == "ref") {
  if (!identical(ED$VerifyM4PatchMD5, "83de2c8b2de1902630eddf00c334b508") || !identical(ED$P6RecertEvolveTrialSHA, "7565a17c625573d2a697d0d4a85efb2fa6c58ad9"))
    stop("library ", LIB, " is not the M4 build (P6 7565a17 + m4_nodropout.patch)")
} else {
  if (!identical(ED$P6RecertEvolveTrialSHA, "7565a17c625573d2a697d0d4a85efb2fa6c58ad9") || !identical(ED$EngineCppSHA, ENGINE_SHA) ||
      !is.null(ED$VerifyM4PatchMD5))
    stop("library ", LIB, " is not the build_dev_tree.sh build of ", ENGINE_SHA, " (EngineCppSHA ", ED$EngineCppSHA %||% "NA", ")")
  if (!exists(".simulate_chunk_cpp", envir = asNamespace("evolveTrial"), inherits = FALSE)) stop("no C++ chunk engine in ", LIB)
  if (FORCE_R) options(evolveTrial.use_cpp_engine = FALSE)                                          # L2-f
}
if (dirname(getNamespaceInfo("evolveTrial", "path")) != normalizePath(LIB)) stop("plan assertion: evolveTrial path")
if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))) stop("RNGkind is ", paste(RNGkind(), collapse = "/"))

# ---- L3-a (second half): the rds is read after the libraries are loaded, so that reading it loads no package early
rd <- readRDS(RDSF)
if (IS_BATON) {
  h <- rd$history_combined; k <- which(h$fidelity == "low" & h$eval_id == EVI)
  if (length(k) != 1L) stop("rds: not exactly one low-fidelity row for eval ", EVI)
  th <- as.list(h$theta[[k]])[PN]
  EXP <- list(objective = h$objective[k], EN_null = h$EN_null[k], EN_alt = h$EN_alt[k], power = h$power[k], type1 = h$type1[k])
} else {
  h <- rd$history; k <- which(h$eval_index == EVI)
  if (length(k) != 1L) stop("rds: not exactly one history row for eval ", EVI)
  th <- as.list(h[k, PN])
  EXP <- list(objective = h$objective[k], power = h$power[k], type1 = h$type1[k])
}
lth <- unlist(lr[paste0("theta_", PN)]); names(lth) <- PN
if (any(abs(unlist(th) - lth) > 1e-13 * pmax(1, abs(lth)))) stop("rds theta differs from ledger theta beyond 1e-13")
LED_COUNTS <- c(power_count = round(lr$power * 3000), type1_count = round(lr$type1 * 3000))
if (abs(lr$power * 3000 - LED_COUNTS[[1]]) > 1e-6 || abs(lr$type1 * 3000 - LED_COUNTS[[2]]) > 1e-6) stop("ledger counts not integral")
IN_SHA <- c(ledger = sha256(LEDF), rds = sha256(RDSF))

# ---- L3-c: wrapper built as E2c built it (recon Section 2b) --------------------------------------------------------
old <- setwd(F0)
suppressMessages(suppressWarnings(source("R/scenario_parser.R", local = FALSE)))
if (IS_BATON) suppressMessages(suppressWarnings(source("R/warmstart_core.R", local = FALSE)))
sc <- suppressMessages(parse_warmstart_scenarios(CSVREL, verbose = FALSE))
spec <- sc[[SID]]
setwd(old)
if (is.null(spec) || !identical(spec$scenario_id, SID)) stop("scenario ", SID, " not parsed")
W0 <- spec$wrapper; we <- environment(W0)
if (!identical(as.integer(we$cached_num_workers), WORKERS)) stop("wrapper cached ", we$cached_num_workers, " workers, expected ", WORKERS)
if (!identical(attr(W0, "trial_type"), "multi_arm")) stop("wrapper trial_type ", attr(W0, "trial_type"), " != multi_arm")
if (!isTRUE(all.equal(c(we$w_total_n, we$w_null, we$w_alt), c(0.5, 1, 0), tolerance = 1e-12)) || !isTRUE(all.equal(we$w_EN, 1 - we$w_total_n)))
  stop("wrapper weights differ from 0.5/1/0")
if (!identical(spec$objective, "weighted_EN_N")) stop("objective ", spec$objective, " != weighted_EN_N")
FL <- c(low = 3000, med = 5000, high = 10000)
if (!isTRUE(all.equal(unlist(we$fidelity_levels)[c("low", "med", "high")], FL, check.attributes = FALSE))) stop("F0 wrapper fidelity is not 3000/5000/10000")
CSV_MD5 <- md5f(file.path(F0, CSVREL))

# ---- L3-d: capture hook innermost, then the two engine switches -----------------------------------------------------
ns <- asNamespace("evolveTrial")
.cap <- new.env(parent = emptyenv()); .cap$calls <- list()
.cap$orig_rs <- get("run_scenarios", envir = ns)
utils::assignInNamespace("run_scenarios", function(base_args, scens, ...) {
  .cap$calls[[length(.cap$calls) + 1L]] <- list(base_args = base_args, scens = scens, dots = list(...),
                                               looks_stop_at_final_option = getOption("evolveTrial.looks_stop_at_final"),
                                               seed_before = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
  .cap$orig_rs(base_args, scens, ...)
}, ns = "evolveTrial")
HOOK <- get("run_scenarios", envir = ns)
source(M4S); source(P6R)
M4LOG <- file.path(OUT, "m4_switch_calls.tsv"); M4PROBE <- file.path(OUT, "m4_probe"); P6EV <- file.path(OUT, "p6_events")
verify_m4_install("off", M4LOG, M4PROBE)
p6r_install("p6", "probe", ev_dir = file.path(P6EV, "_init"), expected_lib = LIB)
if (!identical(.vm4$orig_rs, HOOK)) stop("M4 switch does not wrap the capture hook")
CPPD <- file.path(OUT, "cpp_probe"); dir.create(CPPD)                                              # L2-d
if (MODE == "dev") {
  .l2 <- new.env(parent = emptyenv()); .l2$orig_cpp <- get(".simulate_chunk_cpp", envir = ns)
  if (!identical(names(formals(.l2$orig_cpp)), c("sim_indices", "cfg"))) stop(".simulate_chunk_cpp formals changed")
  utils::assignInNamespace(".simulate_chunk_cpp", function(sim_indices, cfg) {
    cat(sprintf("C\t%d\t%d\t%s\t%s\t%s\t%s\n", Sys.getpid(), length(sim_indices), format(cfg$random_dropout),
                format(cfg$min_follow_up_at_final), format(cfg$looks_stop_at_final), format(cfg$mode)),
        file = file.path(CPPD, sprintf("cpp_%d.tsv", Sys.getpid())), append = TRUE)
    .l2$orig_cpp(sim_indices, cfg)
  }, ns = "evolveTrial")
}
cpp_lines <- function() unlist(lapply(sort(list.files(CPPD, full.names = TRUE)), function(f) readLines(f, warn = FALSE)))
rd_lines <- function(f) if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
probe_lines <- function() unlist(lapply(sort(list.files(M4PROBE, full.names = TRUE)), rd_lines))
EVF <- file.path(OUT, "evals.csv"); .k <- 0L
engine_fail <- function(msg) { writeLines(c(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), msg), file.path(OUT, "ENGINE-FAIL.txt"))
  cat("[layer3] ENGINE FAIL:", msg, "\n"); quit(save = "no", status = 3) }
W <- function(theta, fidelity = "low", seed = NULL, ...) {
  .k <<- .k + 1L; t0 <- Sys.time()
  n4 <- length(rd_lines(M4LOG)); np <- length(probe_lines()); n6 <- p6r_ncalls(); ncp <- length(cpp_lines()); ncap <- length(.cap$calls)
  p6r_set_ev_dir(file.path(P6EV, sprintf("c%05d", .k)))
  res <- tryCatch(W0(theta, fidelity = fidelity, seed = seed, ...), error = function(err) err)
  ok <- !inherits(res, "error")
  m4n <- rd_lines(M4LOG); m4n <- if (length(m4n) > n4) m4n[(n4 + 1):length(m4n)] else character(0)
  pr <- probe_lines(); pr <- if (length(pr) > np) pr[(np + 1):length(pr)] else character(0)
  nc <- p6r_ncalls() - n6; ck <- p6r_check_calls(n6, nc)
  ev <- tryCatch(p6r_read_events(.p6r$ev_dir)$W, error = function(e) NULL)
  m4sw <- unique(vapply(strsplit(m4n, "\t"), `[`, "", 3)); m4pass <- unique(vapply(strsplit(m4n, "\t"), `[`, "", 4))
  prv <- unique(vapply(strsplit(pr, "\t"), `[`, "", 3))
  p6F <- if (is.null(ev)) character(0) else unique(ev$frame_F); p6o <- if (is.null(ev)) character(0) else unique(ev$frame_opt)
  p6b <- if (is.null(ev)) character(0) else unique(ev$build)
  echo <- sprintf("p6[calls=%d F=%s opt=%s workers=%d frameF=%s frameOpt=%s build=%s] dropout[calls=%d switch=%s passed=%s workers=%d seen=%s]",
                  nc, paste(ck$F_values, collapse = "/"), paste(ck$option_values, collapse = "/"), if (is.null(ev)) 0L else length(unique(ev$pid)),
                  paste(p6F, collapse = "/"), paste(p6o, collapse = "/"), substr(paste(p6b, collapse = "/"), 1, 7),
                  length(m4n), paste(m4sw, collapse = "/"), paste(m4pass, collapse = "/"), length(pr), paste(prv, collapse = "/"))
  cp <- cpp_lines(); cp <- if (length(cp) > ncp) cp[(ncp + 1):length(cp)] else character(0)
  cps <- strsplit(cp, "\t", fixed = TRUE); cpf <- function(i) unique(vapply(cps, `[`, "", i))
  cp_trials <- sum(as.integer(vapply(cps, `[`, "", 3))); cp_pids <- length(cpf(2))
  Rf <- tryCatch(unname(environment(W0)$fidelity_levels[[as.character(fidelity)]]), error = function(e) NA_real_)
  nscen <- sum(vapply(.cap$calls[seq_along(.cap$calls) > ncap], function(z) length(z$scens), 0L))
  echo <- paste0(echo, sprintf(" cpp[lines=%d pids=%d trials=%d expect=%s rd=%s F=%s lsf=%s mode=%s] force_r=%s", length(cp), cp_pids,
                               cp_trials, format(nscen * Rf), paste(cpf(4), collapse = "/"), paste(cpf(5), collapse = "/"),
                               paste(cpf(6), collapse = "/"), paste(cpf(7), collapse = "/"), FORCE_R))
  eng_ok <- if (MODE == "ref") {
    !ok || (nc >= 1L && ck$ok && length(m4n) == nc && identical(m4sw, "off") && identical(m4pass, "FALSE") &&
            length(pr) >= 1L && identical(prv, "FALSE") && identical(p6F, "6") && identical(p6o, "TRUE") &&
            identical(p6b, P6R_CONST$et_sha))
  } else {
    !ok || (nc >= 1L && ck$ok && length(m4n) == nc && identical(m4sw, "off") && identical(m4pass, "FALSE") &&
            length(pr) == 0L && NROW(ev) == 0L && nscen >= 1L && is.finite(Rf) && cp_trials == nscen * Rf &&
            identical(cpf(4), "FALSE") && identical(cpf(5), "6") && identical(cpf(6), "TRUE"))
  }
  row <- data.frame(call = .k, fidelity = as.character(fidelity), R = Rf, seed = seed %||% NA, status = if (ok) "ok" else paste("error:", conditionMessage(res)),
                    engine_ok = eng_ok, engine_echo = echo, wall_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  utils::write.table(row, EVF, sep = ",", row.names = FALSE, col.names = !file.exists(EVF), append = file.exists(EVF), qmethod = "double")
  if (!eng_ok) engine_fail(sprintf("call %d: %s", .k, echo))
  if (!ok) stop(res)
  res
}

# ---- L3-e: the call, in the F0 working directory --------------------------------------------------------------------
setwd(F0)
t0 <- Sys.time()
m <- W(th, fidelity = "low", seed = SEED)
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
setwd(old)
RW <- 3000L
cnt <- function(p) { x <- p * RW; if (is.finite(x) && abs(x - round(x)) < 1e-6) round(x) else NA_real_ }
cc <- c(power_count = cnt(m[["power"]]), type1_count = cnt(m[["type1"]]))

# ---- L3-f: checks ----------------------------------------------------------------------------------------------------
if (length(.cap$calls) != 1L) stop("expected exactly one run_scenarios call, saw ", length(.cap$calls))
cap1 <- .cap$calls[[1]]; ba <- cap1$base_args
K <- as.integer(th$n_intervals)
chk <- c(cutpoints_eq_n_intervals = length(ba$interval_cutpoints_sim) - 1L == K,
         random_dropout_FALSE = identical(ba$random_dropout, FALSE),
         min_follow_up_at_final_6 = identical(ba$min_follow_up_at_final, 6),
         looks_stop_at_final_option = isTRUE(getOption("evolveTrial.looks_stop_at_final")) && isTRUE(cap1$looks_stop_at_final_option),
         num_workers_16 = identical(as.numeric(ba$num_workers), 16),
         num_simulations_R = identical(as.numeric(ba$num_simulations), as.numeric(RW)),
         two_scenarios = length(cap1$scens) == 2L,
         n_rep_R = identical(as.numeric(attr(m, "n_rep") %||% NA), as.numeric(RW)),
         counts_integral = !anyNA(cc))
evrow <- utils::read.csv(EVF, stringsAsFactors = FALSE)
OBJ <- m[["weighted_EN_N"]]
# L3-g
saveRDS(list(mode = MODE, force_r = FORCE_R, lib = LIB, engine_sha = ED$EngineCppSHA %||% NA_character_,
             p6_build = ED$P6RecertEvolveTrialSHA %||% NA_character_, m4_patch = ED$VerifyM4PatchMD5 %||% NA_character_,
             engine_ok = isTRUE(evrow$engine_ok), engine_echo = evrow$engine_echo, cpp_probe = cpp_lines(),
             run_id = RUN, eval_index = EVI, is_baton = IS_BATON, seed = SEED, theta = th, K = K,
             counts = cc, power = m[["power"]], type1 = m[["type1"]], objective = OBJ,
             EN_null = m[["EN_null"]], EN_alt = m[["EN_alt"]],
             expected = EXP, ledger_counts = LED_COUNTS, input_sha256 = IN_SHA, csv_md5 = CSV_MD5,
             r_version = R.version.string, host = Sys.info()[["nodename"]], R = RW, workers = WORKERS,
             checks = chk, evals_row = evrow, wall_s = wall, env = ENV, sessionInfo = sessionInfo()), OUTF)
cat(sprintf("[layer3] %s eval %d %s seed %d K %d power_count %d type1_count %d objective %.17g E0N %.17g E1N %.17g wall %.1f s\n",
            RUN, EVI, MODE, SEED, K, as.integer(cc[1]), as.integer(cc[2]), OBJ, m[["EN_null"]], m[["EN_alt"]], wall))
cat("[layer3] engine echo:", evrow$engine_echo, "engine_ok", evrow$engine_ok, "\n")
cat("[layer3] checks:", paste(names(chk), chk, sep = "=", collapse = " "), "\n")
if (!all(chk)) { cat("[layer3] CHECK FAILED\n"); quit(save = "no", status = 4) }
