#' Summarise simulation output by scenario and arm
#'
#' Aggregates the per-arm results returned by `run_scenarios()` and pivots them
#' into a wide table (one row per scenario) for downstream reporting.
#'
#' @param results_df Data frame/data.table produced by `run_scenarios()`
#'   containing at least `scenario` and `Arm_Name`.
#'
#' @return A data.frame with one row per scenario and columns for each arm's key
#'   operating characteristics.
#' @export
pretty_scenario_matrix <- function(results_df) {
  # Validate
  if (!"scenario" %in% names(results_df)) stop("results_df must include 'scenario'")
  if (!"Arm_Name" %in% names(results_df)) stop("results_df must include 'Arm_Name'")
  
  # Aggregate summary by scenario + arm
  tbl <- data.table::as.data.table(results_df)[, .(
    True_Median         = mean(True_Median, na.rm = TRUE),
    Max_Planned_N       = max(Exp_N, na.rm = TRUE),
    Exp_N               = mean(Exp_N, na.rm = TRUE),
    Exp_Events          = mean(Exp_Events, na.rm = TRUE),
    Pr_Reach_Max_N      = mean(Pr_Reach_Max_N, na.rm = TRUE),
    Type_I_Error_or_Power = mean(Type_I_Error_or_Power, na.rm = TRUE),
    PET_Efficacy        = mean(PET_Efficacy, na.rm = TRUE),
    PET_Futility        = mean(PET_Futility, na.rm = TRUE),
    Pr_Final_Efficacy   = mean(Pr_Final_Efficacy, na.rm = TRUE),
    Pr_Final_Futility   = mean(Pr_Final_Futility, na.rm = TRUE)
  ), by = .(scenario, Arm_Name)]
  
  # Pivot wide: one row per scenario, arms side-by-side
  wide <- data.table::dcast(
    tbl,
    scenario ~ Arm_Name,
    value.var = c("True_Median", "Exp_N", "Exp_Events", "Pr_Reach_Max_N",
                  "Type_I_Error_or_Power", "PET_Efficacy", "PET_Futility",
                  "Pr_Final_Efficacy", "Pr_Final_Futility"),
    fun.aggregate = mean
  )
  
  # Optional rounding for cleaner display
  num_cols <- names(wide)[sapply(wide, is.numeric)]
  wide[, (num_cols) := lapply(.SD, function(x) round(x, 2)), .SDcols = num_cols]
  
  # Return formatted data.frame (for printing)
  as.data.frame(wide)
}

#' Export a scenario summary table to Excel
#'
#' @param pretty_tbl Data frame produced by `pretty_scenario_matrix()`.
#' @param file_path Output path for the Excel workbook.
#'
#' @return Invisibly returns `file_path`.  Writes an `.xlsx` file to disk.
#' @export
export_scenario_table_to_excel <- function(pretty_tbl, file_path = "scenario_summary.xlsx") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Install with install.packages('openxlsx')")
  }
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Scenario Summary")
  openxlsx::writeDataTable(wb, sheet = 1, x = pretty_tbl, tableStyle = "TableStyleMedium9")
  openxlsx::setColWidths(wb, sheet = 1, cols = 1:ncol(pretty_tbl), widths = "auto")
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message("[OK] Exported formatted table to: ", normalizePath(file_path))
}

#' Render the scenario summary table to a PNG image
#'
#' @param results_df Data frame returned by `run_scenarios()`.
#' @param file_path Output path for the PNG image.
#' @param title Main title for the table.
#' @param subtitle Optional subtitle.
#' @param highlight_arm Arm name whose columns should be highlighted.
#' @param snapshot_engine Rendering backend; one of `"auto"`, `"webshot2"`,
#'   or `"webshot"`.
#' @param vwidth Viewport width passed to the renderer.
#' @param vheight Viewport height passed to the renderer.
#' @param zoom Zoom factor passed to the renderer.
#'
#' @return Invisibly returns `file_path`.  Writes a PNG (and temporary HTML if
#'   `snapshot_engine = "webshot"`).
#' @export
export_scenario_table_to_png <- function(results_df,
                                         file_path = "scenario_summary.png",
                                         title = "Bayesian Adaptive Design Summary",
                                         subtitle = NULL,
                                         highlight_arm = "Triplet",
                                         snapshot_engine = c("auto", "webshot2", "webshot"),
                                         vwidth = 1400,
                                         vheight = 600,
                                         zoom = 1) {
  if (!requireNamespace("gt", quietly = TRUE)) stop("Please install.packages('gt')")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  
  snapshot_engine <- match.arg(snapshot_engine)
  
  # engine auto-detect
  if (snapshot_engine == "auto") {
    if (requireNamespace("webshot2", quietly = TRUE)) {
      snapshot_engine <- "webshot2"
    } else if (requireNamespace("webshot", quietly = TRUE)) {
      snapshot_engine <- "webshot"
    } else {
      stop("Install either 'webshot2' (Chrome/Chromium) or 'webshot' (PhantomJS).")
    }
  }
  
  # build table data (expects you already defined pretty_scenario_matrix)
  tbl <- pretty_scenario_matrix(results_df)

  # minimal formatting
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required. Install with install.packages('gt')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with install.packages('dplyr')")
  }


  
  gt_tbl <- gt::gt(tbl) %>%
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) %>%
    gt::fmt_number(columns = where(is.numeric), decimals = 2) %>%
    gt::opt_align_table_header("left") %>%
    gt::tab_options(
      table.font.size = 12,
      data_row.padding = gt::px(4),
      heading.background.color = "#f0f0f0",
      column_labels.background.color = "#f7f7f7",
      table.border.top.width = gt::px(1),
      table.border.bottom.width = gt::px(1),
      table.border.bottom.color = "gray"
    )
  
  if (!is.null(highlight_arm)) {
    highlight_cols <- grep(highlight_arm, names(tbl), value = TRUE)
    if (length(highlight_cols) > 0) {
      gt_tbl <- gt_tbl %>%
        gt::tab_style(
          style = gt::cell_fill(color = "#e8f4f8"),
          locations = gt::cells_body(columns = dplyr::all_of(highlight_cols))
        )
    }
  }
  
  # branch by engine
  if (snapshot_engine == "webshot2") {
    # Directly save PNG (gt uses webshot2/chromote internally)
    gt::gtsave(
      data = gt_tbl,
      filename = file_path,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
  } else {
    # Save HTML then rasterize with webshot (PhantomJS)
    if (!requireNamespace("webshot", quietly = TRUE)) {
      stop("snapshot_engine='webshot' requires the 'webshot' package.")
    }
    html_tmp <- sub("\\.png$", ".html", file_path, ignore.case = TRUE)
    gt::gtsave(
      data = gt_tbl,
      filename = html_tmp,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
    webshot::webshot(
      url = html_tmp,
      file = file_path,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
  }
  
  message("[OK] Image saved: ", normalizePath(file_path))
}

#' Calibrate interim and final thresholds for single-arm designs
#'
#' Sweeps candidate interim and final posterior probability thresholds and
#' returns the best combination achieving the desired type I error under the
#' null scenario(s).
#'
#' @param base_args Baseline argument list passed to `run_scenarios()`.
#' @param scens_null Scenario list representing null hypotheses.
#' @param thr_grid_interim Numeric vector of interim success thresholds to try.
#' @param thr_grid_final Numeric vector of final success thresholds to try.
#' @param sims Number of simulations per candidate setting.
#'
#' @return A list containing the chosen thresholds and corresponding estimated
#'   type I error.
#' @export
calibrate_alpha <- function(base_args, scens_null, thr_grid_interim = c(0.9, 0.95, 0.975),
                            thr_grid_final = c(0.95, 0.975, 0.99), sims = 300) {
  best <- NULL
  base_args$num_simulations <- sims
  for (ti in thr_grid_interim) for (tf in thr_grid_final) {
    args <- base_args
    args$efficacy_threshold_hc_prob <- ti
    args$final_success_posterior_prob_threshold <- tf
    res <- run_scenarios(args, scens_null, parallel = TRUE, seed = 123)
    # Type I = early+final efficacy for the experimental arm(s) under null scenario
    typeI <- res[res$Arm_Name != args$reference_arm_name & res$scenario == 1, "Type_I_Error_or_Power"]
    alpha_hat <- mean(typeI)
    cand <- list(ti = ti, tf = tf, alpha = alpha_hat)
    if (is.null(best) || (alpha_hat < best$alpha && alpha_hat <= 0.10)) best <- cand
    cat(sprintf("Interim=%.3f Final=%.3f => alpha=%.3f\n", ti, tf, alpha_hat))
  }
  best
}

# === Grid search for Type I / Power vs thresholds & margin ===
# - Sweeps interim success threshold (vs HC), final success threshold, and Triplet success margin
# - Evaluates Type I under all-null (Doublet=6, Triplet=6) and Power under alt (Doublet=6, Triplet=9)
# - Ranks designs that satisfy alpha <= target_alpha by highest power, then lowest Exp_N under alt

#' Calibrate historical-control thresholds over a grid
#'
#' Evaluates a grid of interim/final thresholds and superiority margins under
#' null and alternative scenarios, returning full operating characteristics for
#' each combination.
#'
#' @param base_args Baseline argument list passed to `run_scenarios()`.
#' @param null_med Control-arm median under the null hypothesis.
#' @param alt_med Experimental median under the alternative.
#' @param margins_abs Numeric vector of absolute superiority margins (months).
#' @param interim_thr_grid Interim success probabilities to evaluate.
#' @param final_thr_grid Final posterior success probabilities to evaluate.
#' @param sims Number of simulations per grid point.
#' @param target_alpha Target type I error used when ranking feasible designs.
#' @param seed RNG seed.
#' @param parallel Logical; run scenarios in parallel.
#'
#' @return A list with components `all`, `feasible`, and `top` summarising the
#'   grid results.
#' @export
grid_calibrate <- function(base_args,
                           null_med = 6,
                           alt_med  = 9,
                           margins_abs = c(1, 2, 3),          # success margin for Triplet in months
                           interim_thr_grid = c(0.90, 0.95),  # Pr(success|data) at interim vs HC
                           final_thr_grid   = c(0.95, 0.975, 0.99), # Pr(success|final posterior)
                           sims = 400,                         # increase for stability
                           target_alpha = 0.10,
                           seed = 123,
                           parallel = TRUE) {
  
  # Build two scenarios: (1) all-null, (2) Triplet=alt
  scens <- scenarios_from_grid(list(
    weibull_median_true_arms = list(
      c(Doublet = null_med, Triplet = null_med),
      c(Doublet = null_med, Triplet = alt_med)
    )
  ))
  
  combos <- CJ(margin = margins_abs,
               interim_thr = interim_thr_grid,
               final_thr   = final_thr_grid)
  
  results <- vector("list", nrow(combos))
  
  for (i in seq_len(nrow(combos))) {
    m  <- combos$margin[i]
    ti <- combos$interim_thr[i]
    tf <- combos$final_thr[i]
    
    args_i <- base_args
    
    # --- thresholds being calibrated ---
    args_i$efficacy_threshold_hc_prob <- ti
    args_i$final_success_posterior_prob_threshold <- tf
    
    # --- require superiority margin for Triplet only (Doublet stays at null) ---
    args_i$median_pfs_success_threshold_arms <- c(Doublet = null_med, Triplet = null_med + m)
    
    # --- keep predictive disabled during calibration (safer for Type I) ---
    args_i$predictive_fast <- FALSE
    
    # --- set sims & seed ---
    args_i$num_simulations <- sims
    
    # run
    res <- run_scenarios(args_i, scens, parallel = parallel, seed = seed)
    
    # Pull alpha (scenario 1, Triplet arm) & power (scenario 2, Triplet arm)
    r_null <- res[res$scenario == 1 & res$Arm_Name == "Triplet", ]
    r_alt  <- res[res$scenario == 2 & res$Arm_Name == "Triplet", ]
    
    alpha_hat <- mean(r_null$Type_I_Error_or_Power)
    power_hat <- mean(r_alt$Type_I_Error_or_Power)
    
    out <- data.table(
      margin_abs = m,
      interim_thr = ti,
      final_thr = tf,
      alpha = alpha_hat,
      power = power_hat,
      ExpN_null = mean(r_null$Exp_N),
      ExpN_alt  = mean(r_alt$Exp_N),
      ExpEvents_null = mean(r_null$Exp_Events),
      ExpEvents_alt  = mean(r_alt$Exp_Events),
      PET_Eff_null = mean(r_null$PET_Efficacy),
      PET_Eff_alt  = mean(r_alt$PET_Efficacy),
      PET_Fut_null = mean(r_null$PET_Futility),
      PET_Fut_alt  = mean(r_alt$PET_Futility)
    )
    
    results[[i]] <- out
  }
  
  grid <- rbindlist(results)
  
  # Designs meeting the alpha target
  ok <- grid[alpha <= target_alpha]
  setorder(ok, -power, ExpN_alt, margin_abs, interim_thr, final_thr)
  
  list(all = grid[order(margin_abs, interim_thr, final_thr)],
       feasible = ok,
       top = head(ok, 10))
}


#' Evaluate PH-based grid of designs
#'
#' Helper that sweeps probability thresholds, gate settings, and HR margins,
#' returning operating characteristics plus expected sample sizes per arm.
#' @param base_args List of defaults passed to `run_scenarios()`.
#' @param grid data.table/data.frame with at least the columns
#'   `label`, `thr_eff`, `thr_fut`, `margin`, `min_ev`, `min_pt`,
#'   and `hr_margin`.
#' @param scens Scenario list (from `scenarios_from_grid()`).
#' @param sims Simulations per grid row.
#' @param seed RNG seed.
#' @param parallel Whether to use `parallel::mclapply`.
#' @return data.table summarising Type I error, power, PETs, expected N per arm,
#'   control-arm expectations, and total expected N under null/alt.
#' Evaluate proportional-hazards vs-reference designs over a grid
#'
#' Runs the supplied grid of thresholds/margins against a scenario list,
#' collecting operating characteristics for each design.
#'
#' @param base_args Baseline argument list passed to `run_scenarios()`.
#' @param grid data.table/data.frame describing the design grid.  Must include
#'   columns `label`, `thr_eff`, `thr_fut`, `margin`, `min_ev`, `min_pt`, and
#'   optionally `hr_margin`.
#' @param scens Scenario list (e.g., from `scenarios_from_grid()`).
#' @param sims Number of simulations per grid row.
#' @param seed RNG seed.
#' @param parallel Logical; use parallel execution.
#'
#' @return data.table summarising alpha, power, PETs, and expected sample sizes
#'   for each labelled design.
#' @export
evaluate_ph_grid <- function(base_args, grid, scens, sims = 2000,
                             seed = 4242, parallel = TRUE) {
  run_one <- function(i) {
    row <- grid[i]
    args_i <- base_args
    args_i$num_simulations                <- sims
    args_i$efficacy_threshold_vs_ref_prob <- row$thr_eff
    args_i$futility_threshold_vs_ref_prob <- row$thr_fut
    args_i$compare_arms_futility_margin   <- row$margin
    args_i$compare_arms_hr_margin         <- row$hr_margin
    args_i$min_events_per_arm             <- row$min_ev
    args_i$min_person_time_frac_per_arm   <- row$min_pt
    args_i$diagnostics                    <- FALSE

    res <- run_scenarios(args_i, scens, parallel = FALSE, seed = seed)
    res_dt <- data.table::as.data.table(res)

    ref_arm <- args_i$reference_arm_name
    if (is.null(ref_arm) || !ref_arm %in% args_i$arm_names) {
      ref_arm <- args_i$arm_names[1]
    }

    summary <- res_dt[Arm_Name == "Triplet",
                      .(alpha = mean(Type_I_Error_or_Power[scenario == 1]),
                        power = mean(Type_I_Error_or_Power[scenario == 2]),
                        PET_Eff_null = mean(PET_Efficacy[scenario == 1]),
                        PET_Eff_alt  = mean(PET_Efficacy[scenario == 2]),
                        PET_Fut_null = mean(PET_Futility[scenario == 1]),
                        PET_Fut_alt  = mean(PET_Futility[scenario == 2]),
                        ExpN_null    = mean(Exp_N[scenario == 1]),
                        ExpN_alt     = mean(Exp_N[scenario == 2]),
                        ExpEvents_null = mean(Exp_Events[scenario == 1]),
                        ExpEvents_alt  = mean(Exp_Events[scenario == 2]))]

    ctrl_null <- res_dt[scenario == 1 & Arm_Name == ref_arm, mean(Exp_N)]
    ctrl_alt  <- res_dt[scenario == 2 & Arm_Name == ref_arm, mean(Exp_N)]
    total_null <- res_dt[scenario == 1, sum(Exp_N)]
    total_alt  <- res_dt[scenario == 2, sum(Exp_N)]

    summary[, `:=`(
      ExpN_ctrl_null = ctrl_null %||% NA_real_,
      ExpN_ctrl_alt  = ctrl_alt %||% NA_real_,
      ExpN_total_null = total_null %||% NA_real_,
      ExpN_total_alt  = total_alt %||% NA_real_
    )]

    cbind(row, summary)
  }

  idx <- seq_len(nrow(grid))
  out <- if (parallel) {
    parallel::mclapply(idx, run_one, mc.cores = max(1L, parallel::detectCores() - 1L))
  } else {
    lapply(idx, run_one)
  }
  data.table::rbindlist(out, use.names = TRUE)
}

# =========================
# Power vs Type I plots
# =========================
# Requires: data.table, ggplot2
# Optional: ggrepel (for nicer labels), scales (percent axes)


#' Plot power versus type I error for calibration grids
#'
#' Visualises the output of `grid_calibrate()` by plotting power against type I
#' error, optionally highlighting Pareto frontiers and feasible designs.
#'
#' @param cal List returned by `grid_calibrate()`.
#' @param target_alpha Type I error threshold to display.
#' @param label_top_n Number of feasible designs to annotate.
#'
#' @return A ggplot object.
#' @export
plot_calibration <- function(cal,
                             target_alpha = 0.10,
                             label_top_n = 3) {
  stopifnot(is.list(cal), !is.null(cal$all))
  df <- as.data.table(cal$all)
  if (nrow(df) == 0L) stop("cal$all is empty")
  
  # Facet label helpers
  df[, margin_lab := sprintf("Delta = %s", format(margin_abs, trim = TRUE))]
  df[, final_lab  := paste0("Final thr = ", final_thr)]
  df[, interim_lab := paste0("Interim thr=", interim_thr)]
  
  # Pareto frontier within each facet (margin x final_thr)
  setorder(df, margin_abs, final_lab, alpha)
  df[, frontier := power == cummax(power), by = .(margin_abs, final_lab)]
  
  # best feasible (alpha <= target) to annotate
  best_tbl <- if (!is.null(cal$feasible) && nrow(cal$feasible) > 0L) {
    bt <- as.data.table(cal$feasible)
    setorder(bt, -power, ExpN_alt)
    bt[, margin_lab := sprintf("Delta = %s", format(margin_abs, trim = TRUE))]
    bt[, final_lab  := paste0("Final thr = ", final_thr)]
    bt[1:min(label_top_n, .N)]
  } else data.table()
  
  p <- ggplot(df, aes(x = alpha, y = power,
                      color = margin_lab,
                      shape = factor(interim_thr))) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = target_alpha, linetype = 2, linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.9) +
    # frontier lines
    geom_line(data = df[frontier == TRUE],
              aes(group = interaction(margin_lab, final_lab)),
              linewidth = 0.7, alpha = 0.8) +
    facet_wrap(~ final_lab) +
    labs(
      x = "Type I error (alpha)",
      y = "Power",
      color = "Success margin",
      shape = "Interim threshold",
      title = "Design calibration: Power vs Type I",
      subtitle = sprintf("Vertical dashed line = target alpha (%g)", target_alpha)
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  # Annotate top feasible designs (if any)
  if (nrow(best_tbl) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        data = best_tbl,
        aes(x = alpha, y = power,
            label = paste0("Delta=", margin_abs,
                           "; Int=", interim_thr,
                           "; Fin=", final_thr,
                           "; N~", round(ExpN_alt, 1))),
        size = 3, label.size = 0.15, alpha = 0.9,
        fill = "white", label.padding = unit(0.15, "lines"),
        max.overlaps = Inf
      )
    } else {
      p <- p + geom_text(
        data = best_tbl,
        aes(x = alpha, y = power,
            label = paste0("Delta=", margin_abs,
                           ", Int=", interim_thr,
                           ", Fin=", final_thr)),
        size = 3, vjust = -0.8
      )
    }
  }
  
  # Optional percent formatting if 'scales' is available
  if (requireNamespace("scales", quietly = TRUE)) {
    p <- p +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, 1))
  } else {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }
  
  p
}

# --- Pick best calibration row and bake into args ---
#' Adopt a calibrated design configuration
#'
#' Picks a selected row from `grid_calibrate()` and returns updated arguments
#' plus a two-scenario list (null vs alternative) for subsequent exploration.
#'
#' @param cal Output from `grid_calibrate()`.
#' @param base_args Baseline argument list.
#' @param null_med Control-arm median under the null.
#' @param alt_med Experimental median under the alternative.
#' @param which Integer index specifying which row of `cal$top` to adopt.
#'
#' @return A list containing the updated arguments (`args_star`), the selected
#'   row (`pick`), and a two-scenario list (`scens2`).
#' @export
adopt_calibration <- function(cal, base_args, null_med, alt_med, which = 1L) {
  stopifnot(is.list(cal), !is.null(cal$top), nrow(cal$top) >= which)
  best <- data.table::as.data.table(cal$top)[which]
  
  args_star <- base_args
  # success thresholds from calibration
  args_star$efficacy_threshold_hc_prob     <- best$interim_thr
  args_star$final_success_posterior_prob_threshold <- best$final_thr
  # set Triplet superiority margin (absolute months)
  args_star$median_pfs_success_threshold_arms <- c(
    Doublet = null_med,
    Triplet = null_med + best$margin_abs
  )
  # keep predictive off during calibration/exploration (safer for Type I)
  args_star$predictive_fast <- FALSE
  
  # 2-scenario grid: all-null vs Triplet=alt
  scens2 <- scenarios_from_grid(list(
    weibull_median_true_arms = list(
      c(Doublet = null_med, Triplet = null_med),
      c(Doublet = null_med, Triplet = alt_med)
    )
  ))
  
  list(args_star = args_star, pick = best, scens2 = scens2)
}

# --- Explore early-stopping & information gates around a calibrated design ---
# Sweeps: futility threshold, min events, min median FU, and calendar look beat
# Explore early stopping while varying the futility comparator for *all* experimental arms
# Comparator options:
#   base = "null+delta": uses P(median < null + delta | data) >= fut_thr  => stop
#   base = "alt"       : uses P(median < alt         | data) >= fut_thr  => stop
#
# Notes:
# - Uses adopt_calibration() to import the calibrated success/margin from `cal` (your prior step).
# - Keeps predictive_fast = FALSE while exploring early-stopping knobs (clean Type I accounting).
#
# --- Explore early-stopping & information gates around a calibrated design ---
# Adds sweeps over:
#   * per-arm gates (min_events_per_arm, min_median_followup_per_arm, min_person_time_frac_per_arm)
#   * schedule type: "calendar" (beats) OR "persontime" (milestones with a calendar backstop)
#
# For person-time schedules, pass a LIST of milestone vectors via `pt_milestones_choices`,
# e.g. list(c(0.30,0.45,0.60,0.80,1.00), c(0.30,0.60,0.90)).
#
#' Explore early stopping knobs around a calibrated design
#'
#' Sweeps futility thresholds, information gates, and interim schedules around a
#' calibrated design to characterise trade-offs between alpha, power, and PETs.
#'
#' @param cal Calibration output from `grid_calibrate()` or similar.
#' @param base_args Baseline argument list.
#' @param null_med Control-arm median under null.
#' @param alt_med Experimental median under alternative.
#' @param base Futility baseline mode: "null+delta" or "alt".
#' @param futility_delta_grid Numeric vector of deltas added to null medians.
#' @param fut_thr_grid Numeric vector of futility probability thresholds.
#' @param min_events_grid Global minimum event counts to consider.
#' @param min_medFU_grid Global minimum median follow-up values.
#' @param schedule_modes Character vector: "calendar", "persontime", or both.
#' @param beat_grid Calendar beat schedules to evaluate (months).
#' @param pt_milestones_choices List of person-time milestone vectors.
#' @param latest_calendar_look_grid Backstop calendar times for person-time schedules.
#' @param min_events_per_arm_grid Per-arm minimum events for gating.
#' @param min_median_followup_per_arm_grid Per-arm minimum median follow-up.
#' @param min_person_time_frac_per_arm_grid Per-arm minimum person-time fractions.
#' @param sims Number of simulations per configuration.
#' @param seed RNG seed.
#' @param parallel Logical; use parallel processing.
#'
#' @return data.table of operating characteristics for each configuration.
#' @export
explore_early_stopping_from_cal <- function(
    cal,
    base_args,
    null_med,
    alt_med,
    base                  = c("null+delta","alt"),
    futility_delta_grid   = c(0, 1, 2, 3),    # only used when base = "null+delta"
    fut_thr_grid          = c(0.6, 0.7, 0.8, 0.9),
    
    # legacy/global gates
    min_events_grid       = c(12, 18),
    min_medFU_grid        = c(3, 4.5),
    
    # schedule sweep
    schedule_modes        = c("calendar","persontime"),
    beat_grid             = c(3, 6),          # used when schedule="calendar"
    
    # person-time schedule knobs
    pt_milestones_choices = list(c(0.30,0.45,0.60,0.80,1.00)),
    latest_calendar_look_grid = c(Inf),       # e.g., c(12, 18) months as a backstop
    
    # per-arm gates (NEW)
    min_events_per_arm_grid          = c(8, 12),
    min_median_followup_per_arm_grid = c(0, 4.5),
    min_person_time_frac_per_arm_grid= c(0.00, 0.25),
    
    sims                 = 400,
    seed                 = 123,
    parallel             = (.Platform$OS.type == "unix")
) {
  base <- match.arg(base)
  stopifnot(is.list(cal), !is.null(cal$top) || !is.null(cal$feasible))
  
  # 1) import calibrated success thresholds/margin & the two scenarios (null, alt)
  adopted  <- adopt_calibration(cal, base_args, null_med = null_med, alt_med = alt_med, which = 1L)
  args0    <- adopted$args_star
  scens2   <- adopted$scens2
  exp_arms <- exp_arms_from_args(args0)
  
  # helper to stringify milestone vectors for the output table
  fmt_frac_vec <- function(x) if (is.null(x)) "NULL" else paste0(sprintf("%.2f", x), collapse = ",")
  
  # build the list of combinations manually (because milestones is a list-column)
  combos <- list()
  for (ft in fut_thr_grid) {
    # choose the futility comparator set for exp arms
    fut_deltas <- if (base == "null+delta") futility_delta_grid else 0
    for (fd in fut_deltas) {
      for (ev in min_events_grid) {
        for (mu in min_medFU_grid) {
          for (me in min_events_per_arm_grid) {
            for (mfu in min_median_followup_per_arm_grid) {
              for (mpt in min_person_time_frac_per_arm_grid) {
                
                # schedule: calendar beats
                if ("calendar" %in% schedule_modes) {
                  for (bt in beat_grid) {
                    combos[[length(combos) + 1L]] <- data.table::data.table(
                      schedule = "calendar",
                      fut_base = base,
                      fut_delta = fd,
                      fut_thr = ft,
                      min_events = ev,
                      min_medFU = mu,
                      beat = bt,
                      pt_milestones = NA_character_,
                      latest_calendar_look = NA_real_,
                      min_events_per_arm = me,
                      min_median_followup_per_arm = mfu,
                      min_person_time_frac_per_arm = mpt
                    )
                  }
                }
                
                # schedule: person-time milestones
                if ("persontime" %in% schedule_modes) {
                  for (ml in seq_along(pt_milestones_choices)) {
                    mlv <- pt_milestones_choices[[ml]]
                    stopifnot(is.numeric(mlv), all(mlv > 0 & mlv <= 1))
                    for (lc in latest_calendar_look_grid) {
                      combos[[length(combos) + 1L]] <- data.table::data.table(
                        schedule = "persontime",
                        fut_base = base,
                        fut_delta = fd,
                        fut_thr = ft,
                        min_events = ev,
                        min_medFU = mu,
                        beat = NA_real_,
                        pt_milestones = fmt_frac_vec(mlv),
                        latest_calendar_look = lc,
                        min_events_per_arm = me,
                        min_median_followup_per_arm = mfu,
                        min_person_time_frac_per_arm = mpt
                      )
                    }
                  }
                }
                
              } # mpt
            } # mfu
          } # me
        } # mu
      } # ev
    } # fd
  } # ft
  
  combos_dt <- data.table::rbindlist(combos, use.names = TRUE)
  out <- vector("list", nrow(combos_dt))
  
  # 2) iterate & simulate
  for (i in seq_len(nrow(combos_dt))) {
    row <- combos_dt[i]
    
    args_i <- args0
    
    # set the *interim futility* comparator for experimental arms
    args_i <- set_futility_medians(
      args    = args_i,
      null_med = null_med,
      alt_med  = alt_med,
      base     = row$fut_base,
      delta    = row$fut_delta
    )
    args_i$futility_threshold_hc_prob <- row$fut_thr
    
    # global gates
    args_i$min_events_hc <- row$min_events
    args_i$min_median_followup_hc <- row$min_medFU
    
    # per-arm gates
    args_i$min_events_per_arm             <- row$min_events_per_arm
    args_i$min_median_followup_per_arm    <- row$min_median_followup_per_arm
    args_i$min_person_time_frac_per_arm   <- row$min_person_time_frac_per_arm
    
    # schedule
    if (row$schedule == "calendar") {
      args_i$interim_calendar_beat   <- row$beat
      args_i$person_time_milestones  <- NULL
      args_i$latest_calendar_look    <- Inf
    } else {
      # person-time milestones
      mlv <- as.numeric(strsplit(row$pt_milestones, ",")[[1]])
      args_i$person_time_milestones <- mlv
      args_i$latest_calendar_look   <- row$latest_calendar_look
      # keep a (small) calendar beat as a guard if you like, but we'll let the backstop rule dominate
      args_i$interim_calendar_beat  <- args0$interim_calendar_beat %||% 2
    }
    
    # clean exploration (predictive off)
    args_i$predictive_fast <- FALSE
    args_i$num_simulations <- sims
    
    # run null vs alt
    res <- run_scenarios(args_i, adopted$scens2, parallel = parallel, seed = seed)
    
    r_null <- res[res$scenario == 1 & res$Arm_Name %in% exp_arms, ]
    r_alt  <- res[res$scenario == 2 & res$Arm_Name %in% exp_arms, ]
    
    out[[i]] <- data.table::data.table(
      schedule     = row$schedule,
      fut_base     = row$fut_base,
      fut_delta    = row$fut_delta,
      fut_thr      = row$fut_thr,
      min_events   = row$min_events,
      min_medFU    = row$min_medFU,
      beat         = ifelse(row$schedule == "calendar", row$beat, NA_real_),
      pt_milestones = ifelse(row$schedule == "persontime", row$pt_milestones, NA_character_),
      latest_calendar_look = ifelse(row$schedule == "persontime", row$latest_calendar_look, NA_real_),
      
      min_events_per_arm           = row$min_events_per_arm,
      min_median_followup_per_arm  = row$min_median_followup_per_arm,
      min_person_time_frac_per_arm = row$min_person_time_frac_per_arm,
      
      alpha        = mean(r_null$Type_I_Error_or_Power),
      power        = mean(r_alt$Type_I_Error_or_Power),
      ExpN_null    = mean(r_null$Exp_N),
      ExpN_alt     = mean(r_alt$Exp_N),
      PET_Eff_null = mean(r_null$PET_Efficacy),
      PET_Eff_alt  = mean(r_alt$PET_Efficacy),
      PET_Fut_null = mean(r_null$PET_Futility),
      PET_Fut_alt  = mean(r_alt$PET_Futility)
    )
  }
  
  early <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  data.table::setorder(
    early,
    schedule, fut_base, fut_delta, fut_thr,
    min_events, min_medFU,
    min_events_per_arm, min_median_followup_per_arm, min_person_time_frac_per_arm,
    beat, pt_milestones, latest_calendar_look
  )
  early[]
}


# --- Quick plot (optional) ---
#' Plot early-stopping trade-offs
#'
#' Visualises power versus type I error for the early-stopping exploration grid.
#'
#' @param early_df data.table/data.frame returned by
#'   `explore_early_stopping_from_cal()`.
#' @param target_alpha Type I error target to show as a reference line.
#' @param fix_min_ev Optional scalar to filter `min_events`.
#' @param fix_mfu Optional scalar to filter `min_medFU`.
#' @param fix_beat Optional scalar to filter calendar beats.
#'
#' @return A ggplot object (or the filtered data when `ggplot2` is unavailable).
#' @export
plot_early_tradeoff <- function(early_df,
                                target_alpha = 0.10,
                                fix_min_ev = NULL,
                                fix_mfu    = NULL,
                                fix_beat   = NULL) {
  df <- data.table::as.data.table(early_df)
  if (!is.null(fix_min_ev)) df <- df[min_events == fix_min_ev]
  if (!is.null(fix_mfu))    df <- df[min_medFU  == fix_mfu]
  if (!is.null(fix_beat))   df <- df[beat       == fix_beat]
  if (nrow(df) == 0L) stop("No rows after filtering--relax the fixed settings.")
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(df, ggplot2::aes(alpha, power, color = factor(fut_thr))) +
      ggplot2::geom_vline(xintercept = target_alpha, linetype = 2, linewidth = 0.4) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_path(ggplot2::aes(group = fut_thr), linewidth = 0.6, alpha = 0.7) +
      ggplot2::facet_grid(min_events ~ min_medFU, labeller = "label_both") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "Type I error", y = "Power", color = "Futility thr")
  } else {
    message("Install ggplot2 to plot. Returning data.")
    df
  }
}


# --- Filter + rank early-stopping designs by constraints ---
#' Filter early-stopping designs by operating targets
#'
#' @param early_df data.table/data.frame produced by
#'   `explore_early_stopping_from_cal()`.
#' @param alpha_cap Maximum acceptable type I error.
#' @param power_floor Minimum acceptable power.
#'
#' @return data.table containing the subset that meets the supplied criteria.
#' @export
filter_early_grid <- function(early_df,
                              alpha_cap   = 0.10,
                              power_floor = 0.70) {
  df <- data.table::as.data.table(early_df)
  ok <- df[alpha <= alpha_cap & power >= power_floor]
  if (nrow(ok) == 0L) return(ok)
  
  # Rank: smallest N under alt, then more futility under null, then less efficacy under null
  data.table::setorder(ok, ExpN_alt, -PET_Fut_null, PET_Eff_null, fut_thr, min_events, min_medFU, beat)
  ok[]
}


# --- Pick the single recommended design under constraints ---
#' Recommend a single early-stopping design
#'
#' Selects the top-performing row from an early-stopping grid subject to
#' alpha/power (and optionally PET) constraints.
#'
#' @param df data.table/data.frame from `explore_early_stopping_from_cal()`.
#' @param alpha_cap Maximum acceptable type I error.
#' @param power_floor Minimum acceptable power.
#' @param pet_fut_cap Optional cap on alternative PET for futility.
#'
#' @return data.table row describing the recommended design.
#' @export
recommend_design_from_early <- function(df,
                                        alpha_cap = 0.10,
                                        power_floor = 0.80,
                                        pet_fut_cap = NULL) {
  df_filt <- df[alpha <= alpha_cap & power >= power_floor]
  if (!is.null(pet_fut_cap)) {
    df_filt <- df_filt[PET_Fut_alt <= pet_fut_cap]
  }
  if (nrow(df_filt) == 0)
    stop("No designs meet specified criteria.")
  
  df_filt[order(-power, PET_Fut_alt, ExpN_alt)][1]
}


# --- Bake the recommended early-stopping knobs back into args ---
#' Apply a recommended early-stopping configuration to the argument list
#'
#' @param args_star Baseline argument list (typically from `adopt_calibration()`).
#' @param rec_row Single-row data.table produced by
#'   `recommend_design_from_early()`.
#'
#' @return Modified argument list with the recommended early-stopping settings.
#' @export
apply_recommended_to_args <- function(args_star, rec_row) {
  stopifnot(nrow(rec_row) == 1L)
  a <- args_star
  a$futility_threshold_hc_prob <- rec_row$fut_thr
  a$min_events_hc         <- rec_row$min_events
  a$min_median_followup_hc             <- rec_row$min_medFU
  a$interim_calendar_beat           <- rec_row$beat
  
  # NEW: carry per-arm gates if present
  if ("min_events_per_arm" %in% names(rec_row)) {
    a$min_events_per_arm <- rec_row$min_events_per_arm
  }
  if ("min_median_followup_per_arm" %in% names(rec_row)) {
    a$min_median_followup_per_arm <- rec_row$min_median_followup_per_arm
  }
  if ("min_person_time_frac_per_arm" %in% names(rec_row)) {
    a$min_person_time_frac_per_arm <- rec_row$min_person_time_frac_per_arm
  }
  a
}


# All non-reference arms are considered "experimental"
#' Extract experimental arm names from an argument list
#'
#' @param args Argument list containing `arm_names` and `reference_arm_name`.
#'
#' @return Character vector of experimental arm names.
#' @export
exp_arms_from_args <- function(args) {
  setdiff(args$arm_names, args$reference_arm_name)
}


# Modify args so interim *futility* compares to either:
#   - null + delta  (base = "null+delta"; delta is a single number for all exp arms)
#   - alt           (base = "alt")
# The HC/reference arm keeps its "null" futility threshold.
#' Adjust futility medians for experimental arms
#'
#' @param args Argument list whose `futility_median_arms` entry will be updated.
#' @param null_med Numeric scalar representing the null median.
#' @param alt_med Numeric scalar representing the alternative median.
#' @param base Character scalar selecting `"null+delta"` or `"alt"`.
#' @param delta Numeric offset added when `base = "null+delta"`.
#'
#' @return Modified argument list with updated `futility_median_arms`.
#' @export
set_futility_medians <- function(args, null_med, alt_med, base = c("null+delta","alt"), delta = 0) {
  base <- match.arg(base)
  arms_exp <- exp_arms_from_args(args)
  fut <- args$futility_median_arms
  if (base == "null+delta") {
    fut[arms_exp] <- null_med + delta
  } else if (base == "alt") {
    fut[arms_exp] <- alt_med
  }
  args$futility_median_arms <- fut
  args
}
