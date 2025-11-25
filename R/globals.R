#' @keywords internal
#' @noRd
utils::globalVariables(c(
  # data.table special symbols
  ".N", ".SD", ".",

  # columns you create via := or use in i/j
  "event_status", "observed_time", "interval_num", "events", "person_time",
  "alpha", "power", "ExpN_alt", "PET_Fut_null", "PET_Fut_alt", "PET_Eff_null", "fut_thr",
  "min_events", "min_medFU", "beat", "pt_milestones", "latest_calendar_look",
  "margin_abs", "interim_thr", "final_thr", "interim_lab", "final_lab",
  "frontier", "True_Median", "Exp_N", "Exp_Events", "Pr_Reach_Max_N", "Type_I_Error_or_Power",
  "PET_Efficacy", "PET_Futility", "Pr_Final_Efficacy", "Pr_Final_Futility",
  "scenario", "Arm_Name", "margin_lab", "min_events_per_arm",
  "min_median_followup_per_arm", "min_person_time_frac_per_arm",
  "schedule", "fut_base", "fut_delta"
))
