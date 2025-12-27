# Create hybrid design parameter structure

Create hybrid design parameter structure

## Usage

``` r
create_hybrid_theta(
  eff_sa = 0.9,
  fut_sa = 0.1,
  hr_threshold_sa = 0.8,
  ev_sa = 15,
  nmax_sa = 40,
  conversion_trigger = c("any_single_success", "all_single_success", "k_of_K"),
  k_required = 1,
  pp_go = 0.7,
  pp_nogo = 0.2,
  ss_method = c("predictive", "posterior"),
  max_additional_n = 60,
  n_add_candidates = seq(10, 100, by = 10),
  eff_ba = 0.975,
  fut_ba = 0.05,
  ev_ba = 15,
  nmax_ba = 80,
  futility_action = c("drop_arm", "stop_trial", "continue"),
  prior_strength = 0.5,
  n_outer = 1000,
  n_inner = 250
)
```

## Arguments

- eff_sa:

  SA efficacy threshold (default 0.90)

- fut_sa:

  SA futility threshold (default 0.10)

- hr_threshold_sa:

  Target HR vs historical (default 0.80)

- ev_sa:

  Minimum events for SA interim (default 15)

- nmax_sa:

  Maximum N in SA phase per arm (default 40)

- conversion_trigger:

  Trigger type: "any_single_success", "all_single_success", "k_of_K"

- k_required:

  For k_of_K trigger, number required (default 1)

- pp_go:

  PP threshold to proceed to BA (default 0.70)

- pp_nogo:

  PP threshold to stop (default 0.20)

- ss_method:

  SSR method: "predictive" or "posterior" (default "predictive")

- max_additional_n:

  Maximum additional patients for BA (default 60)

- n_add_candidates:

  Candidate N values for PP curve (default seq(10, 100, 10))

- eff_ba:

  BA efficacy threshold (default 0.975)

- fut_ba:

  BA futility threshold (default 0.05)

- ev_ba:

  Minimum events for BA interim (default 15)

- nmax_ba:

  Maximum N per arm in BA phase (default 80)

- futility_action:

  Action on SA futility (default "drop_arm")

- prior_strength:

  Gamma prior concentration (default 0.5)

- n_outer:

  MC outer samples for PP (default 1000)

- n_inner:

  MC inner samples for PP (default 250)

## Value

Named list of hybrid design parameters

## Note

LIMITATION: This implementation assumes exactly 2 arms (1 reference + 1
experimental). Key places that assume 2 arms:

- Event gate uses ev_ba \* 2 (line ~422)

- BA comparison takes first non-reference arm only (lines ~652, ~849)

- See docs/TWO_ARM_ASSUMPTIONS.md in adaptive-trial-bo-paper for full
  audit.
