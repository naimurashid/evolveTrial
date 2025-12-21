# Apply futility action to trial state

Apply futility action to trial state

## Usage

``` r
apply_futility_action(state, arm, action = "drop_arm")
```

## Arguments

- state:

  Current trial state

- arm:

  Arm that hit futility

- action:

  Action to take: "stop_trial", "drop_arm", or "continue"

## Value

Updated state and action taken
