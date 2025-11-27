# Resolve gate parameter vector for specified arms with optional scaling

Internal helper that converts gate parameters (which may be scalar,
named vector, or positional vector) into a named numeric vector for the
specified arms. Optionally applies proportional scaling based on
randomization probabilities.

## Usage

``` r
resolve_gate_vec(
  raw,
  target_arms,
  all_arm_names,
  randomization_probs = NULL,
  default = 0,
  scale = FALSE
)
```

## Arguments

- raw:

  Raw gate parameter value (scalar, named vector, or positional vector).

- target_arms:

  Character vector of arm names to resolve gates for.

- all_arm_names:

  Character vector of all arm names in the trial (for positional
  matching).

- randomization_probs:

  Optional named numeric vector of randomization probabilities.

- default:

  Default value if `raw` is NULL (default 0).

- scale:

  Logical; if TRUE and `raw` is scalar/NULL, apply proportional scaling
  by dividing each arm's randomization probability by the maximum
  probability (default FALSE).

## Value

Named numeric vector with one element per arm in `target_arms`.
