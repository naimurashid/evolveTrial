# Compute single-arm posterior probability

P(HR \< c \| data) where HR = lambda_arm / lambda_hist

## Usage

``` r
compute_p_single_arm(post_a, post_b, hist_hazard, hr_threshold, base_args)
```

## Arguments

- post_a:

  Posterior shape parameters

- post_b:

  Posterior rate parameters

- hist_hazard:

  Historical hazard

- hr_threshold:

  Target HR threshold

- base_args:

  evolveTrial base configuration

## Value

Posterior probability
