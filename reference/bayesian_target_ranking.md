# Bayesian Target Ranking

Bayesian approach to target ranking that models uncertainty in evidence
and computes posterior probabilities of being a true target.

## Usage

``` r
bayesian_target_ranking(
  scored_targets,
  prior_probability = 0.01,
  noise_sd = 0.1
)
```

## Arguments

- scored_targets:

  Output from composite_target_score()

- prior_probability:

  Prior probability of being a true target (default: 0.01)

- noise_sd:

  Assumed noise level in scores (default: 0.1)

## Value

data.frame with posterior probabilities
