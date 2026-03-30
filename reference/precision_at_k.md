# Precision at K

Computes precision for the top K predicted targets.

## Usage

``` r
precision_at_k(scored_targets, known_targets, k = 10)
```

## Arguments

- scored_targets:

  data.frame with columns: gene, composite_score

- known_targets:

  Character vector of known gold standard genes

- k:

  Number of top predictions to consider

## Value

Precision value (0-1)
