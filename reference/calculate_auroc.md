# Calculate AUROC

Computes Area Under the Receiver Operating Characteristic curve for
target prediction.

## Usage

``` r
calculate_auroc(scored_targets, known_targets)
```

## Arguments

- scored_targets:

  data.frame with columns: gene, composite_score

- known_targets:

  Character vector of known gold standard genes

## Value

AUROC value (0-1)
