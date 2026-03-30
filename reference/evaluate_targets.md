# Evaluation Module

Benchmarking functions for comparing predicted targets against gold
standard gene sets. Computes AUROC, precision@K, recall, and generates
comprehensive evaluation reports.

Computes precision of predicted targets against a known set.

## Usage

``` r
evaluate_targets(predicted, known_targets)
```

## Arguments

- predicted:

  Character vector of predicted gene symbols

- known_targets:

  Character vector of known gold standard genes

## Value

Precision value (0-1)
