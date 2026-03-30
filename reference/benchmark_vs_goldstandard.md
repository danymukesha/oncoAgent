# Benchmark Against Gold Standard

Comprehensive benchmarking against multiple gold standard gene sets with
stratification.

## Usage

``` r
benchmark_vs_goldstandard(scored_targets, gold_standard, n_permutations = 1000)
```

## Arguments

- scored_targets:

  Output from composite_target_score or rank_targets

- gold_standard:

  Named list with \$genes and optionally \$per_source

- n_permutations:

  Number of permutation tests for significance

## Value

List with precision, recall, F1, permutation p-values
