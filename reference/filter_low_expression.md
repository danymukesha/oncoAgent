# Filter Low-Expression Genes

Removes genes with low expression across samples to reduce noise and
improve statistical power.

## Usage

``` r
filter_low_expression(count_matrix, min_count = 10, min_samples = 0.1)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix

- min_count:

  Minimum count threshold (default: 10)

- min_samples:

  Minimum proportion of samples that must exceed min_count (default:
  0.1)

## Value

Filtered count matrix
