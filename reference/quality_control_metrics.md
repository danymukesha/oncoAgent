# Compute Quality Control Metrics

Calculates per-sample QC metrics for RNA-seq data including library
size, detection rate, and outlier flags.

## Usage

``` r
quality_control_metrics(count_matrix, col_data = NULL)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix

- col_data:

  Sample metadata data.frame

## Value

data.frame with QC metrics per sample
