# Run edgeR Differential Expression

Performs differential expression using edgeR's quasi-likelihood pipeline
with TMM normalization.

## Usage

``` r
run_edger(
  count_matrix,
  col_data,
  contrast_levels = NULL,
  padj_threshold = 0.05
)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix (integer)

- col_data:

  Sample metadata with 'condition' column

- contrast_levels:

  Character vector: c("tumor", "normal")

- padj_threshold:

  Adjusted p-value threshold (default: 0.05)

## Value

data.frame with DE results
