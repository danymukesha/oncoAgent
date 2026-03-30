# Run limma-voom Differential Expression

Performs differential expression using limma-voom with sample quality
weights.

## Usage

``` r
run_limma_voom(
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
