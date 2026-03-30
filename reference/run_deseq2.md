# Differential Expression Module

Unified interface for differential expression analysis using DESeq2,
edgeR, and limma-voom. Includes consensus ranking across methods.

Performs differential expression analysis using DESeq2 with shrinkage
estimation and multiple testing correction.

## Usage

``` r
run_deseq2(
  count_matrix,
  col_data,
  contrast = NULL,
  lfc_threshold = 0,
  padj_threshold = 0.05,
  shrink = TRUE
)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix (integer)

- col_data:

  Sample metadata with 'condition' column

- contrast:

  Character vector of length 3: c("condition", "tumor", "normal")

- lfc_threshold:

  Log2 fold-change threshold for significance (default: 0)

- padj_threshold:

  Adjusted p-value threshold (default: 0.05)

- shrink:

  Logical, apply apeglm shrinkage (default: TRUE)

## Value

data.frame with DE results sorted by adjusted p-value
