# Preprocessing Module

Quality control, normalization, batch correction, and data harmonization
functions for multi-omics preprocessing.

Applies normalization to RNA-seq count matrices. Supports DESeq2
variance-stabilizing transformation (VST), regularized log (rlog), and
edgeR TMM normalization.

## Usage

``` r
normalize_counts(count_matrix, method = "vst", col_data = NULL, blind = TRUE)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix (integer)

- method:

  Normalization method: "vst", "rlog", "tmm", "cpm"

- col_data:

  Sample metadata data.frame with condition column

- blind:

  Logical, blind VST/rlog to experimental design (default: TRUE)

## Value

Normalized expression matrix
