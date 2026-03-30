# Consensus Differential Expression

Runs multiple DE methods and computes a consensus ranking. Genes
detected by multiple methods are ranked higher.

## Usage

``` r
consensus_de(
  count_matrix,
  col_data,
  methods = c("deseq2", "edger", "limma"),
  padj_threshold = 0.05,
  min_methods = 2
)
```

## Arguments

- count_matrix:

  Gene x Sample count matrix

- col_data:

  Sample metadata with 'condition' column

- methods:

  Character vector of methods: "deseq2", "edger", "limma"

- padj_threshold:

  Adjusted p-value threshold

- min_methods:

  Minimum number of methods that must agree (default: 2)

## Value

data.frame with consensus results
