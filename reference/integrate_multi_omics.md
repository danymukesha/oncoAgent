# Multi-Omics Integration

Integrates gene expression, mutation, methylation, and CRISPR data into
a unified feature matrix for target scoring.

## Usage

``` r
integrate_multi_omics(
  expression,
  mutations = NULL,
  methylation = NULL,
  crispr = NULL,
  integration_method = "weighted_average"
)
```

## Arguments

- expression:

  Gene x Sample expression matrix

- mutations:

  Mutation data.frame (gene, sample, consequence)

- methylation:

  Gene x Sample methylation matrix (beta values)

- crispr:

  Gene dependency scores

- integration_method:

  "concatenate", "weighted_average", or "factor"

## Value

Integrated gene x feature matrix
