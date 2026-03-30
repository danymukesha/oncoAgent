# Calculate Hazard Ratios for Multiple Genes

Runs univariate Cox regression for each gene and returns a summary table
of hazard ratios.

## Usage

``` r
calculate_hazard_ratios(expr_matrix, time, status, genes = NULL, n_cores = 1)
```

## Arguments

- expr_matrix:

  Gene x Sample expression matrix

- time:

  Numeric vector of survival times

- status:

  Integer vector of event status

- genes:

  Character vector of genes to test (NULL = all)

- n_cores:

  Number of parallel cores (default: 1)

## Value

data.frame with hazard ratios, CIs, and p-values per gene
