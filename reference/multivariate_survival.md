# Multivariate Survival Analysis

Fits a multivariate Cox model with gene expression and clinical
covariates.

## Usage

``` r
multivariate_survival(expr_matrix, clinical_df, clinical_vars = NULL)
```

## Arguments

- expr_matrix:

  Gene x Sample expression matrix (top genes)

- clinical_df:

  data.frame with time, status, and clinical covariates

- clinical_vars:

  Character vector of clinical variable names to include

## Value

List with: model, coefficients, forest_plot_data
