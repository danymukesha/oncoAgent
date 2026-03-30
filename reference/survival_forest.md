# Survival Random Forest

Runs Cox regression across multiple genes simultaneously using penalized
regression (LASSO) for variable selection.

## Usage

``` r
survival_forest(expr_matrix, time, status, alpha = 1, n_folds = 10)
```

## Arguments

- expr_matrix:

  Gene x Sample expression matrix

- time:

  Numeric vector of survival times

- status:

  Integer vector of event status

- alpha:

  Elastic net mixing parameter (1 = LASSO, 0 = Ridge)

- n_folds:

  Cross-validation folds (default: 10)

## Value

List with: selected_genes, coefficients, lambda, cv_error
