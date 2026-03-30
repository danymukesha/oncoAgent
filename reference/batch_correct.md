# Batch Correction

Applies batch correction using ComBat (sva) or limma's
removeBatchEffect.

## Usage

``` r
batch_correct(expr_matrix, batch_factor, method = "combat", covariates = NULL)
```

## Arguments

- expr_matrix:

  Normalized expression matrix (genes x samples)

- batch_factor:

  Factor or vector indicating batch membership

- method:

  "combat" or "limma"

- covariates:

  Optional data.frame of biological covariates to preserve

## Value

Batch-corrected expression matrix
