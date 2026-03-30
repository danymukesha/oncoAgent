# Survival Analysis Module

Cox regression, Kaplan-Meier, and multivariate survival analysis for
linking gene expression to patient outcomes.

Fits univariate Cox PH model for a single gene's expression against
survival outcomes.

## Usage

``` r
run_cox_regression(
  gene_expression,
  time,
  status,
  gene_name = NULL,
  continuous = TRUE,
  median_split = FALSE
)
```

## Arguments

- gene_expression:

  Numeric vector of expression values for one gene

- time:

  Numeric vector of survival times

- status:

  Integer vector of event status (1 = event, 0 = censored)

- gene_name:

  Gene symbol for labeling

- continuous:

  Logical, treat expression as continuous (default: TRUE)

- median_split:

  Logical, dichotomize at median (default: FALSE)

## Value

List with: model, summary, hazard_ratio, pvalue, concordance
