# Kaplan-Meier Survival Analysis

Fits Kaplan-Meier survival curves and performs log-rank test.

## Usage

``` r
run_kaplan_meier(
  gene_expression,
  time,
  status,
  gene_name = NULL,
  split_method = "median"
)
```

## Arguments

- gene_expression:

  Numeric vector of expression values

- time:

  Numeric vector of survival times

- status:

  Integer vector of event status

- gene_name:

  Gene symbol

- split_method:

  Dichotomization method: "median", "quartile", "optimal"

## Value

List with: fit, pvalue, group_survival, plot_data
