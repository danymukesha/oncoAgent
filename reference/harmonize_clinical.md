# Harmonize Clinical Data

Standardizes clinical data from TCGA into a consistent format for
downstream survival and subgroup analyses.

## Usage

``` r
harmonize_clinical(
  clinical,
  vital_status_col = NULL,
  time_col = NULL,
  time_unit = "days"
)
```

## Arguments

- clinical:

  Raw clinical data.frame (from TCGAbiolinks or similar)

- vital_status_col:

  Column name for vital status

- time_col:

  Column name for follow-up time

- time_unit:

  Unit of time: "days", "months", "years"

## Value

Harmonized data.frame with standardized column names
