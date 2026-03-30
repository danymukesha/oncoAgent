# Fetch GDSC Drug Sensitivity Data

Downloads IC50 and AUC drug sensitivity data from the Genomics of Drug
Sensitivity in Cancer (GDSC) database.

## Usage

``` r
fetch_gdsc_drug_sensitivity(
  dataset = "GDSC2",
  metric = "IC50",
  cache_dir = NULL
)
```

## Arguments

- dataset:

  Which GDSC dataset: "GDSC1", "GDSC2", or "both"

- metric:

  Sensitivity metric: "IC50" or "AUC"

- cache_dir:

  Directory for caching

## Value

data.frame with drug sensitivity measurements per cell line
