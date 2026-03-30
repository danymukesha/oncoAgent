# Data Ingestion Module

Functions for fetching multi-omics data from public repositories (TCGA,
DepMap, GDSC, STRING, PubMed) with caching, error handling, and
standardized output formats.

Downloads and prepares TCGA data including gene expression, clinical,
and mutation data for a specified cancer project.

## Usage

``` r
fetch_tcga_data(
  project = "TCGA-BRCA",
  data_types = c("expression", "clinical", "mutation"),
  cache_dir = NULL,
  force_download = FALSE
)
```

## Arguments

- project:

  TCGA project identifier (e.g., "TCGA-BRCA", "TCGA-LUAD")

- data_types:

  Character vector of data types to fetch. Options: "expression",
  "clinical", "mutation", "methylation". Default: all.

- cache_dir:

  Directory for caching downloaded data

- force_download:

  Logical, force re-download even if cached

## Value

A named list containing requested data modalities as
SummarizedExperiment or data.frame objects
