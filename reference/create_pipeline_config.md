# Orchestration Module

Pipeline orchestration with adaptive agentic loop. Manages the full
workflow from data ingestion through target scoring with gap-driven
re-execution.

Creates a standardized configuration list for the oncoAgent pipeline.

## Usage

``` r
create_pipeline_config(
  project = "TCGA-BRCA",
  cancer_type = "breast cancer",
  de_methods = c("deseq2", "edger", "limma"),
  padj_threshold = 0.05,
  max_iterations = 5,
  output_dir = "results"
)
```

## Arguments

- project:

  TCGA project identifier

- cancer_type:

  Cancer type label for reference panel

- de_methods:

  DE methods to use

- padj_threshold:

  Adjusted p-value threshold

- max_iterations:

  Maximum agentic loop iterations

- output_dir:

  Directory for results

## Value

Configuration list
