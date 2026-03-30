# Build Reference Panel from Multiple Gold Standards

Aggregates known cancer gene sets from COSMIC Cancer Gene Census, Open
Targets, and other sources into a unified reference for benchmarking.

## Usage

``` r
build_reference_panel(
  sources = c("cosmic", "cgc"),
  cancer_type = NULL,
  cache_dir = NULL
)
```

## Arguments

- sources:

  Character vector: "cosmic", "opentargets", "intogen", "cgc"

- cancer_type:

  Filter to specific cancer type (NULL for all)

- cache_dir:

  Directory for caching

## Value

A list with: \$genes (character vector), \$metadata (data.frame),
\$per_source (named list of gene sets)
