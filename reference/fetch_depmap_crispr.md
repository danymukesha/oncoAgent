# Fetch DepMap CRISPR Screen Data

Downloads CRISPR dependency scores from the DepMap portal. Returns gene
effect scores indicating essentiality across cell lines.

## Usage

``` r
fetch_depmap_crispr(release_version = NULL, cache_dir = NULL)
```

## Arguments

- release_version:

  DepMap release version (e.g., "24Q2"). NULL for latest.

- cache_dir:

  Directory for caching

## Value

data.frame with cell lines as rows, genes as columns, values as
dependency scores (more negative = more essential)
