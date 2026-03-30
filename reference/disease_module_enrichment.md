# Disease Module Enrichment

Tests whether network modules are enriched for known disease genes using
Fisher's exact test.

## Usage

``` r
disease_module_enrichment(modules, disease_genes, universe = NULL)
```

## Arguments

- modules:

  Module detection output from detect_network_modules()

- disease_genes:

  Character vector of known disease gene symbols

- universe:

  Character vector of all genes considered (background)

## Value

data.frame with enrichment results per module
