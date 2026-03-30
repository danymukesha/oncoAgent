# Detect Network Modules (Communities)

Identifies functional modules in the PPI network using Louvain community
detection.

## Usage

``` r
detect_network_modules(g, method = "louvain", resolution = 1)
```

## Arguments

- g:

  igraph network object

- method:

  Community detection: "louvain", "walktrap", "infomap"

- resolution:

  Resolution parameter for Louvain (default: 1.0)

## Value

List with: membership, modularity, module_genes, module_sizes
