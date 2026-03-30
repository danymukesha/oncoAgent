# Network Analysis Module

Protein-protein interaction network construction, community detection,
centrality analysis, and disease module identification.

Constructs an igraph network from STRING PPI data or a custom edge list.
Filters by score and removes disconnected components.

## Usage

``` r
build_ppi_network(
  ppi_data,
  score_threshold = 400,
  remove_isolates = TRUE,
  directed = FALSE
)
```

## Arguments

- ppi_data:

  data.frame with columns: from, to, score

- score_threshold:

  Minimum interaction score (default: 400)

- remove_isolates:

  Logical, remove unconnected nodes (default: TRUE)

- directed:

  Logical, create directed graph (default: FALSE)

## Value

igraph object
