# Calculate Network Centrality Metrics

Computes degree, betweenness, closeness, eigenvector, and PageRank
centrality for all nodes.

## Usage

``` r
calculate_network_centrality(
  g,
  metrics = c("degree", "betweenness", "closeness", "pagerank")
)
```

## Arguments

- g:

  igraph network object

- metrics:

  Character vector of metrics to compute

## Value

data.frame with centrality scores per node
