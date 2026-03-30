# Export Network for Visualization

Exports the network in formats compatible with Cytoscape or other
network visualization tools.

## Usage

``` r
export_network(g, output_file, node_attributes = NULL, format = "graphml")
```

## Arguments

- g:

  igraph network object

- output_file:

  Output file path (.graphml, .sif, or .csv)

- node_attributes:

  data.frame of node attributes to include

- format:

  Export format: "graphml", "sif", "edgelist"

## Value

Invisibly returns the file path
