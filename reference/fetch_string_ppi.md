# Fetch Protein-Protein Interaction Network from STRING

Retrieves PPI data from the STRING database for a set of genes.

## Usage

``` r
fetch_string_ppi(
  genes,
  species = 9606,
  score_threshold = 400,
  cache_dir = NULL
)
```

## Arguments

- genes:

  Character vector of gene symbols

- species:

  Species taxon ID (default: 9606 for human)

- score_threshold:

  Minimum combined score (0-1000, default: 400)

- cache_dir:

  Directory for caching

## Value

data.frame with columns: from, to, score, and sub-scores
