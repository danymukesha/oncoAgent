# Fetch PubMed Abstracts

Queries PubMed for abstracts related to specified genes and cancer type.
Returns structured metadata and abstract text.

## Usage

``` r
fetch_pubmed_abstracts(
  genes,
  cancer_type = NULL,
  max_results = 100,
  date_from = "2015/01/01",
  cache_dir = NULL
)
```

## Arguments

- genes:

  Character vector of gene symbols

- cancer_type:

  Cancer type for focused search (e.g., "breast cancer")

- max_results:

  Maximum number of results per gene (default: 100)

- date_from:

  Start date filter (YYYY/MM/DD format)

- cache_dir:

  Directory for caching

## Value

data.frame with columns: pmid, title, abstract, genes, date
