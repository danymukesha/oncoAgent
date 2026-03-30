# Map Gene Identifiers

Converts between gene identifiers (Ensembl, Entrez, symbol) using
biomaRt.

## Usage

``` r
map_gene_ids(
  ids,
  from_type = "hgnc_symbol",
  to_type = "ensembl_gene_id",
  mart_dataset = "hsapiens_gene_ensembl"
)
```

## Arguments

- ids:

  Character vector of gene identifiers

- from_type:

  Source ID type: "ensembl_gene_id", "entrezgene_id", "hgnc_symbol"

- to_type:

  Target ID type

- mart_dataset:

  Ensembl dataset (default: "hsapiens_gene_ensembl")

## Value

data.frame with original and mapped IDs
