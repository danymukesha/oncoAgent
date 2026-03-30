# Rank DE Genes by Combined Score

Ranks DE results by a combination of statistical significance and effect
size.

## Usage

``` r
rank_de_genes(de_results, w_pvalue = 0.5, w_lfc = 0.5)
```

## Arguments

- de_results:

  data.frame with columns: gene, log2FoldChange, padj

- w_pvalue:

  Weight for p-value component (default: 0.5)

- w_lfc:

  Weight for log-fold-change component (default: 0.5)

## Value

data.frame with added rank and combined_score columns
