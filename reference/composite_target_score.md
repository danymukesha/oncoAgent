# Target Scoring Module

Multi-dimensional target scoring integrating differential expression,
survival, network centrality, druggability, and CRISPR dependency data.

Computes a weighted composite score for each gene across multiple
evidence dimensions.

## Usage

``` r
composite_target_score(
  de_results,
  survival_hr = NULL,
  centrality = NULL,
  crispr_scores = NULL,
  weights = NULL
)
```

## Arguments

- de_results:

  data.frame with columns: gene, log2FoldChange, padj

- survival_hr:

  data.frame with columns: gene, hr, pvalue (from
  calculate_hazard_ratios)

- centrality:

  data.frame with columns: gene, composite_centrality

- crispr_scores:

  Optional data.frame with columns: gene, dependency_score

- weights:

  Named list of weights for each evidence type: de, survival,
  centrality, crispr (default: equal weights)

## Value

data.frame with composite scores and per-dimension breakdowns
