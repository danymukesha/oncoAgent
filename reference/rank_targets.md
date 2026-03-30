# Rank Targets

Final ranking combining composite scores with additional filters and
annotations.

## Usage

``` r
rank_targets(
  scored_targets,
  top_n = 50,
  min_de_padj = 0.05,
  min_evidence_dims = 1
)
```

## Arguments

- scored_targets:

  Output from composite_target_score()

- top_n:

  Return top N targets (default: 50)

- min_de_padj:

  Minimum DE significance (default: 0.05)

- min_evidence_dims:

  Minimum evidence dimensions required (default: 1)

## Value

data.frame of ranked targets
