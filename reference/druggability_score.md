# Druggability Score

Estimates druggability of target based on protein class, known drug
interactions, and structural features.

## Usage

``` r
druggability_score(genes, drug_data = NULL)
```

## Arguments

- genes:

  Character vector of gene symbols

- drug_data:

  GDSC or similar drug sensitivity data

## Value

data.frame with druggability scores
