# Generate Evaluation Report

Creates a comprehensive markdown evaluation report comparing predicted
targets against gold standards.

## Usage

``` r
generate_evaluation_report(scored_targets, gold_standard, output_file = NULL)
```

## Arguments

- scored_targets:

  Output from composite_target_score

- gold_standard:

  Named list with \$genes and \$per_source

- output_file:

  Optional file path to write report

## Value

Character string with the report
