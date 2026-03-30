# Run oncoAgent Pipeline

Executes the full agentic pipeline for oncology target discovery. This
is the main entry point for running the system.

## Usage

``` r
run_oncopipeline(config, existing_data = NULL, gold_standard = NULL)
```

## Arguments

- config:

  Pipeline configuration (from create_pipeline_config)

- existing_data:

  Optional pre-loaded data to skip ingestion

- gold_standard:

  Optional gold standard gene set for evaluation

## Value

List with: targets, gaps, evaluation, agent_states, config
