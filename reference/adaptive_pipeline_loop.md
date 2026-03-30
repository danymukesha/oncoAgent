# Adaptive Pipeline Loop

Runs the pipeline iteratively, adjusting configuration based on gap
detection results. This is the "agentic iteration" loop.

## Usage

``` r
adaptive_pipeline_loop(
  config,
  max_iterations = 5,
  improvement_threshold = 0.01
)
```

## Arguments

- config:

  Initial pipeline configuration

- max_iterations:

  Maximum iterations (default: 5)

- improvement_threshold:

  Minimum improvement to continue (default: 0.01)

## Value

List with final results and iteration history
