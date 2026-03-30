# Set oncoAgent Configuration

Sets global configuration options for the oncoAgent package.

## Usage

``` r
set_oncoagent_config(
  cache_dir = NULL,
  n_cores = 1,
  verbose = TRUE,
  log_level = "INFO"
)
```

## Arguments

- cache_dir:

  Default cache directory

- n_cores:

  Number of parallel cores to use

- verbose:

  Logical, enable verbose logging

- log_level:

  Logging level: "INFO", "DEBUG", "WARN", "ERROR"

## Value

Invisibly returns the configuration list
