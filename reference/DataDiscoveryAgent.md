# Data Discovery Agent

R6 class that checks data availability, identifies missing modalities,
and recommends data sources for a given cancer project.

## Public fields

- `state`:

  Current agent state

## Methods

### Public methods

- [`DataDiscoveryAgent$new()`](#method-DataDiscoveryAgent-new)

- [`DataDiscoveryAgent$assess()`](#method-DataDiscoveryAgent-assess)

- [`DataDiscoveryAgent$fetch_missing()`](#method-DataDiscoveryAgent-fetch_missing)

- [`DataDiscoveryAgent$status()`](#method-DataDiscoveryAgent-status)

- [`DataDiscoveryAgent$clone()`](#method-DataDiscoveryAgent-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new DataDiscoveryAgent

#### Usage

    DataDiscoveryAgent$new(
      project = "TCGA-BRCA",
      required_modalities = c("expression", "clinical", "mutation", "ppi"),
      cache_dir = NULL
    )

#### Arguments

- `project`:

  TCGA project identifier

- `required_modalities`:

  Character vector of required data types

- `cache_dir`:

  Cache directory for data

------------------------------------------------------------------------

### Method `assess()`

Assess what data is available and what is missing

#### Usage

    DataDiscoveryAgent$assess(existing_data = list())

#### Arguments

- `existing_data`:

  Named list of already-fetched data

#### Returns

Updated state with availability assessment

------------------------------------------------------------------------

### Method `fetch_missing()`

Fetch missing data automatically

#### Usage

    DataDiscoveryAgent$fetch_missing(existing_data = list())

#### Arguments

- `existing_data`:

  Named list of current data

#### Returns

Named list with all fetched data

------------------------------------------------------------------------

### Method `status()`

Get agent status report

#### Usage

    DataDiscoveryAgent$status()

#### Returns

List with current status

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DataDiscoveryAgent$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
