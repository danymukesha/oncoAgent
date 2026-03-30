# Analysis Agent

R6 class that executes the analytical pipeline (DE, survival, network
analysis) and manages analysis state.

## Public fields

- `state`:

  A list storing the internal state of the agent, including:

  - `de_methods`: Character vector of differential expression methods.

  - `padj_threshold`: Adjusted p-value threshold.

  - `lfc_threshold`: Log-fold-change threshold.

  - `results`: List of results from different analysis steps.

  - `execution_log`: List of log messages with timestamps.

  - `completed_steps`: Character vector of completed steps.

## Methods

### Public methods

- [`AnalysisAgent$new()`](#method-AnalysisAgent-new)

- [`AnalysisAgent$run()`](#method-AnalysisAgent-run)

- [`AnalysisAgent$status()`](#method-AnalysisAgent-status)

- [`AnalysisAgent$clone()`](#method-AnalysisAgent-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new AnalysisAgent

#### Usage

    AnalysisAgent$new(
      de_methods = c("deseq2", "edger", "limma"),
      padj_threshold = 0.05,
      lfc_threshold = 0
    )

#### Arguments

- `de_methods`:

  Character vector of DE methods

- `padj_threshold`:

  Adjusted p-value threshold

- `lfc_threshold`:

  Log-fold-change threshold

------------------------------------------------------------------------

### Method `run()`

Run full analysis pipeline

#### Usage

    AnalysisAgent$run(data)

#### Arguments

- `data`:

  Named list of data (expression, clinical, ppi, etc.)

#### Returns

Named list of all analysis results

------------------------------------------------------------------------

### Method `status()`

Get analysis status

#### Usage

    AnalysisAgent$status()

#### Returns

List with analysis status

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnalysisAgent$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
