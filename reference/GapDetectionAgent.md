# Gap Detection Agent

R6 class that identifies data gaps, statistical weaknesses, and
conflicting signals across evidence dimensions.

## Public fields

- `state`:

  A list storing the internal state of the agent, including:

  - `thresholds`: List of minimum thresholds for analysis (samples,
    DEGs, modules)

  - `gaps`: List of detected gaps with type, severity, and description

  - `recommendations`: Character vector of recommendations based on gaps

  - `severity`: Overall severity of detected gaps
    (critical/warning/info)

## Methods

### Public methods

- [`GapDetectionAgent$new()`](#method-GapDetectionAgent-new)

- [`GapDetectionAgent$detect()`](#method-GapDetectionAgent-detect)

- [`GapDetectionAgent$summary()`](#method-GapDetectionAgent-summary)

- [`GapDetectionAgent$needs_rerun()`](#method-GapDetectionAgent-needs_rerun)

- [`GapDetectionAgent$clone()`](#method-GapDetectionAgent-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new GapDetectionAgent

#### Usage

    GapDetectionAgent$new(min_samples = 50, min_degs = 100, min_modules = 3)

#### Arguments

- `min_samples`:

  Minimum sample size for adequate power

- `min_degs`:

  Minimum DEGs expected

- `min_modules`:

  Minimum network modules expected

------------------------------------------------------------------------

### Method `detect()`

Detect gaps in analysis results

#### Usage

    GapDetectionAgent$detect(results, data)

#### Arguments

- `results`:

  Named list of pipeline results

- `data`:

  Named list of input data

#### Returns

List of detected gaps with severity and recommendations

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Get gap summary

#### Usage

    GapDetectionAgent$summary()

#### Returns

List with gap statistics

------------------------------------------------------------------------

### Method `needs_rerun()`

Check if pipeline should be re-run

#### Usage

    GapDetectionAgent$needs_rerun()

#### Returns

Logical, TRUE if critical gaps exist

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    GapDetectionAgent$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
