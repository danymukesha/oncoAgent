# Evaluation Agent

R6 class that benchmarks predicted targets against gold standard gene
sets and computes evaluation metrics.

## Public fields

- `state`:

  A list storing the internal state of the agent, including:

  - `thresholds`: List of minimum thresholds (samples, DEGs, modules).

  - `gaps`: List of detected gaps with type, severity, and description.

  - `recommendations`: Character vector of recommendations based on
    gaps.

  - `severity`: Overall severity of detected gaps
    (critical/warning/info).

## Methods

### Public methods

- [`EvaluationAgent$new()`](#method-EvaluationAgent-new)

- [`EvaluationAgent$evaluate()`](#method-EvaluationAgent-evaluate)

- [`EvaluationAgent$report()`](#method-EvaluationAgent-report)

- [`EvaluationAgent$clone()`](#method-EvaluationAgent-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new EvaluationAgent

#### Usage

    EvaluationAgent$new(gold_standard = NULL, k_values = c(10, 20, 50))

#### Arguments

- `gold_standard`:

  Named list of gene sets for benchmarking

- `k_values`:

  K values for Precision@K (default: c(10, 20, 50))

------------------------------------------------------------------------

### Method `evaluate()`

Evaluate targets against gold standards

#### Usage

    EvaluationAgent$evaluate(scored_targets, gold_standard = NULL)

#### Arguments

- `scored_targets`:

  Output from composite_target_score or rank_targets

- `gold_standard`:

  Named list of gene sets (if NULL, uses constructor value)

#### Returns

List of evaluation metrics

------------------------------------------------------------------------

### Method `report()`

Generate evaluation report

#### Usage

    EvaluationAgent$report()

#### Returns

Formatted evaluation report

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    EvaluationAgent$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
