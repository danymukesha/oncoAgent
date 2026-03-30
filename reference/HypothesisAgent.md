# Hypothesis Agent

R6 class that generates and ranks hypotheses about target mechanisms,
combining multi-omics evidence with literature signals.

## Public fields

- `state`:

  A list storing the internal state of the agent, including:

  - `top_n`: Maximum number of hypotheses to generate.

  - `min_score`: Minimum evidence score threshold.

  - `hypotheses`: List of generated hypotheses.

  - `evidence_log`: List of evidence collected for each hypothesis.

## Methods

### Public methods

- [`HypothesisAgent$new()`](#method-HypothesisAgent-new)

- [`HypothesisAgent$generate()`](#method-HypothesisAgent-generate)

- [`HypothesisAgent$top()`](#method-HypothesisAgent-top)

- [`HypothesisAgent$clone()`](#method-HypothesisAgent-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new HypothesisAgent

#### Usage

    HypothesisAgent$new(top_n_hypotheses = 20, min_evidence_score = 0.1)

#### Arguments

- `top_n_hypotheses`:

  Maximum hypotheses to generate

- `min_evidence_score`:

  Minimum evidence score threshold

------------------------------------------------------------------------

### Method `generate()`

Generate hypotheses from analysis results

#### Usage

    HypothesisAgent$generate(
      scored_targets,
      literature = NULL,
      network_results = NULL
    )

#### Arguments

- `scored_targets`:

  Output from composite_target_score

- `literature`:

  Optional literature data.frame

- `network_results`:

  Optional network analysis results

#### Returns

data.frame of ranked hypotheses

------------------------------------------------------------------------

### Method `top()`

Get top hypotheses

#### Usage

    HypothesisAgent$top(n = 10)

#### Arguments

- `n`:

  Number to return

#### Returns

data.frame of top hypotheses

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    HypothesisAgent$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
