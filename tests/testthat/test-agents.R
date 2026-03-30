# test-agents.R
# Tests for the agentic layer

test_that("DataDiscoveryAgent assesses available data", {
  agent <- DataDiscoveryAgent$new(
    project = "TCGA-BRCA",
    required_modalities = c("expression", "clinical", "mutation")
  )

  mock_data <- list(
    expression = matrix(rnorm(1000), nrow = 100, ncol = 10),
    clinical = data.frame(sample_id = paste0("S", 1:10))
  )

  state <- agent$assess(mock_data)

  expect_true("available_modalities" %in% names(state))
  expect_true("missing_modalities" %in% names(state))
  expect_true("expression" %in% state$available_modalities)
  expect_true("mutation" %in% state$missing_modalities)
  expect_true(length(state$recommendations) > 0)
})

test_that("AnalysisAgent runs pipeline and tracks state", {
  agent <- AnalysisAgent$new(de_methods = c("deseq2", "limma"))

  set.seed(42)
  n_genes <- 100
  n_samples <- 20
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)

  counts[1:10, 1:10] <- rnbinom(100, mu = 500, size = 10)

  col_data <- data.frame(
    condition = factor(rep(c("tumor", "normal"), each = 10)),
    row.names = colnames(counts)
  )

  data <- list(expression = counts, col_data = col_data)
  results <- agent$run(data)

  expect_true(is.list(results))
  expect_true("de" %in% names(results))
  expect_true("scoring" %in% names(results))

  status <- agent$status()
  expect_true("de" %in% status$completed_steps)
})

test_that("GapDetectionAgent identifies critical gaps", {
  agent <- GapDetectionAgent$new()

  results <- list(de = NULL, survival = NULL, network = NULL, scoring = NULL)
  data <- list()

  gaps <- agent$detect(results, data)

  expect_true(is.list(gaps))
  expect_true(length(gaps) > 0)
  expect_true(any(sapply(gaps, function(g) g$severity == "critical")))

  summary <- agent$summary()
  expect_true(summary$severity == "critical")
  expect_true(agent$needs_rerun())
})

test_that("GapDetectionAgent reports no critical gaps with complete data", {
  agent <- GapDetectionAgent$new(min_samples = 10, min_degs = 50, min_modules = 2)

  de <- data.frame(
    gene = paste0("Gene_", 1:200),
    padj = runif(200, 0, 0.05),
    log2FoldChange = rnorm(200, 0, 2),
    n_methods = rep(3, 200),
    stringsAsFactors = FALSE
  )

  scoring <- data.frame(
    gene = paste0("Gene_", 1:200),
    composite_score = runif(200, 0, 1),
    stringsAsFactors = FALSE
  )

  surv <- data.frame(
    gene = paste0("Gene_", 1:100),
    hr = runif(100, 0.5, 2),
    pvalue = runif(100, 0, 0.05),
    stringsAsFactors = FALSE
  )

  net <- list(
    network = list(n = 50),
    centrality = data.frame(gene = paste0("Gene_", 1:50),
                            composite_centrality = runif(50), stringsAsFactors = FALSE)
  )

  results <- list(de = de, survival = surv, network = net, scoring = scoring)
  data <- list(
    expression = matrix(rnorm(2000), nrow = 200, ncol = 10),
    clinical = data.frame(sample_id = paste0("S", 1:10))
  )

  gaps <- agent$detect(results, data)
  gap_summary <- agent$summary()

  expect_true(is.list(gap_summary))
  expect_true("n_gaps" %in% names(gap_summary))
  expect_true("severity" %in% names(gap_summary))
})

test_that("HypothesisAgent generates hypotheses from scored targets", {
  agent <- HypothesisAgent$new(top_n_hypotheses = 10)

  scored <- data.frame(
    gene = paste0("Gene_", 1:50),
    composite_score = runif(50, 0, 1),
    de_score = runif(50, 0, 1),
    hr_score = runif(50, 0, 1),
    centrality_score = runif(50, 0, 1),
    crispr_score = runif(50, 0, 1),
    stringsAsFactors = FALSE
  )

  hypotheses <- agent$generate(scored)

  expect_true(is.data.frame(hypotheses))
  expect_true(nrow(hypotheses) <= 10)
  expect_true("mechanism" %in% names(hypotheses))
  expect_true("confidence" %in% names(hypotheses))
})

test_that("EvaluationAgent computes metrics against gold standard", {
  agent <- EvaluationAgent$new()

  scored <- data.frame(
    gene = c("TP53", "BRCA1", "EGFR", "KRAS", "PIK3CA",
             paste0("Gene_", 6:50)),
    composite_score = seq(1, 0, length.out = 50),
    stringsAsFactors = FALSE
  )

  gs <- list(
    genes = c("TP53", "BRCA1", "EGFR", "KRAS", "BRAF", "PTEN", "RB1"),
    per_source = list(
      cosmic = c("TP53", "BRCA1", "EGFR", "KRAS"),
      cgc = c("BRAF", "PTEN", "RB1")
    )
  )

  metrics <- agent$evaluate(scored, gs)

  expect_true(is.list(metrics))
  expect_true("overall" %in% names(metrics))
  expect_true(metrics$overall > 0)

  report <- agent$report()
  expect_true(is.character(report))
  expect_true(grepl("Precision", report))
})
