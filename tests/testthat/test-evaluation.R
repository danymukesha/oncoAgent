# test-evaluation.R
# Tests for the evaluation module

test_that("evaluate_targets computes precision correctly", {
  predicted <- c("TP53", "BRCA1", "GENE1", "GENE2")
  known <- c("TP53", "BRCA1", "EGFR", "KRAS")

  precision <- evaluate_targets(predicted, known)

  expect_equal(precision, 0.5)
})

test_that("precision_at_k returns correct values", {
  scored <- data.frame(
    gene = c("TP53", "BRCA1", "GENE1", "EGFR", "KRAS", "GENE2", "BRAF"),
    composite_score = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
    stringsAsFactors = FALSE
  )
  known <- c("TP53", "BRCA1", "EGFR", "KRAS", "BRAF")

  p3 <- precision_at_k(scored, known, k = 3)
  p5 <- precision_at_k(scored, known, k = 5)

  expect_equal(p3, 2/3)
  expect_equal(p5, 4/5)
})

test_that("recall_of_known_targets computes correctly", {
  predicted <- c("TP53", "BRCA1", "EGFR")
  known <- c("TP53", "BRCA1", "EGFR", "KRAS", "BRAF")

  recall <- recall_of_known_targets(predicted, known)

  expect_equal(recall, 3/5)
})

test_that("calculate_auroc returns valid value", {
  set.seed(42)
  scored <- data.frame(
    gene = paste0("Gene_", 1:100),
    composite_score = c(runif(10, 0.5, 1), runif(90, 0, 0.5)),
    stringsAsFactors = FALSE
  )

  known <- paste0("Gene_", 1:10)

  auroc <- calculate_auroc(scored, known)

  expect_true(is.numeric(auroc))
  expect_true(auroc >= 0 & auroc <= 1)
  expect_true(auroc > 0.5)
})

test_that("benchmark_vs_goldstandard produces complete metrics", {
  scored <- data.frame(
    gene = c("TP53", "BRCA1", "EGFR", "KRAS", paste0("G", seq_len(46))),
    composite_score = seq(1, 0, length.out = 50),
    stringsAsFactors = FALSE
  )

  gs <- list(genes = c("TP53", "BRCA1", "EGFR", "KRAS", "BRAF"))

  bench <- benchmark_vs_goldstandard(scored, gs, n_permutations = 100)

  expect_true(is.list(bench))
  expect_true("precision" %in% names(bench))
  expect_true("recall" %in% names(bench))
  expect_true("f1" %in% names(bench))
  expect_true("permutation_pvalue" %in% names(bench))
  expect_true(bench$precision > 0)
})
