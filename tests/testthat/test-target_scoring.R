# test-target_scoring.R
# Tests for the target scoring module

test_that("composite_target_score integrates multiple evidence types", {
  de <- data.frame(
    gene = paste0("Gene_", 1:50),
    log2FoldChange = rnorm(50, 0, 2),
    padj = c(runif(25, 0, 0.01), runif(25, 0.01, 0.1)),
    stringsAsFactors = FALSE
  )

  surv <- data.frame(
    gene = paste0("Gene_", 1:30),
    hr = runif(30, 0.5, 2),
    pvalue = runif(30, 0, 0.05),
    stringsAsFactors = FALSE
  )

  cent <- data.frame(
    gene = paste0("Gene_", 1:40),
    composite_centrality = runif(40, 0, 1),
    stringsAsFactors = FALSE
  )

  scores <- composite_target_score(de, surv, cent)

  expect_true(is.data.frame(scores))
  expect_true("gene" %in% names(scores))
  expect_true("composite_score" %in% names(scores))
  expect_true("de_score" %in% names(scores))
  expect_true("hr_score" %in% names(scores))
  expect_true("centrality_score" %in% names(scores))
  finite_scores <- scores$composite_score[is.finite(scores$composite_score)]
  expect_true(length(finite_scores) > 0)
  expect_true(all(finite_scores >= 0))
  expect_true(all(finite_scores <= 1))
})

test_that("rank_targets returns proper tiers", {
  scored <- data.frame(
    gene = paste0("Gene_", 1:100),
    composite_score = runif(100, 0, 1),
    de_score = runif(100, 0, 1),
    hr_score = runif(100, 0, 1),
    centrality_score = runif(100, 0, 1),
    crispr_score = runif(100, 0, 1),
    stringsAsFactors = FALSE
  )

  ranked <- rank_targets(scored, top_n = 30)

  expect_true("target_tier" %in% names(ranked))
  expect_true("Tier 1" %in% ranked$target_tier)
  expect_true(nrow(ranked) <= 30)
})

test_that("druggability_score classifies known druggable genes", {
  genes <- c("EGFR", "BRAF", "TP53", "BRCA1", "PIK3CA", "EZH2", "UNKNOWN")

  drug_scores <- druggability_score(genes)

  expect_true(is.data.frame(drug_scores))
  expect_equal(nrow(drug_scores), 7)
  expect_true(drug_scores$druggability_score[drug_scores$gene == "EGFR"] > 0.5)
  expect_true(drug_scores$is_kinase[drug_scores$gene == "BRAF"])
})
