# test-differential_expression.R
# Tests for the differential expression module

test_that("run_deseq2 detects differential expression", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 20

  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)

  counts[1:20, 1:10] <- rnbinom(200, mu = 500, size = 10)

  col_data <- data.frame(
    condition = factor(rep(c("tumor", "normal"), each = 10)),
    row.names = colnames(counts)
  )

  de <- run_deseq2(counts, col_data)

  expect_true(is.data.frame(de))
  expect_true("gene" %in% names(de))
  expect_true("log2FoldChange" %in% names(de))
  expect_true("padj" %in% names(de))
  expect_true(nrow(de) > 0)
  expect_true(sum(de$padj < 0.05, na.rm = TRUE) > 0)
})

test_that("run_edger produces valid results", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 20

  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)

  counts[1:20, 1:10] <- rnbinom(200, mu = 500, size = 10)

  col_data <- data.frame(
    condition = factor(rep(c("tumor", "normal"), each = 10)),
    row.names = colnames(counts)
  )

  de <- run_edger(counts, col_data)

  expect_true(is.data.frame(de))
  expect_true("log2FoldChange" %in% names(de))
  expect_true("padj" %in% names(de))
})

test_that("run_limma_voom produces valid results", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 20

  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)

  counts[1:20, 1:10] <- rnbinom(200, mu = 500, size = 10)

  col_data <- data.frame(
    condition = factor(rep(c("tumor", "normal"), each = 10)),
    row.names = colnames(counts)
  )

  de <- run_limma_voom(counts, col_data)

  expect_true(is.data.frame(de))
  expect_true("log2FoldChange" %in% names(de) || "logFC" %in% names(de))
})

test_that("rank_de_genes assigns proper rankings", {
  de_results <- data.frame(
    gene = paste0("Gene_", 1:100),
    log2FoldChange = rnorm(100, 0, 2),
    padj = runif(100, 0, 0.1),
    stringsAsFactors = FALSE
  )

  ranked <- rank_de_genes(de_results)

  expect_true("combined_score" %in% names(ranked))
  expect_true("rank" %in% names(ranked))
  expect_equal(min(ranked$rank), 1)
  expect_equal(max(ranked$rank), 100)
})
