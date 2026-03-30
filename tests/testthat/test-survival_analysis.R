# test-survival_analysis.R
# Tests for the survival analysis module

test_that("run_cox_regression fits Cox PH model", {
  set.seed(42)
  n <- 200

  gene_expr <- rnorm(n, 5, 2)
  time <- rexp(n, rate = 0.01)
  status <- rbinom(n, 1, 0.5)

  result <- run_cox_regression(gene_expr, time, status, gene_name = "TP53")

  expect_true(is.list(result))
  expect_equal(result$gene, "TP53")
  expect_true("hazard_ratio" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("concordance" %in% names(result))
  expect_true(result$n_samples == n)
})

test_that("run_kaplan_meier produces survival curves", {
  set.seed(42)
  n <- 200

  gene_expr <- rnorm(n, 5, 2)
  time <- rexp(n, rate = 0.01)
  status <- rbinom(n, 1, 0.5)

  result <- run_kaplan_meier(gene_expr, time, status, gene_name = "BRCA1")

  expect_true(is.list(result))
  expect_true("fit" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("plot_data" %in% names(result))
  expect_true("cutoff" %in% names(result))
})

test_that("calculate_hazard_ratios processes multiple genes", {
  set.seed(42)
  n_genes <- 20
  n_samples <- 100

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene_", 1:n_genes)
  colnames(expr) <- paste0("Sample_", 1:n_samples)

  time <- rexp(n_samples, 0.01)
  status <- rbinom(n_samples, 1, 0.5)

  hr_table <- calculate_hazard_ratios(expr, time, status)

  expect_true(is.data.frame(hr_table))
  expect_true("gene" %in% names(hr_table))
  expect_true("hr" %in% names(hr_table))
  expect_true("pvalue" %in% names(hr_table))
  expect_true("padj" %in% names(hr_table))
  expect_true(nrow(hr_table) > 0)
})
