# test-preprocessing.R
# Tests for the preprocessing module

test_that("normalize_counts handles vst normalization", {
  set.seed(42)
  counts <- matrix(rnbinom(15000, mu = 100, size = 10), nrow = 1500, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:1500)
  colnames(counts) <- paste0("Sample_", 1:10)

  col_data <- data.frame(
    condition = factor(rep(c("tumor", "normal"), each = 5)),
    row.names = colnames(counts)
  )

  norm <- normalize_counts(counts, method = "vst", col_data = col_data)

  expect_true(is.matrix(norm))
  expect_equal(nrow(norm), 1500)
  expect_equal(ncol(norm), 10)
  expect_false(any(is.na(norm)))
})

test_that("normalize_counts handles TMM normalization", {
  set.seed(42)
  counts <- matrix(rnbinom(1000, mu = 100, size = 10), nrow = 100, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:100)
  colnames(counts) <- paste0("Sample_", 1:10)

  norm <- normalize_counts(counts, method = "tmm")

  expect_true(is.matrix(norm))
  expect_equal(nrow(norm), 100)
})

test_that("filter_low_expression removes low-count genes", {
  set.seed(42)
  counts <- matrix(0, nrow = 50, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:50)
  colnames(counts) <- paste0("Sample_", 1:10)

  counts[1:10, ] <- rnbinom(100, mu = 100, size = 10)
  counts[11:50, ] <- rnbinom(400, mu = 0.1, size = 1)

  filtered <- filter_low_expression(counts, min_count = 10, min_samples = 0.1)

  expect_true(nrow(filtered) <= 50)
  expect_true(nrow(filtered) >= 5)
})

test_that("quality_control_metrics computes correct metrics", {
  set.seed(42)
  counts <- matrix(rnbinom(500, mu = 100, size = 10), nrow = 50, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:50)
  colnames(counts) <- paste0("Sample_", 1:10)

  qc <- quality_control_metrics(counts)

  expect_equal(nrow(qc), 10)
  expect_true("library_size" %in% names(qc))
  expect_true("n_genes_detected" %in% names(qc))
  expect_true("is_outlier" %in% names(qc))
  expect_true(all(qc$library_size > 0))
})

test_that("harmonize_clinical standardizes columns", {
  clinical <- data.frame(
    bcr_patient_barcode = paste0("TCGA-AB-", 1:10),
    vital_status = c(rep("Dead", 5), rep("Alive", 5)),
    days_to_last_follow_up = runif(10, 100, 3000),
    stringsAsFactors = FALSE
  )

  harmonized <- harmonize_clinical(clinical)

  expect_true("status" %in% names(harmonized))
  expect_true("time" %in% names(harmonized))
  expect_true("sample_id" %in% names(harmonized))
  expect_equal(sum(harmonized$status), 5)
})
