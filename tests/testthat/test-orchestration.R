# test-orchestration.R
# Tests for the orchestration module

test_that("create_pipeline_config returns valid config", {
  config <- create_pipeline_config(
    project = "TCGA-BRCA",
    cancer_type = "breast cancer"
  )

  expect_true(is.list(config))
  expect_equal(config$project, "TCGA-BRCA")
  expect_equal(config$cancer_type, "breast cancer")
  expect_true("de_methods" %in% names(config))
  expect_true("max_iterations" %in% names(config))
  expect_true(config$max_iterations > 0)
})

test_that("run_oncopipeline executes with synthetic data", {
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

  clinical <- data.frame(
    sample_id = colnames(counts),
    time = rexp(n_samples, 0.01),
    status = rbinom(n_samples, 1, 0.5),
    stringsAsFactors = FALSE
  )

  existing_data <- list(
    expression = counts,
    col_data = col_data,
    clinical = clinical
  )

  config <- create_pipeline_config(
    project = "TCGA-BRCA",
    output_dir = tempdir()
  )

  result <- run_oncopipeline(config, existing_data = existing_data)

  expect_true(is.list(result))
  expect_true("targets" %in% names(result))
  expect_true("gaps" %in% names(result))
  expect_true("config" %in% names(result))
  expect_true(!is.null(result$targets))
  expect_true(nrow(result$targets) > 0)
})

test_that("export_results creates output files", {
  result <- list(
    targets = data.frame(gene = paste0("Gene_", 1:10),
                         composite_score = runif(10),
                         stringsAsFactors = FALSE),
    hypotheses = data.frame(gene = paste0("Gene_", 1:5),
                            stringsAsFactors = FALSE),
    evaluation = list(precision = 0.5, recall = 0.3),
    config = list(project = "test"),
    timestamp = Sys.time()
  )

  out_dir <- file.path(tempdir(), "test_export")
  files <- export_results(result, out_dir)

  expect_true(length(files) > 0)
  expect_true(file.exists(file.path(out_dir, "ranked_targets.csv")))
  expect_true(file.exists(file.path(out_dir, "hypotheses.csv")))
  expect_true(file.exists(file.path(out_dir, "config.json")))

  unlink(out_dir, recursive = TRUE)
})
