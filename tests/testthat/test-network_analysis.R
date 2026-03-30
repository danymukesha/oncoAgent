# test-network_analysis.R
# Tests for the network analysis module

test_that("build_ppi_network creates valid igraph object", {
  ppi <- data.frame(
    from = c("A", "B", "C", "D", "E"),
    to = c("B", "C", "D", "E", "A"),
    score = c(900, 800, 700, 600, 500),
    stringsAsFactors = FALSE
  )

  g <- build_ppi_network(ppi, score_threshold = 400)

  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), 5)
  expect_equal(igraph::ecount(g), 5)
})

test_that("detect_network_modules finds communities", {
  set.seed(42)
  edges <- data.frame(
    from = c("A", "B", "C", "D", "E", "F", "G", "H"),
    to = c("B", "C", "A", "E", "F", "D", "H", "G"),
    score = rep(800, 8),
    stringsAsFactors = FALSE
  )

  g <- build_ppi_network(edges, score_threshold = 400)
  modules <- detect_network_modules(g)

  expect_true(is.list(modules))
  expect_true("membership" %in% names(modules))
  expect_true("modularity" %in% names(modules))
  expect_true("module_genes" %in% names(modules))
  expect_true(modules$n_modules >= 1)
})

test_that("calculate_network_centrality computes all metrics", {
  set.seed(42)
  edges <- data.frame(
    from = c("A", "A", "A", "B", "C", "D"),
    to = c("B", "C", "D", "C", "D", "E"),
    score = rep(800, 6),
    stringsAsFactors = FALSE
  )

  g <- build_ppi_network(edges, score_threshold = 400)
  centrality <- calculate_network_centrality(g)

  expect_true(is.data.frame(centrality))
  expect_true("gene" %in% names(centrality))
  expect_true("degree" %in% names(centrality))
  expect_true("betweenness" %in% names(centrality))
  expect_true("composite_centrality" %in% names(centrality))
})

test_that("disease_module_enrichment computes Fisher test", {
  set.seed(42)
  edges <- data.frame(
    from = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
    to = c("B", "C", "A", "E", "F", "D", "H", "I", "J", "G"),
    score = rep(800, 10),
    stringsAsFactors = FALSE
  )

  g <- build_ppi_network(edges, score_threshold = 400)
  modules <- detect_network_modules(g)

  disease_genes <- c("A", "B", "C", "X", "Y")
  enrichment <- disease_module_enrichment(modules, disease_genes)

  expect_true(is.data.frame(enrichment))
  expect_true("module" %in% names(enrichment))
  expect_true("fold_enrichment" %in% names(enrichment))
  expect_true("pvalue" %in% names(enrichment))
})
