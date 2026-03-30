#' @title Network Analysis Module
#' @description Protein-protein interaction network construction, community
#'   detection, centrality analysis, and disease module identification.

#' @name build_ppi_network
#' @title Build PPI Network
#' @description Constructs an igraph network from STRING PPI data or a custom
#'   edge list. Filters by score and removes disconnected components.
#' @param ppi_data data.frame with columns: from, to, score
#' @param score_threshold Minimum interaction score (default: 400)
#' @param remove_isolates Logical, remove unconnected nodes (default: TRUE)
#' @param directed Logical, create directed graph (default: FALSE)
#' @return igraph object
#' @export
build_ppi_network <- function(ppi_data,
                              score_threshold = 400,
                              remove_isolates = TRUE,
                              directed = FALSE) {

  filtered <- ppi_data[ppi_data$score >= score_threshold, ]

  if (nrow(filtered) == 0) {
    warning("No edges above score threshold ", score_threshold)
    return(igraph::make_empty_graph(n = 0, directed = directed))
  }

  g <- igraph::graph_from_data_frame(
    d = filtered[, c("from", "to", "score")],
    directed = directed
  )

  if (remove_isolates) {
    isolates <- which(igraph::degree(g) == 0)
    if (length(isolates) > 0) {
      g <- igraph::delete_vertices(g, isolates)
    }
  }

  igraph::E(g)$weight <- igraph::E(g)$score / 1000

  message("PPI network: ", igraph::vcount(g), " nodes, ",
          igraph::ecount(g), " edges",
          if (remove_isolates) " (isolates removed)")

  g
}

#' @name detect_network_modules
#' @title Detect Network Modules (Communities)
#' @description Identifies functional modules in the PPI network using
#'   Louvain community detection.
#' @param g igraph network object
#' @param method Community detection: "louvain", "walktrap", "infomap"
#' @param resolution Resolution parameter for Louvain (default: 1.0)
#' @return List with: membership, modularity, module_genes, module_sizes
#' @export
detect_network_modules <- function(g,
                                   method = "louvain",
                                   resolution = 1.0) {

  if (igraph::vcount(g) < 3) {
    warning("Network too small for community detection")
    return(list(membership = rep(1, igraph::vcount(g)),
                modularity = 0,
                module_genes = list("1" = igraph::V(g)$name)))
  }

  if (method == "louvain") {
    comm <- igraph::cluster_louvain(g, resolution = resolution)
  } else if (method == "walktrap") {
    comm <- igraph::cluster_walktrap(g)
  } else if (method == "infomap") {
    comm <- igraph::cluster_infomap(g)
  } else {
    stop("Unknown community detection method: ", method)
  }

  membership <- igraph::membership(comm)
  modularity <- igraph::modularity(comm)

  module_genes <- split(igraph::V(g)$name, membership)

  result <- list(
    membership = membership,
    modularity = modularity,
    algorithm = method,
    module_genes = module_genes,
    module_sizes = sapply(module_genes, length),
    n_modules = length(module_genes),
    object = comm
  )

  message("Detected ", result$n_modules, " modules (modularity = ",
          round(modularity, 3), ")")

  result
}

#' @name calculate_network_centrality
#' @title Calculate Network Centrality Metrics
#' @description Computes degree, betweenness, closeness, eigenvector, and
#'   PageRank centrality for all nodes.
#' @param g igraph network object
#' @param metrics Character vector of metrics to compute
#' @return data.frame with centrality scores per node
#' @export
calculate_network_centrality <- function(g,
                                         metrics = c("degree", "betweenness",
                                                     "closeness", "pagerank")) {

  nodes <- igraph::V(g)$name
  centrality <- data.frame(gene = nodes, stringsAsFactors = FALSE)

  if ("degree" %in% metrics) {
    centrality$degree <- igraph::degree(g)
    centrality$degree_norm <- centrality$degree / (igraph::vcount(g) - 1)
  }

  if ("betweenness" %in% metrics) {
    centrality$betweenness <- igraph::betweenness(g, normalized = TRUE)
  }

  if ("closeness" %in% metrics) {
    centrality$closeness <- igraph::closeness(g, normalized = TRUE)
  }

  if ("eigenvector" %in% metrics) {
    centrality$eigenvector <- tryCatch(
      igraph::eigen_centrality(g)$vector,
      error = function(e) rep(0, igraph::vcount(g))
    )
  }

  if ("pagerank" %in% metrics) {
    pr <- igraph::page_rank(g)
    centrality$pagerank <- pr$vector
  }

  if (length(metrics) > 1) {
    numeric_cols <- sapply(centrality, is.numeric)
    centrality$composite_centrality <- rowMeans(
      scale(centrality[, numeric_cols]),
      na.rm = TRUE
    )
  }

  centrality <- centrality[order(-centrality$composite_centrality), ]
  centrality
}

#' @name disease_module_enrichment
#' @title Disease Module Enrichment
#' @description Tests whether network modules are enriched for known disease
#'   genes using Fisher's exact test.
#' @param modules Module detection output from detect_network_modules()
#' @param disease_genes Character vector of known disease gene symbols
#' @param universe Character vector of all genes considered (background)
#' @return data.frame with enrichment results per module
#' @export
disease_module_enrichment <- function(modules,
                                      disease_genes,
                                      universe = NULL) {

  if (is.null(universe)) {
    universe <- unique(unlist(modules$module_genes))
  }

  enrichment <- purrr::map_dfr(names(modules$module_genes), function(mod_id) {
    module_genes <- modules$module_genes[[mod_id]]

    in_module <- universe %in% module_genes
    in_disease <- universe %in% disease_genes

    contingency <- table(in_module, in_disease)

    if (nrow(contingency) < 2 || ncol(contingency) < 2) {
      return(data.frame(
        module = mod_id,
        n_genes = length(module_genes),
        n_disease_genes = 0,
        fold_enrichment = 0,
        odds_ratio = NA,
        pvalue = 1,
        stringsAsFactors = FALSE
      ))
    }

    ft <- stats::fisher.test(contingency, alternative = "greater")

    n_overlap <- sum(module_genes %in% disease_genes)
    expected <- length(module_genes) * length(disease_genes) / length(universe)
    fold_enrich <- if (expected > 0) n_overlap / expected else 0

    data.frame(
      module = mod_id,
      n_genes = length(module_genes),
      n_disease_genes = n_overlap,
      overlap_genes = paste(intersect(module_genes, disease_genes), collapse = ";"),
      expected = round(expected, 2),
      fold_enrichment = round(fold_enrich, 2),
      odds_ratio = round(ft$estimate, 2),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  enrichment$padj <- stats::p.adjust(enrichment$pvalue, method = "BH")
  enrichment <- enrichment[order(enrichment$pvalue), ]

  enrichment
}

#' @name export_network
#' @title Export Network for Visualization
#' @description Exports the network in formats compatible with Cytoscape
#'   or other network visualization tools.
#' @param g igraph network object
#' @param output_file Output file path (.graphml, .sif, or .csv)
#' @param node_attributes data.frame of node attributes to include
#' @param format Export format: "graphml", "sif", "edgelist"
#' @return Invisibly returns the file path
#' @export
export_network <- function(g,
                           output_file,
                           node_attributes = NULL,
                           format = "graphml") {

  if (!is.null(node_attributes)) {
    for (col in names(node_attributes)) {
      if (col == "gene") next
      matched <- match(igraph::V(g)$name, node_attributes$gene)
      igraph::vertex_attr(g, col) <- node_attributes[[col]][matched]
    }
  }

  if (format == "graphml") {
    igraph::write_graph(g, output_file, format = "graphml")
  } else if (format == "sif") {
    edges <- igraph::as_data_frame(g, what = "edges")
    sif <- data.frame(
      source = edges$from,
      interaction = "pp",
      target = edges$to,
      stringsAsFactors = FALSE
    )
    utils::write.table(sif, output_file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
  } else if (format == "edgelist") {
    edges <- igraph::as_data_frame(g, what = "edges")
    utils::write.csv(edges, output_file, row.names = FALSE)
  } else {
    stop("Unknown format: ", format)
  }

  message("Network exported to ", output_file)
  invisible(output_file)
}
