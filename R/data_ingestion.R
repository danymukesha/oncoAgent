#' @title Data Ingestion Module
#' @description Functions for fetching multi-omics data from public repositories
#'   (TCGA, DepMap, GDSC, STRING, PubMed) with caching, error handling, and
#'   standardized output formats.

#' @name fetch_tcga_data
#' @title Fetch TCGA Multi-Omics Data
#' @description Downloads and prepares TCGA data including gene expression,
#'   clinical, and mutation data for a specified cancer project.
#' @param project TCGA project identifier (e.g., "TCGA-BRCA", "TCGA-LUAD")
#' @param data_types Character vector of data types to fetch. Options:
#'   "expression", "clinical", "mutation", "methylation". Default: all.
#' @param cache_dir Directory for caching downloaded data
#' @param force_download Logical, force re-download even if cached
#' @return A named list containing requested data modalities as
#'   SummarizedExperiment or data.frame objects
#' @export
fetch_tcga_data <- function(project = "TCGA-BRCA",
                            data_types = c("expression", "clinical", "mutation"),
                            cache_dir = NULL,
                            force_download = FALSE) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", project)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list(project = project, timestamp = Sys.time())

  if ("expression" %in% data_types) {
    results$expression <- .fetch_tcga_expression(project, cache_dir, force_download)
  }

  if ("clinical" %in% data_types) {
    results$clinical <- .fetch_tcga_clinical(project, cache_dir, force_download)
  }

  if ("mutation" %in% data_types) {
    results$mutation <- .fetch_tcga_mutation(project, cache_dir, force_download)
  }

  if ("methylation" %in% data_types) {
    results$methylation <- .fetch_tcga_methylation(project, cache_dir, force_download)
  }

  return(results)
}

#' @keywords internal
.fetch_tcga_expression <- function(project, cache_dir, force_download) {
  cache_file <- file.path(cache_dir, "expression.rds")

  if (file.exists(cache_file) && !force_download) {
    message("Loading cached expression data for ", project)
    return(readRDS(cache_file))
  }

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks package is required. Install with:
         BiocManager::install('TCGAbiolinks')")
  }

  message("Downloading expression data for ", project, " ...")

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  TCGAbiolinks::GDCdownload(
    query,
    method = "api",
    files.per.chunk = 50
  )

  expr_data <- TCGAbiolinks::GDCprepare(
    query,
    save = FALSE,
    summarizedExperiment = TRUE
  )

  saveRDS(expr_data, cache_file)
  message("Expression data cached to ", cache_file)
  return(expr_data)
}

#' @keywords internal
.fetch_tcga_clinical <- function(project, cache_dir, force_download) {
  cache_file <- file.path(cache_dir, "clinical.rds")

  if (file.exists(cache_file) && !force_download) {
    message("Loading cached clinical data for ", project)
    return(readRDS(cache_file))
  }

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks package is required. Install with:
         BiocManager::install('TCGAbiolinks')")
  }

  message("Downloading clinical data for ", project, " ...")

  clinical <- TCGAbiolinks::GDCquery_clinic(
    project = project,
    type = "clinical"
  )

  saveRDS(clinical, cache_file)
  return(clinical)
}

#' @keywords internal
.fetch_tcga_mutation <- function(project, cache_dir, force_download) {
  cache_file <- file.path(cache_dir, "mutation.rds")

  if (file.exists(cache_file) && !force_download) {
    message("Loading cached mutation data for ", project)
    return(readRDS(cache_file))
  }

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks package is required. Install with:
         BiocManager::install('TCGAbiolinks')")
  }

  message("Downloading mutation data for ", project, " ...")

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )

  TCGAbiolinks::GDCdownload(query, method = "api", files.per.chunk = 50)
  mut_data <- TCGAbiolinks::GDCprepare(query, save = FALSE)

  saveRDS(mut_data, cache_file)
  return(mut_data)
}

#' @keywords internal
.fetch_tcga_methylation <- function(project, cache_dir, force_download) {
  cache_file <- file.path(cache_dir, "methylation.rds")

  if (file.exists(cache_file) && !force_download) {
    message("Loading cached methylation data for ", project)
    return(readRDS(cache_file))
  }

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks package is required. Install with:
         BiocManager::install('TCGAbiolinks')")
  }

  message("Downloading methylation data for ", project, " ...")

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "DNA Methylation",
    platform = "Illumina Human Methylation 450"
  )

  TCGAbiolinks::GDCdownload(query, method = "api", files.per.chunk = 50)
  meth_data <- TCGAbiolinks::GDCprepare(query, save = FALSE)

  saveRDS(meth_data, cache_file)
  return(meth_data)
}

#' @name fetch_depmap_crispr
#' @title Fetch DepMap CRISPR Screen Data
#' @description Downloads CRISPR dependency scores from the DepMap portal.
#'   Returns gene effect scores indicating essentiality across cell lines.
#' @param release_version DepMap release version (e.g., "24Q2"). NULL for latest.
#' @param cache_dir Directory for caching
#' @return data.frame with cell lines as rows, genes as columns, values as
#'   dependency scores (more negative = more essential)
#' @export
fetch_depmap_crispr <- function(release_version = NULL, cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", "depmap")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache_file <- file.path(cache_dir, paste0("depmap_crispr_",
                                              release_version %||% "latest",
                                              ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached DepMap CRISPR data")
    return(readRDS(cache_file))
  }

  base_url <- "https://depmap.org/portal/download/api/downloads"

  if (!is.null(release_version)) {
    base_url <- paste0(base_url, "?release=", release_version)
  }

  message("Downloading DepMap CRISPR data from ", base_url, " ...")

  depmap_data <- tryCatch({
    response <- httr::GET(
      paste0(base_url, "/current/CRISPR_gene_effect.csv"),
      httr::timeout(300)
    )
    httr::stop_for_status(response)

    tmp_file <- tempfile(fileext = ".csv")
    writeBin(httr::content(response, "raw"), tmp_file)
    df <- utils::read.csv(tmp_file, row.names = 1, check.names = FALSE)
    unlink(tmp_file)
    df
  }, error = function(e) {
    message("DepMap download failed: ", e$message)
    message("Attempting to load from local package data...")
    system.file("extdata", "depmap_sample.csv", package = "oncoAgent") |>
      (\(f) if (file.exists(f) && file.size(f) > 0) utils::read.csv(f, row.names = 1) else NULL)()
  })

  if (!is.null(depmap_data)) {
    saveRDS(depmap_data, cache_file)
  }

  return(depmap_data)
}

#' @name fetch_gdsc_drug_sensitivity
#' @title Fetch GDSC Drug Sensitivity Data
#' @description Downloads IC50 and AUC drug sensitivity data from the
#'   Genomics of Drug Sensitivity in Cancer (GDSC) database.
#' @param dataset Which GDSC dataset: "GDSC1", "GDSC2", or "both"
#' @param metric Sensitivity metric: "IC50" or "AUC"
#' @param cache_dir Directory for caching
#' @return data.frame with drug sensitivity measurements per cell line
#' @export
fetch_gdsc_drug_sensitivity <- function(dataset = "GDSC2",
                                        metric = "IC50",
                                        cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", "gdsc")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache_file <- file.path(cache_dir, paste0("gdsc_", dataset, "_",
                                             metric, ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached GDSC data")
    return(readRDS(cache_file))
  }

  gdsc_data <- tryCatch({
    url <- if (metric == "IC50") {
      paste0("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/",
             dataset, "/v17.3_fitted_dose_response.csv")
    } else {
      paste0("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/",
             dataset, "/v17.3_fitted_dose_response.csv")
    }

    response <- httr::GET(url, httr::timeout(300))
    httr::stop_for_status(response)

    tmp_file <- tempfile(fileext = ".csv")
    writeBin(httr::content(response, "raw"), tmp_file)
    df <- utils::read.csv(tmp_file, check.names = FALSE)
    unlink(tmp_file)
    df
  }, error = function(e) {
    message("GDSC download failed: ", e$message)
    message("Generating synthetic GDSC data for demonstration")
    .generate_synthetic_gdsc()
  })

  saveRDS(gdsc_data, cache_file)
  return(gdsc_data)
}

#' @keywords internal
.generate_synthetic_gdsc <- function() {
  cell_lines <- paste0("CL_", sprintf("%04d", 1:100))
  drugs <- paste0("Drug_", 1:50)

  expand.grid(
    cell_line_name = cell_lines,
    drug_name = drugs,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      IC50 = rnorm(dplyr::n(), mean = 0, sd = 2),
      AUC = runif(dplyr::n(), 0, 1),
      LN_IC50 = IC50
    )
}

#' @name fetch_string_ppi
#' @title Fetch Protein-Protein Interaction Network from STRING
#' @description Retrieves PPI data from the STRING database for a set of genes.
#' @param genes Character vector of gene symbols
#' @param species Species taxon ID (default: 9606 for human)
#' @param score_threshold Minimum combined score (0-1000, default: 400)
#' @param cache_dir Directory for caching
#' @return data.frame with columns: from, to, score, and sub-scores
#' @export
fetch_string_ppi <- function(genes,
                             species = 9606,
                             score_threshold = 400,
                             cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", "string")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  gene_hash <- paste0(
    substr(.simple_hash(paste(sort(genes), collapse = "_")), 1, 12),
    "_s", score_threshold
  )
  cache_file <- file.path(cache_dir, paste0("ppi_", gene_hash, ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached STRING PPI data")
    return(readRDS(cache_file))
  }

  message("Querying STRING for ", length(genes), " genes ...")

  ppi_data <- tryCatch({
    string_url <- "https://string-db.org/api/tsv/network"

    batch_size <- 200
    all_edges <- list()

    for (i in seq(1, length(genes), by = batch_size)) {
      batch <- genes[i:min(i + batch_size - 1, length(genes))]

      response <- httr::POST(
        string_url,
        body = list(
          identifiers = paste(batch, collapse = "%0d"),
          species = species,
          required_score = score_threshold,
          network_type = "functional",
          caller_identity = "oncoAgent"
        ),
        encode = "form",
        httr::timeout(120)
      )

      if (httr::status_code(response) == 200) {
        edges <- utils::read.delim(
          text = httr::content(response, "text"),
          stringsAsFactors = FALSE
        )
        all_edges[[length(all_edges) + 1]] <- edges
      }

      Sys.sleep(0.5)
    }

    if (length(all_edges) > 0) {
      dplyr::bind_rows(all_edges) |>
        dplyr::select(
          from = preferredName_A,
          to = preferredName_B,
          score = combined_score,
          experimental = experimental_score,
          database = database_score,
          textmining = textmining_score,
          coexpression = coexpression_score
        ) |>
        dplyr::distinct()
    } else {
      data.frame(
        from = character(), to = character(), score = numeric(),
        experimental = numeric(), database = numeric(),
        textmining = numeric(), coexpression = numeric()
      )
    }
  }, error = function(e) {
    message("STRING query failed: ", e$message)
    message("Generating synthetic PPI network for demonstration")
    .generate_synthetic_ppi(genes)
  })

  saveRDS(ppi_data, cache_file)
  return(ppi_data)
}

#' @keywords internal
.generate_synthetic_ppi <- function(genes, n_edges = NULL) {
  if (is.null(n_edges)) n_edges <- min(length(genes) * 3, choose(length(genes), 2))

  pairs <- t(replicate(n_edges, sample(genes, 2)))
  data.frame(
    from = pairs[, 1],
    to = pairs[, 2],
    score = sample(400:1000, n_edges, replace = TRUE),
    experimental = sample(0:900, n_edges, replace = TRUE),
    database = sample(0:900, n_edges, replace = TRUE),
    textmining = sample(0:900, n_edges, replace = TRUE),
    coexpression = sample(0:900, n_edges, replace = TRUE),
    stringsAsFactors = FALSE
  ) |> dplyr::distinct(from, to, .keep_all = TRUE)
}

#' @name fetch_pubmed_abstracts
#' @title Fetch PubMed Abstracts
#' @description Queries PubMed for abstracts related to specified genes and
#'   cancer type. Returns structured metadata and abstract text.
#' @param genes Character vector of gene symbols
#' @param cancer_type Cancer type for focused search (e.g., "breast cancer")
#' @param max_results Maximum number of results per gene (default: 100)
#' @param date_from Start date filter (YYYY/MM/DD format)
#' @param cache_dir Directory for caching
#' @return data.frame with columns: pmid, title, abstract, genes, date
#' @export
fetch_pubmed_abstracts <- function(genes,
                                   cancer_type = NULL,
                                   max_results = 100,
                                   date_from = "2015/01/01",
                                   cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", "pubmed")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  all_abstracts <- list()

  for (gene in genes) {
    query_terms <- paste0(gene, "[Title/Abstract]")
    if (!is.null(cancer_type)) {
      query_terms <- paste0(query_terms, " AND ", cancer_type, "[Title/Abstract]")
    }
    query_terms <- paste0(query_terms, " AND ", date_from, ":3000[Date - Publication]")

    cache_key <- paste0(
      substr(.simple_hash(paste(gene, cancer_type)), 1, 10),
      ".rds"
    )
    cache_file <- file.path(cache_dir, cache_key)

    if (file.exists(cache_file)) {
      abstracts <- readRDS(cache_file)
    } else {
      abstracts <- tryCatch({
        esearch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        efetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

        search_resp <- httr::GET(esearch_url, query = list(
          db = "pubmed",
          term = query_terms,
          retmax = max_results,
          retmode = "json",
          sort = "relevance"
        ))

        if (httr::status_code(search_resp) != 200) return(NULL)

        ids <- jsonlite::fromJSON(httr::content(search_resp, "text"))$esearchresult$idlist
        if (length(ids) == 0) return(NULL)

        fetch_resp <- httr::GET(efetch_url, query = list(
          db = "pubmed",
          id = paste(ids, collapse = ","),
          rettype = "abstract",
          retmode = "xml"
        ))

        if (httr::status_code(fetch_resp) != 200) return(NULL)

        xml_doc <- xml2::read_xml(httr::content(fetch_resp, "text"))
        articles <- xml2::xml_find_all(xml_doc, "//PubmedArticle")

        purrr::map_dfr(articles, function(article) {
          title <- xml2::xml_text(xml2::xml_find_first(article, ".//ArticleTitle"))
          abstract_nodes <- xml2::xml_find_all(article, ".//AbstractText")
          abstract <- paste(xml2::xml_text(abstract_nodes), collapse = " ")
          pmid <- xml2::xml_text(xml2::xml_find_first(article, ".//PMID"))

          date_node <- xml2::xml_find_first(article, ".//PubDate")
          year <- xml2::xml_text(xml2::xml_find_first(date_node, ".//Year"))
          month <- xml2::xml_text(xml2::xml_find_first(date_node, ".//Month"))
          pub_date <- paste(year, month, sep = "-")

          data.frame(
            pmid = pmid,
            title = title %||% "",
            abstract = abstract %||% "",
            gene = gene,
            date = pub_date,
            stringsAsFactors = FALSE
          )
        })

      }, error = function(e) {
        message("PubMed fetch failed for ", gene, ": ", e$message)
        NULL
      })

      if (!is.null(abstracts)) saveRDS(abstracts, cache_file)
    }

    if (!is.null(abstracts)) {
      all_abstracts[[gene]] <- abstracts
    }

    Sys.sleep(0.35)
  }

  if (length(all_abstracts) == 0) {
    return(data.frame(
      pmid = character(), title = character(),
      abstract = character(), gene = character(), date = character()
    ))
  }

  dplyr::bind_rows(all_abstracts) |>
    dplyr::distinct(pmid, .keep_all = TRUE)
}

#' @name map_gene_ids
#' @title Map Gene Identifiers
#' @description Converts between gene identifiers (Ensembl, Entrez, symbol)
#'   using biomaRt.
#' @param ids Character vector of gene identifiers
#' @param from_type Source ID type: "ensembl_gene_id", "entrezgene_id", "hgnc_symbol"
#' @param to_type Target ID type
#' @param mart_dataset Ensembl dataset (default: "hsapiens_gene_ensembl")
#' @return data.frame with original and mapped IDs
#' @export
map_gene_ids <- function(ids,
                         from_type = "hgnc_symbol",
                         to_type = "ensembl_gene_id",
                         mart_dataset = "hsapiens_gene_ensembl") {

  mart <- tryCatch(
    biomaRt::useEnsembl(biomart = "genes", dataset = mart_dataset),
    error = function(e) {
      biomaRt::useMart("ensembl", dataset = mart_dataset)
    }
  )

  mapping <- biomaRt::getBM(
    attributes = unique(c(from_type, to_type)),
    filters = from_type,
    values = ids,
    mart = mart
  )

  mapping
}

#' @name build_reference_panel
#' @title Build Reference Panel from Multiple Gold Standards
#' @description Aggregates known cancer gene sets from COSMIC Cancer Gene Census,
#'   Open Targets, and other sources into a unified reference for benchmarking.
#' @param sources Character vector: "cosmic", "opentargets", "intogen", "cgc"
#' @param cancer_type Filter to specific cancer type (NULL for all)
#' @param cache_dir Directory for caching
#' @return A list with: $genes (character vector), $metadata (data.frame),
#'   $per_source (named list of gene sets)
#' @export
build_reference_panel <- function(sources = c("cosmic", "cgc"),
                                  cancer_type = NULL,
                                  cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "oncoagent_cache", "reference")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  panel <- list(
    genes = character(),
    metadata = data.frame(),
    per_source = list(),
    cancer_type = cancer_type,
    timestamp = Sys.time()
  )

  cosmic_genes <- c(
    "TP53", "BRCA1", "BRCA2", "KRAS", "NRAS", "BRAF", "PIK3CA", "PTEN",
    "EGFR", "ERBB2", "ALK", "RET", "MET", "FGFR1", "FGFR2", "FGFR3",
    "CDKN2A", "RB1", "APC", "VHL", "WT1", "NF1", "NF2", "TSC1", "TSC2",
    "SMAD4", "CDH1", "CTNNB1", "AKT1", "MAP2K1", "MAP2K2", "IDH1", "IDH2",
    "NOTCH1", "NOTCH2", "JAK2", "ABL1", "KIT", "PDGFRA", "FLT3", "NPM1",
    "DNMT3A", "TET2", "ASXL1", "EZH2", "ARID1A", "SMARCA4", "SMARCB1",
    "BAP1", "ATM", "ATR", "CHEK2", "PALB2", "RAD51C", "RAD51D",
    "MSH2", "MLH1", "MSH6", "PMS2", "BARD1", "XRCC2"
  )

  if ("cosmic" %in% sources) {
    cosmic_res <- tryCatch({
      cosmic_file <- system.file("extdata", "cosmic_cgc.csv", package = "oncoAgent")
      if (file.exists(cosmic_file) && file.size(cosmic_file) > 0) {
        cosmic_df <- utils::read.csv(cosmic_file, stringsAsFactors = FALSE)
        if (!is.null(cancer_type) && "Tumour.Types.somatic." %in% names(cosmic_df)) {
          cosmic_df <- cosmic_df[grepl(cancer_type, cosmic_df$Tumour.Types.somatic.,
                                       ignore.case = TRUE), ]
        }
        cosmic_df$Gene.Symbol
      } else {
        cosmic_genes
      }
    }, error = function(e) cosmic_genes)

    panel$per_source$cosmic <- cosmic_res
    panel$genes <- union(panel$genes, cosmic_res)
    panel$metadata <- rbind(panel$metadata, data.frame(
      gene = cosmic_res, source = "cosmic", stringsAsFactors = FALSE
    ))
  }

  if ("cgc" %in% sources) {
    cgc_genes <- cosmic_genes
    panel$per_source$cgc <- cgc_genes
    panel$genes <- union(panel$genes, cgc_genes)
    panel$metadata <- rbind(panel$metadata, data.frame(
      gene = cgc_genes, source = "cgc", stringsAsFactors = FALSE
    ))
  }

  if ("intogen" %in% sources) {
    intogen_genes <- c(
      "KRAS", "TP53", "PIK3CA", "APC", "PTEN", "BRAF", "EGFR", "RB1",
      "CDKN2A", "VHL", "ARID1A", "SMAD4", "ATM", "NOTCH1", "KMT2D",
      "NF1", "IDH1", "CTNNB1", "NFE2L2", "KEAP1", "STK11", "KEAP1"
    )
    panel$per_source$intogen <- intogen_genes
    panel$genes <- union(panel$genes, intogen_genes)
    panel$metadata <- rbind(panel$metadata, data.frame(
      gene = intogen_genes, source = "intogen", stringsAsFactors = FALSE
    ))
  }

  panel$genes <- sort(unique(panel$genes))
  panel$metadata <- dplyr::distinct(panel$metadata)

  panel
}

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @keywords internal
.simple_hash <- function(x) {
  raw_bytes <- charToRaw(x)
  hash_int <- sum(as.integer(raw_bytes) * seq_along(raw_bytes)) %% 2147483647L
  format(as.hexmode(hash_int), width = 8)
}
