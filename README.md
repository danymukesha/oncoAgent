# oncoAgent  

**Agentic AI-Driven Oncology Target Discovery from Multi-Omics Data** <a href="https://github.com/danymukesha/oncoAgent"><img src="man/figures/oncoAgent_logo_origin.png" align="right" height="120" alt="oncoAgent website" /></a>

![R >= 4.1.0](https://img.shields.io/badge/R-%E2%89%A54.1.0-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Tests: 135 passing](https://img.shields.io/badge/Tests-135%20passing-brightgreen)
[![R-CMD-check](https://github.com/danymukesha/oncoAgent/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danymukesha/oncoAgent/actions/workflows/R-CMD-check.yaml)

---

## Overview

**oncoAgent** is a semi-autonomous research system that integrates multi-omics
data (genomics, transcriptomics, CRISPR dependency screens, drug sensitivity)
with protein interaction networks and clinical outcomes to identify and
prioritize oncology drug targets.

The system features **five specialized AI agents** that collaborate through
adaptive workflows: detecting data gaps, suggesting next analyses, generating
hypotheses, and benchmarking predictions against gold-standard cancer gene sets
(COSMIC, Cancer Gene Census). Unlike static pipelines, oncoAgent iteratively
refines its analysis until evidence quality meets predefined thresholds.

### Why oncoAgent?

| Feature | Traditional Pipeline | oncoAgent |
|---------|---------------------|-----------|
| Data handling | Manual download + preprocessing | Automated ingestion with caching |
| DE analysis | Single method | Consensus across DESeq2, edgeR, limma-voom |
| Gap detection | None | Automated detection of missing data / low power |
| Target scoring | Ad hoc | Multi-dimensional weighted composite |
| Iteration | Fixed pipeline | Adaptive re-analysis driven by gap agent |
| Benchmarking | Manual | Automated AUROC, precision@K vs gold standards |

---

## Scientific Foundation

oncoAgent builds on established methods from the following peer-reviewed work:

### Data Sources

- **The Cancer Genome Atlas (TCGA)**: Pan-cancer genomic characterization
  resource spanning 33 cancer types, 11,000+ tumors [^1][^2]. oncoAgent
  uses TCGAbiolinks for programmatic access [^3].

- **Cancer Dependency Map (DepMap)**: Genome-scale CRISPR-Cas9 essentiality
  screens across 1,000+ cancer cell lines [^4][^5]. Used to identify
  genes with selective dependencies in specific cancer contexts.

- **Genomics of Drug Sensitivity in Cancer (GDSC)**: Systematic drug response
  profiling linking genomic features to pharmacological sensitivity [^6][^7].

- **STRING Database**: Protein-protein interaction network covering >14,000
  organisms with scored evidence from experiments, databases, and text
  mining [^8].

- **COSMIC Cancer Gene Census**: Curated list of genes with mutations causally
  implicated in cancer, used as the primary gold standard for benchmarking [^9].

### Analysis Methods

- **Differential Expression**: Consensus ranking across DESeq2 [^10],
  edgeR [^11], and limma-voom [^12], following best practices from
  systematic benchmarking studies [^13].

- **Survival Analysis**: Cox proportional hazards regression [^14] with
  LASSO-penalized variable selection for multivariate models [^15].

- **Network Analysis**: Louvain community detection on PPI networks [^16]
  with centrality-based target prioritization, following the network
  proximity paradigm [^17].

- **Target Scoring**: Composite scoring integrating differential expression
  significance, survival association, network centrality, and CRISPR
  dependency, with Bayesian posterior ranking [^18].

---

## Installation

```r
# Install Bioconductor dependencies first
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "SummarizedExperiment",
                        "biomaRt", "pROC"))

# Install CRAN dependencies
install.packages(c("survival", "glmnet", "igraph", "data.table",
                    "dplyr", "tidyr", "purrr", "stringr", "jsonlite",
                    "httr", "xml2", "ggplot2", "R6", "parallel"))

# Install oncoAgent from local source
# R CMD INSTALL oncoAgent
```

---

## Quick Start

### 1. Load the package and configure

```r
library(oncoAgent)

config <- create_pipeline_config(
  project = "TCGA-BRCA",
  cancer_type = "breast cancer",
  de_methods = c("deseq2", "edger", "limma"),
  padj_threshold = 0.05,
  output_dir = "results"
)
```

### 2. Run with synthetic demo data

```r
set.seed(42)

# Simulate RNA-seq count matrix (1000 genes x 40 samples)
n_genes <- 1000
n_samples <- 40

counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 200, size = 10),
  nrow = n_genes, ncol = n_samples
)
rownames(counts) <- paste0("Gene_", seq_len(n_genes))
colnames(counts) <- paste0("Sample_", seq_len(n_samples))

# Inject signal in first 50 genes (upregulated in tumor)
counts[1:50, 1:20] <- rnbinom(50 * 20, mu = 1000, size = 10)

# Sample metadata
col_data <- data.frame(
  condition = factor(rep(c("tumor", "normal"), each = 20)),
  row.names = colnames(counts)
)

# Clinical data
clinical <- data.frame(
  sample_id = colnames(counts),
  time = rexp(n_samples, rate = 0.002),
  status = rbinom(n_samples, 1, 0.4),
  stringsAsFactors = FALSE
)

# Run full pipeline
result <- run_oncopipeline(
  config,
  existing_data = list(
    expression = counts,
    col_data = col_data,
    clinical = clinical
  )
)

# View top targets
head(result$targets[, c("gene", "composite_score", "target_tier")], 10)
```

### 3. Run individual modules

#### Differential Expression (Consensus)

```r
de <- consensus_de(counts, col_data,
                   methods = c("deseq2", "edger", "limma"),
                   padj_threshold = 0.05)

# Top DE genes detected by >= 2 methods
head(de[de$n_methods >= 2, ], 10)
```

#### Survival Analysis

```r
# Univariate Cox regression per gene
hr_table <- calculate_hazard_ratios(
  counts,
  time = clinical$time,
  status = clinical$status,
  genes = rownames(counts)[1:100]
)

# Kaplan-Meier for a specific gene
km <- run_kaplan_meier(
  counts["Gene_1", ], clinical$time, clinical$status,
  gene_name = "Gene_1", split_method = "median"
)
```

#### Network Analysis

```r
# Build PPI from synthetic edges (in practice, fetch from STRING)
ppi <- data.frame(
  from = sample(rownames(counts)[1:100], 200, replace = TRUE),
  to = sample(rownames(counts)[1:100], 200, replace = TRUE),
  score = sample(400:900, 200, replace = TRUE),
  stringsAsFactors = FALSE
)

g <- build_ppi_network(ppi, score_threshold = 400)
centrality <- calculate_network_centrality(g)
modules <- detect_network_modules(g)
```

#### Target Scoring & Ranking

```r
scores <- composite_target_score(
  de_results = de,
  survival_hr = hr_table,
  centrality = centrality
)

ranked <- rank_targets(scores, top_n = 20)
print(ranked[, c("gene", "composite_score", "target_tier")])
```

#### Evaluation Against Gold Standards

```r
gs <- build_reference_panel(sources = c("cosmic", "cgc"))

eval_agent <- EvaluationAgent$new()
metrics <- eval_agent$evaluate(ranked, gs)
cat(eval_agent$report())
```

---

## Agentic Architecture

oncoAgent implements five R6 agent classes that collaborate in an adaptive loop:

```
bb