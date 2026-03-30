# oncoAgent  

**Agentic AI-Driven Oncology Target Discovery from Multi-Omics Data** <a href="https://github.com/danymukesha/oncoAgent"><img src="man/figures/oncoAgent_logo_origin.png" align="right" height="50" alt="oncoAgent website" /></a>

![R >= 4.1.0](https://img.shields.io/badge/R-%E2%89%A54.1.0-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Tests: 135 passing](https://img.shields.io/badge/Tests-135%20passing-brightgreen)

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
┌─────────────────────────────────────────────────────────────┐
│                    Orchestration Layer                      │
│              run_oncopipeline(config)                       │
└──────────┬──────────┬─────────┬─────────┬───────────┬───────┘
           │          │         │         │           |
     ┌─────▼──┐ ┌─────▼──┐ ┌────▼───┐ ┌───▼────┐ ┌────▼─────┐
     │ Data   │ │Analy-  │ │Hypo-   │ │ Gap    │ │Evalua-   │
     │ Discov │ │sis     │ │thesis  │ │Detect  │ │tion      │
     │ Agent  │ │Agent   │ │Agent   │ │Agent   │ │Agent     │
     └────────┘ └────────┘ └────────┘ └────────┘ └──────────┘
         │           │          │          │           │
    Fetch data   Run DE/     Generate   Find gaps   Benchmark
    Check gaps   Surv/Net    mechan-    Suggest     vs COSMIC
                 pipelines   isms       fixes       CGC gold
                                                    standards
```

### Adaptive Loop

```r
# The agentic loop automatically re-runs if critical gaps are found
adaptive_result <- adaptive_pipeline_loop(
  config,
  max_iterations = 5,
  improvement_threshold = 0.01
)
```

The gap detection agent identifies:

- Missing data modalities (expression, clinical, PPI, CRISPR)
- Insufficient sample sizes (< 50 samples)
- Low DEG signal (< 100 significant genes)
- Conflicting signals between scoring dimensions

---

## Project Structure

```
oncoAgent/
├── R/
│   ├── data_ingestion.R         # TCGA, DepMap, GDSC, STRING, PubMed
│   ├── preprocessing.R          # Normalization, batch correction, QC
│   ├── differential_expression.R # DESeq2, edgeR, limma-voom, consensus
│   ├── survival_analysis.R      # Cox PH, KM, LASSO-survival
│   ├── network_analysis.R       # PPI, Louvain, centrality
│   ├── target_scoring.R         # Composite, Bayesian, druggability
│   ├── agents.R                 # 5 R6 agent classes
│   ├── orchestration.R          # Pipeline runner, adaptive loop
│   ├── evaluation.R             # AUROC, precision@K, benchmarking
│   └── utils.R                  # Config, version
├── tests/testthat/              # 135 unit tests
├── vignettes/                   # Getting started tutorial
├── _targets.R                   # Reproducible targets pipeline
├── DESCRIPTION
└── NAMESPACE
```

---

## API Reference

### Data Ingestion

| Function | Description |
|----------|-------------|
| `fetch_tcga_data()` | Download expression, clinical, mutation, methylation from TCGA |
| `fetch_depmap_crispr()` | CRISPR dependency scores from DepMap portal |
| `fetch_gdsc_drug_sensitivity()` | Drug IC50/AUC from GDSC |
| `fetch_string_ppi()` | Protein-protein interactions from STRING |
| `fetch_pubmed_abstracts()` | Literature search via PubMed E-utilities |
| `map_gene_ids()` | Cross-reference IDs via biomaRt |
| `build_reference_panel()` | Aggregate COSMIC/CGC gold standards |

### Analysis

| Function | Description |
|----------|-------------|
| `run_deseq2()` | DESeq2 differential expression |
| `run_edger()` | edgeR quasi-likelihood pipeline |
| `run_limma_voom()` | limma-voom with quality weights |
| `consensus_de()` | Multi-method consensus ranking |
| `run_cox_regression()` | Univariate Cox PH model |
| `run_kaplan_meier()` | KM curves with log-rank test |
| `build_ppi_network()` | Construct igraph PPI network |
| `detect_network_modules()` | Louvain community detection |

### Agents

| Class | Role |
|-------|------|
| `DataDiscoveryAgent` | Data availability, missing modality detection |
| `AnalysisAgent` | Execute DE/survival/network pipelines |
| `HypothesisAgent` | Generate mechanistic hypotheses |
| `GapDetectionAgent` | Find data gaps, low power, conflicts |
| `EvaluationAgent` | Benchmark vs gold standards |

### Pipeline

| Function | Description |
|----------|-------------|
| `run_oncopipeline()` | Full 7-phase pipeline execution |
| `create_pipeline_config()` | Standardized configuration |
| `adaptive_pipeline_loop()` | Gap-driven iterative re-analysis |

---

## References

[^1]: Weinstein, J.N. et al. (2013). The Cancer Genome Atlas Pan-Cancer
     analysis project. *Nature Genetics*, 45(10), 1113-1120.
     doi:[10.1038/ng.2764](https://doi.org/10.1038/ng.2764)

[^2]: Hutter, C. & Zenklusen, J.C. (2018). The Cancer Genome Atlas:
     Creating a lasting legacy. *Cell*, 173(2), 283-285.
     doi:[10.1016/j.cell.2018.03.042](https://doi.org/10.1016/j.cell.2018.03.042)

[^3]: Colaprico, A. et al. (2016). TCGAbiolinks: an R/Bioconductor package
     for integrative analysis of TCGA data. *Nucleic Acids Research*,
     44(8), e71. doi:[10.1093/nar/gkv1507](https://doi.org/10.1093/nar/gkv1507)

[^4]: Tsherniak, A. et al. (2017). Defining a cancer dependency map.
     *Cell*, 170(3), 564-576.
     doi:[10.1016/j.cell.2017.06.010](https://doi.org/10.1016/j.cell.2017.06.010)

[^5]: Vazquez, F. et al. (2024). The present and future of the Cancer
     Dependency Map. *Nature Reviews Cancer*, 24, 861-874.
     doi:[10.1038/s41568-024-00763-x](https://doi.org/10.1038/s41568-024-00763-x)

[^6]: Yang, W. et al. (2013). Genomics of Drug Sensitivity in Cancer
     (GDSC): a resource for therapeutic biomarker discovery in cancer
     cells. *Nucleic Acids Research*, 41(D1), D955-D961.
     doi:[10.1093/nar/gks1111](https://doi.org/10.1093/nar/gks1111)

[^7]: Iorio, F. et al. (2016). A landscape of pharmacogenomic interactions
     in cancer. *Cell*, 166(3), 740-754.
     doi:[10.1016/j.cell.2016.06.017](https://doi.org/10.1016/j.cell.2016.06.017)

[^8]: Szklarczyk, D. et al. (2023). STRING v12: protein-protein
     association networks with increased coverage, supporting functional
     discovery in genome-wide experimental datasets. *Nucleic Acids
     Research*, 51(D1), D605-D612.
     doi:[10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000)

[^9]: Sondka, Z. et al. (2018). The COSMIC Cancer Gene Census: describing
     genetic dysfunction across all human cancers. *Nature Reviews
     Cancer*, 18, 696-705.
     doi:[10.1038/s41568-018-0060-1](https://doi.org/10.1038/s41568-018-0060-1)

[^10]: Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation
      of fold change and dispersion for RNA-seq data with DESeq2.
      *Genome Biology*, 15, 550.
      doi:[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

[^11]: Robinson, M.D., McCarthy, D.J. & Smyth, G.K. (2010). edgeR: a
      Bioconductor package for differential expression analysis of
      digital gene expression data. *Bioinformatics*, 26(1), 139-140.
      doi:[10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

[^12]: Law, C.W. et al. (2014). voom: precision weights unlock linear
      model analysis tools for RNA-seq read counts. *Genome Biology*,
      15, R29. doi:[10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)

[^13]: Schurch, N.J. et al. (2016). How many biological replicates are
      needed in an RNA-seq experiment and which differential expression
      tool should you use? *RNA*, 22(6), 839-851.
      doi:[10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)

[^14]: Cox, D.R. (1972). Regression models and life-tables. *Journal of
      the Royal Statistical Society: Series B*, 34(2), 187-202.
      doi:[10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)

[^15]: Simon, N. et al. (2011). Regularization paths for Cox's
      proportional hazards model via coordinate descent. *Journal of
      Statistical Software*, 39(5), 1-13.
      doi:[10.18637/jss.v039.i05](https://doi.org/10.18637/jss.v039.i05)

[^16]: Blondel, V.D. et al. (2008). Fast unfolding of communities in
      large networks. *Journal of Statistical Mechanics: Theory and
      Experiment*, 2008(10), P10008.
      doi:[10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

[^17]: Guney, E., Menche, J., Vidal, M. & Barabasi, A.L. (2016).
      Network-based in silico drug efficacy screening. *Nature
      Communications*, 7, 10331.
      doi:[10.1038/ncomms10331](https://doi.org/10.1038/ncomms10331)

[^18]: Mitsos, A. et al. (2009). Identifying drug effects via pathway
      alterations using an integer linear programming optimization
      formulation on phosphoproteomic data. *PLoS Computational
      Biology*, 5(12), e1000591.
      doi:[10.1371/journal.pcbi.1000591](https://doi.org/10.1371/journal.pcbi.1000591)

---

## Citation

```bibtex
@software{oncoagent2026,
  title  = {oncoAgent: Agentic AI-Driven Oncology Target Discovery
            from Multi-Omics Data},
  author = {Dany Mukesha},
  year   = {2026},
  url    = {https://github.com/oncoagent/oncoAgent}
}
```

---

## Demo

A complete end-to-end demo is included at `inst/scripts/run_demo.R`. It:

- Builds a reference panel from COSMIC Cancer Gene Census
- Fetches PPI data from STRING (with synthetic fallback)
- Simulates realistic TCGA-BRCA-like RNA-seq data
- Runs consensus DE, survival analysis, network analysis
- Scores and ranks targets
- Benchmarks against gold standards
- Generates gap reports and hypotheses

```r
# After installing oncoAgent:
source(system.file("scripts", "run_demo.R", package = "oncoAgent"))
```

Runtime: ~2-5 minutes.

---

## License

MIT
