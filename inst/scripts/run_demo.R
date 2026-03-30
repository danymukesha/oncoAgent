#!/usr/bin/env Rscript
# ============================================================================
# oncoAgent Demo: End-to-End Oncology Target Discovery
# ============================================================================
#
# This script demonstrates the full oncoAgent pipeline:
#   1. Build reference panel from COSMIC Cancer Gene Census
#   2. Fetch real PPI data from STRING (for known breast cancer genes)
#   3. Simulate realistic TCGA-BRCA-like RNA-seq data
#   4. Run consensus differential expression
#   5. Survival analysis
#   6. Network-based target prioritization
#   7. Composite scoring and ranking
#   8. Evaluation against gold standards
#   9. Gap detection and adaptive recommendations
#
# Requirements: oncoAgent installed (R CMD INSTALL oncoAgent)
# Runtime: ~2-5 minutes depending on network speed
# ============================================================================

cat("================================================================\n")
cat("  oncoAgent Demo: Oncology Target Discovery Pipeline\n")
cat("  Date:", as.character(Sys.Date()), "\n")
cat("================================================================\n\n")

# --- Load oncoAgent ----------------------------------------------------------
library(oncoAgent)

# --- Step 1: Build Reference Panel from COSMIC Cancer Gene Census ------------
cat("--- Step 1: Build Reference Panel ---\n")
#
# The COSMIC Cancer Gene Census (Sondka et al., Nat Rev Cancer 2018) is the
# most widely used gold standard for cancer gene identification. We use it
# to benchmark our target predictions.
#
# Reference: Sondka, Z. et al. (2018). The COSMIC Cancer Gene Census:
# describing genetic dysfunction across all human cancers. Nature Reviews
# Cancer, 18, 696-705. doi:10.1038/s41568-018-0060-1

# Curated set of well-established breast cancer driver genes
# Sources: COSMIC CGC, TCGA BRCA (Cancer Genome Atlas Network, Nature 2012),
# and Open Targets (Ochoa et al., Nucleic Acids Res 2022)
known_breast_cancer_genes <- c(
  # Tier 1: Well-established drivers
  "TP53", "PIK3CA", "ERBB2", "GATA3", "CDH1", "MAP3K1", "MAP2K4",
  "AKT1", "FOXA1", "SF3B1", "CBFB", "RUNX1", "PTEN", "RB1", "BRCA1",
  "BRCA2", "ATM", "PALB2", "CHEK2", "CDK4", "CDK6", "CCND1",
  # Tier 2: Strong evidence
  "MYC", "ESR1", "ARID1A", "NF1", "MED12", "KRAS", "NF2",
  "NOTCH1", "NOTCH2", "TBX3", "NCOR1", "SMARCD1", "PIK3R1",
  # Tier 3: Emerging targets
  "MDM2", "EGFR", "FGFR1", "FGFR2", "IGF1R", "MET", "ALK",
  "MTOR", "TSC1", "TSC2", "CTNNB1", "APC", "SMAD4", "BAP1",
  "IDH1", "IDH2", "EZH2", "KMT2C", "ARID1B", "KDM6A"
)

# Build the reference panel using COSMIC and CGC sources
gs <- build_reference_panel(sources = c("cosmic", "cgc"))
cat("  Gold standard genes:", length(gs$genes), "\n")
cat("  COSMIC genes:", length(gs$per_source$cosmic), "\n")
cat("  CGC genes:", length(gs$per_source$cgc), "\n\n")

# --- Step 2: Fetch Real PPI Data from STRING ---------------------------------
cat("--- Step 2: Fetch PPI Network from STRING ---\n")
#
# STRING database (Szklarczyk et al., Nucleic Acids Res 2023) provides
# protein-protein interaction evidence from experiments, databases, and
# text mining. We query for our known breast cancer genes.
#
# Reference: Szklarczyk, D. et al. (2023). STRING v12: protein-protein
# association networks with increased coverage. Nucleic Acids Research,
# 51(D1), D605-D612. doi:10.1093/nar/gkac1000

# Query STRING for interactions among breast cancer genes
query_genes <- known_breast_cancer_genes[1:30]

cat("  Querying STRING for", length(query_genes), "genes...\n")
ppi_data <- tryCatch(
  fetch_string_ppi(query_genes, score_threshold = 400),
  error = function(e) {
    cat("  STRING query failed:", e$message, "\n")
    cat("  Generating representative PPI network...\n")
    # Generate a biologically-plausible synthetic network
    # based on known interaction patterns
    edges <- data.frame(
      from = c(
        "TP53", "TP53", "TP53", "BRCA1", "BRCA2", "ATM", "ATM",
        "ERBB2", "PIK3CA", "PIK3CA", "AKT1", "PTEN", "PTEN",
        "CDH1", "GATA3", "FOXA1", "ESR1", "ESR1", "CCND1",
        "CDK4", "CDK6", "RB1", "RB1", "MYC", "MYC",
        "EGFR", "EGFR", "KRAS", "KRAS", "MAP3K1",
        "BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
        "NOTCH1", "NOTCH2", "NF1", "PTEN", "PIK3R1"
      ),
      to = c(
        "BRCA1", "MDM2", "RB1", "BRCA2", "PALB2", "CHEK2", "BRCA1",
        "PIK3CA", "AKT1", "PTEN", "MTOR", "PIK3CA", "AKT1",
        "CTNNB1", "FOXA1", "ESR1", "GATA3", "FOXA1", "CDK4",
        "CDK6", "RB1", "TP53", "EZH2", "TP53", "CDK4",
        "ERBB2", "PIK3CA", "BRAF", "MAP3K1", "MAP2K4",
        "PALB2", "CHEK2", "BRCA1", "TP53", "BRCA1",
        "NOTCH2", "NOTCH1", "KRAS", "PIK3CA", "PIK3CA"
      ),
      score = sample(400:999, 40, replace = TRUE),
      experimental = sample(0:500, 40, replace = TRUE),
      database = sample(0:500, 40, replace = TRUE),
      textmining = sample(0:500, 40, replace = TRUE),
      coexpression = sample(0:500, 40, replace = TRUE),
      stringsAsFactors = FALSE
    )
    # Remove self-loops
    edges <- edges[edges$from != edges$to, ]
    edges
  }
)

cat("  PPI edges retrieved:", nrow(ppi_data), "\n")

# Build the network
g <- build_ppi_network(ppi_data, score_threshold = 400)
cat("  Network: ", igraph::vcount(g), "nodes,", igraph::ecount(g), "edges\n")

# Compute centrality
centrality <- calculate_network_centrality(g)
cat("  Top hub genes:", paste(head(centrality$gene, 5), collapse = ", "), "\n\n")

# --- Step 3: Simulate Realistic TCGA-BRCA-like Data -------------------------
cat("--- Step 3: Simulate TCGA-BRCA-like Expression Data ---\n")
#
# We simulate RNA-seq count data that mimics the statistical properties
# of real TCGA breast cancer data:
#   - 5000 genes, 60 samples (30 tumor, 30 normal)
#   - ~200 truly differentially expressed genes
#   - Negative binomial distribution (standard for RNA-seq)
#   - Log-fold changes in the range 0.5-3.0 for DE genes
#
# Reference: Schurch, N.J. et al. (2016). How many biological replicates
# are needed in an RNA-seq experiment? RNA, 22(6), 839-851.
# doi:10.1261/rna.053959.115

set.seed(2026)

n_genes <- 5000
n_samples <- 60
n_tumor <- 30
n_normal <- 30

# Base expression: negative binomial with mean ~200, dispersion ~0.1
base_counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 200, size = 10),
  nrow = n_genes, ncol = n_samples
)
rownames(base_counts) <- paste0("Gene_", seq_len(n_genes))
colnames(base_counts) <- c(
  paste0("TCGA-Tumor-", sprintf("%03d", seq_len(n_tumor))),
  paste0("TCGA-Normal-", sprintf("%03d", seq_len(n_normal)))
)

# Inject DE signal in first 200 genes (upregulated in tumor)
# Log-fold changes range from 0.5 to 3.0
de_genes_idx <- 1:200
lfc <- runif(200, 0.5, 3.0)
fold_change <- 2^lfc

for (i in de_genes_idx) {
  base_counts[i, 1:n_tumor] <- rnbinom(
    n_tumor,
    mu = 200 * fold_change[i],
    size = 10
  )
}

# Rename first 20 DE genes to known cancer genes for realistic evaluation
cancer_names <- c(
  "TP53", "PIK3CA", "ERBB2", "GATA3", "CDH1", "MAP3K1", "AKT1",
  "FOXA1", "PTEN", "RB1", "BRCA1", "BRCA2", "ATM", "CCND1",
  "CDK4", "CDK6", "MYC", "ESR1", "KRAS", "EGFR"
)
rownames(base_counts)[1:20] <- cancer_names

# Sample metadata
col_data <- data.frame(
  condition = factor(rep(c("tumor", "normal"), c(n_tumor, n_normal))),
  row.names = colnames(base_counts)
)

# Clinical data with survival endpoints
# Tumor patients have worse prognosis (realistic for breast cancer)
clinical <- data.frame(
  sample_id = colnames(base_counts),
  time = c(
    rexp(n_tumor, rate = 0.003),   # Tumor: median ~330 days
    rexp(n_normal, rate = 0.001)   # Normal: median ~1000 days
  ),
  status = c(
    rbinom(n_tumor, 1, 0.5),       # Tumor: 50% events
    rbinom(n_normal, 1, 0.2)       # Normal: 20% events
  ),
  age = round(rnorm(n_samples, 55, 12)),
  stage = sample(c("I", "II", "III", "IV"), n_samples,
                 replace = TRUE, prob = c(0.15, 0.40, 0.30, 0.15)),
  stringsAsFactors = FALSE
)

cat("  Expression matrix:", n_genes, "genes x", n_samples, "samples\n")
cat("  Tumor samples:", n_tumor, "| Normal samples:", n_normal, "\n")
cat("  True DE genes (injected):", length(de_genes_idx), "\n")
cat("  Known cancer genes in DE set: ~20\n\n")

# --- Step 4: Run Consensus Differential Expression ---------------------------
cat("--- Step 4: Consensus Differential Expression ---\n")
#
# We run three independent DE methods and compute consensus:
#   - DESeq2 (Love et al., Genome Biology 2014)
#   - edgeR (Robinson et al., Bioinformatics 2010)
#   - limma-voom (Law et al., Genome Biology 2014)
#
# Consensus reduces method-specific biases (Schurch et al., RNA 2016).

de_results <- consensus_de(
  base_counts, col_data,
  methods = c("deseq2", "edger", "limma"),
  padj_threshold = 0.05,
  min_methods = 2
)

n_consensus <- sum(de_results$n_methods >= 2)
n_deseq2 <- if ("sig_deseq2" %in% names(de_results)) sum(de_results$sig_deseq2) else 0
n_edger <- if ("sig_edger" %in% names(de_results)) sum(de_results$sig_edger) else 0
n_limma <- if ("sig_limma" %in% names(de_results)) sum(de_results$sig_limma) else 0

cat("  DESeq2 significant:", n_deseq2, "\n")
cat("  edgeR significant:", n_edger, "\n")
cat("  limma-voom significant:", n_limma, "\n")
cat("  Consensus (>=2 methods):", n_consensus, "\n\n")

# --- Step 5: Survival Analysis -----------------------------------------------
cat("--- Step 5: Survival Analysis ---\n")
#
# Cox proportional hazards regression for top DE genes.
# Reference: Cox, D.R. (1972). Regression models and life-tables.
# Journal of the Royal Statistical Society: Series B, 34(2), 187-202.

top_survival_genes <- head(de_results$gene, 100)
common_samples <- intersect(colnames(base_counts), clinical$sample_id)

surv_results <- calculate_hazard_ratios(
  base_counts[top_survival_genes, common_samples],
  time = clinical$time[match(common_samples, clinical$sample_id)],
  status = clinical$status[match(common_samples, clinical$sample_id)],
  genes = top_survival_genes
)

n_sig_surv <- sum(surv_results$padj < 0.05, na.rm = TRUE)
cat("  Genes tested:", nrow(surv_results), "\n")
cat("  Significant (padj < 0.05):", n_sig_surv, "\n")

# Top prognostic genes
if (nrow(surv_results) > 0) {
  top_surv <- head(surv_results[order(surv_results$pvalue), ], 5)
  cat("  Top prognostic genes:\n")
  for (i in seq_len(nrow(top_surv))) {
    cat(sprintf("    %s: HR=%.2f (95%% CI: %.2f-%.2f), p=%.2e\n",
                top_surv$gene[i], top_surv$hr[i],
                top_surv$hr_lower[i], top_surv$hr_upper[i],
                top_surv$pvalue[i]))
  }
}
cat("\n")

# --- Step 6: Target Scoring & Ranking ---------------------------------------
cat("--- Step 6: Composite Target Scoring ---\n")
#
# Multi-dimensional scoring integrating:
#   - DE significance (-log10 padj)
#   - Survival association (-log10 Cox p-value)
#   - Network centrality (composite centrality metric)
#
# Clamped to [0, 1] per dimension.

scored_targets <- composite_target_score(
  de_results = de_results,
  survival_hr = surv_results,
  centrality = centrality
)

ranked_targets <- rank_targets(scored_targets, top_n = 30)

cat("  Total targets scored:", nrow(scored_targets), "\n")
cat("  Top 30 ranked\n\n")

cat("  Top 10 Targets:\n")
cat("  ", paste(rep("-", 60), collapse = ""), "\n")
top10 <- head(ranked_targets, 10)
for (i in seq_len(nrow(top10))) {
  cat(sprintf("  %2d. %-12s  Score: %.3f  Tier: %s\n",
              i, top10$gene[i], top10$composite_score[i], top10$target_tier[i]))
}
cat("\n")

# --- Step 7: Evaluate Against Gold Standards --------------------------------
cat("--- Step 7: Benchmarking Against COSMIC/CGC ---\n")
#
# We evaluate predictions against COSMIC Cancer Gene Census (Sondka 2018)
# using precision, recall, AUROC, and precision@K metrics.

eval_agent <- EvaluationAgent$new()
metrics <- eval_agent$evaluate(ranked_targets, gs)

cat("\n", eval_agent$report(), "\n\n")

# Per-source breakdown
if (!is.null(metrics$per_source)) {
  cat("  Per-source breakdown:\n")
  for (src in names(metrics$per_source)) {
    m <- metrics$per_source[[src]]
    cat(sprintf("    %s: precision=%.3f, recall=%.3f (%d/%d recovered)\n",
                src, m$precision, m$recall, m$n_recovered, m$n_known))
  }
}
cat("\n")

# --- Step 8: Gap Detection ---------------------------------------------------
cat("--- Step 8: Gap Detection ---\n")
#
# The GapDetectionAgent identifies missing data, low statistical power,
# and conflicting signals that could compromise target identification.

gap_agent <- GapDetectionAgent$new(
  min_samples = 30,
  min_degs = 50,
  min_modules = 2
)

analysis_results <- list(
  de = de_results,
  survival = surv_results,
  network = list(network = g, centrality = centrality, modules = detect_network_modules(g)),
  scoring = scored_targets
)

analysis_data <- list(
  expression = base_counts,
  col_data = col_data,
  clinical = clinical,
  ppi = ppi_data
)

gaps <- gap_agent$detect(analysis_results, analysis_data)
gap_summary <- gap_agent$summary()

cat("  Gaps detected:", gap_summary$n_gaps, "\n")
cat("  Overall severity:", gap_summary$severity, "\n")

if (gap_summary$n_gaps > 0) {
  cat("  Gap details:\n")
  for (gap in gaps) {
    cat(sprintf("    [%s] %s\n", toupper(gap$severity), gap$description))
  }
}

cat("\n  Recommendations:\n")
for (rec in gap_summary$recommendations) {
  cat("    -", rec, "\n")
}
cat("\n")

# --- Step 9: Hypothesis Generation ------------------------------------------
cat("--- Step 9: Hypothesis Generation ---\n")

hyp_agent <- HypothesisAgent$new(top_n_hypotheses = 10)
hypotheses <- hyp_agent$generate(scored_targets)

cat("  Generated", nrow(hypotheses), "hypotheses\n\n")
top_hyp <- head(hypotheses, 5)
for (i in seq_len(nrow(top_hyp))) {
  cat(sprintf("  %s (confidence: %.0f%%): %s\n",
              top_hyp$gene[i], top_hyp$confidence[i] * 100, top_hyp$mechanism[i]))
}
cat("\n")

# --- Step 10: Export Results ------------------------------------------------
cat("--- Step 10: Export Results ---\n")

output_dir <- file.path("demo", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- create_pipeline_config(
  project = "TCGA-BRCA-Demo",
  cancer_type = "breast cancer",
  output_dir = output_dir
)

pipeline_result <- list(
  targets = ranked_targets,
  hypotheses = hypotheses,
  gaps = gap_summary,
  evaluation = metrics,
  results = analysis_results,
  config = config,
  timestamp = Sys.time()
)

exported_files <- export_results(pipeline_result, output_dir)
cat("  Exported", length(exported_files), "files to", output_dir, "\n")
for (f in exported_files) {
  cat("    -", basename(f), "\n")
}
cat("\n")

# --- Summary -----------------------------------------------------------------
cat("================================================================\n")
cat("  Demo Complete\n")
cat("================================================================\n")
cat(sprintf("  DE genes (consensus):     %d\n", n_consensus))
cat(sprintf("  Survival-significant:     %d\n", n_sig_surv))
cat(sprintf("  Targets ranked:           %d\n", nrow(ranked_targets)))
cat(sprintf("  COSMIC genes recovered:   %d\n",
            sum(ranked_targets$gene %in% gs$genes)))
cat(sprintf("  Gaps detected:            %d (%s)\n",
            gap_summary$n_gaps, gap_summary$severity))
cat(sprintf("  Hypotheses generated:     %d\n", nrow(hypotheses)))
cat("================================================================\n")
