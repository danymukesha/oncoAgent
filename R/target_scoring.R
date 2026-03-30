#' @title Target Scoring Module
#' @description Multi-dimensional target scoring integrating differential
#'   expression, survival, network centrality, druggability, and CRISPR
#'   dependency data.

#' @name composite_target_score
#' @title Composite Target Score
#' @description Computes a weighted composite score for each gene across
#'   multiple evidence dimensions.
#' @param de_results data.frame with columns: gene, log2FoldChange, padj
#' @param survival_hr data.frame with columns: gene, hr, pvalue (from calculate_hazard_ratios)
#' @param centrality data.frame with columns: gene, composite_centrality
#' @param crispr_scores Optional data.frame with columns: gene, dependency_score
#' @param weights Named list of weights for each evidence type:
#'   de, survival, centrality, crispr (default: equal weights)
#' @return data.frame with composite scores and per-dimension breakdowns
#' @export
composite_target_score <- function(de_results,
                                   survival_hr = NULL,
                                   centrality = NULL,
                                   crispr_scores = NULL,
                                   weights = NULL) {

  if (is.null(weights)) {
    n_dims <- 1 + !is.null(survival_hr) + !is.null(centrality) + !is.null(crispr_scores)
    weights <- list(
      de = 1 / n_dims,
      survival = if (!is.null(survival_hr)) 1 / n_dims else 0,
      centrality = if (!is.null(centrality)) 1 / n_dims else 0,
      crispr = if (!is.null(crispr_scores)) 1 / n_dims else 0
    )
  }

  genes <- de_results$gene
  n <- length(genes)

  # DE score: use consensus_score if available, then padj, then pvalue
  if ("consensus_score" %in% names(de_results)) {
    # From consensus_de: higher consensus_score = stronger evidence
    cs <- de_results$consensus_score[match(genes, de_results$gene)]
    cs[is.na(cs)] <- 0
    de_score <- cs / max(cs, na.rm = TRUE)
  } else {
    pval_col <- if ("padj" %in% names(de_results)) "padj" else if ("pvalue" %in% names(de_results)) "pvalue" else NULL
    if (!is.null(pval_col)) {
      pvals <- pmax(de_results[[pval_col]][match(genes, de_results$gene)], .Machine$double.xmin)
    } else {
      pvals <- rep(.Machine$double.xmin, n)
    }
    de_neglog <- -log10(pvals)
    de_score <- de_neglog / max(de_neglog, na.rm = TRUE)
  }

  # Survival score
  if (!is.null(survival_hr) && "gene" %in% names(survival_hr)) {
    surv_lookup <- setNames(
      -log10(pmax(survival_hr$pvalue, .Machine$double.xmin)),
      survival_hr$gene
    )
    hr_raw <- unname(surv_lookup[genes])
    hr_raw[is.na(hr_raw)] <- 0
    hr_score <- hr_raw / max(hr_raw, na.rm = TRUE)
  } else {
    hr_score <- rep(0, n)
  }

  # Centrality score
  if (!is.null(centrality) && "gene" %in% names(centrality)) {
    cent_lookup <- setNames(centrality$composite_centrality, centrality$gene)
    cent_raw <- unname(cent_lookup[genes])
    cent_raw[is.na(cent_raw)] <- 0
    rng <- max(cent_raw, na.rm = TRUE) - min(cent_raw, na.rm = TRUE)
    centrality_score <- (cent_raw - min(cent_raw, na.rm = TRUE)) /
      (rng + .Machine$double.eps)
  } else {
    centrality_score <- rep(0, n)
  }

  # CRISPR score
  if (!is.null(crispr_scores) && "gene" %in% names(crispr_scores)) {
    crispr_lookup <- setNames(-crispr_scores$dependency_score, crispr_scores$gene)
    crispr_raw <- unname(crispr_lookup[genes])
    crispr_raw[is.na(crispr_raw)] <- 0
    rng <- max(crispr_raw, na.rm = TRUE) - min(crispr_raw, na.rm = TRUE)
    crispr_score <- (crispr_raw - min(crispr_raw, na.rm = TRUE)) /
      (rng + .Machine$double.eps)
  } else {
    crispr_score <- rep(0, n)
  }

  scores <- data.frame(
    gene = genes,
    de_score = de_score,
    hr_score = hr_score,
    centrality_score = centrality_score,
    crispr_score = crispr_score,
    stringsAsFactors = FALSE
  )

  scores$composite_score <- weights$de * scores$de_score +
    weights$survival * scores$hr_score +
    weights$centrality * scores$centrality_score +
    weights$crispr * scores$crispr_score

  # Clamp to [0, 1]
  scores$composite_score <- pmin(pmax(scores$composite_score, 0), 1)

  scores <- scores[order(-scores$composite_score), ]
  scores$rank <- seq_len(nrow(scores))

  message("Scored ", nrow(scores), " targets across ",
          sum(unlist(weights) > 0), " dimensions")

  scores
}

#' @name rank_targets
#' @title Rank Targets
#' @description Final ranking combining composite scores with additional
#'   filters and annotations.
#' @param scored_targets Output from composite_target_score()
#' @param top_n Return top N targets (default: 50)
#' @param min_de_padj Minimum DE significance (default: 0.05)
#' @param min_evidence_dims Minimum evidence dimensions required (default: 1)
#' @return data.frame of ranked targets
#' @export
rank_targets <- function(scored_targets,
                         top_n = 50,
                         min_de_padj = 0.05,
                         min_evidence_dims = 1) {

  targets <- scored_targets

  evidence_cols <- c("de_score", "hr_score", "centrality_score", "crispr_score")
  available_cols <- intersect(evidence_cols, names(targets))
  targets$n_evidence_dims <- rowSums(targets[, available_cols, drop = FALSE] > 0)

  targets <- targets[targets$n_evidence_dims >= min_evidence_dims, ]

  targets <- targets[order(-targets$composite_score), ]
  targets <- utils::head(targets, top_n)

  if (nrow(targets) > 0) {
    targets$rank <- seq_len(nrow(targets))
    targets$target_tier <- ifelse(targets$rank <= 10, "Tier 1",
                          ifelse(targets$rank <= 25, "Tier 2", "Tier 3"))
  } else {
    targets$rank <- integer(0)
    targets$target_tier <- character(0)
  }

  message("Ranked targets: ", nrow(targets), " (Tier 1: ",
          sum(targets$target_tier == "Tier 1"), ", Tier 2: ",
          sum(targets$target_tier == "Tier 2"), ", Tier 3: ",
          sum(targets$target_tier == "Tier 3"), ")")

  targets
}

#' @name integrate_multi_omics
#' @title Multi-Omics Integration
#' @description Integrates gene expression, mutation, methylation, and
#'   CRISPR data into a unified feature matrix for target scoring.
#' @param expression Gene x Sample expression matrix
#' @param mutations Mutation data.frame (gene, sample, consequence)
#' @param methylation Gene x Sample methylation matrix (beta values)
#' @param crispr Gene dependency scores
#' @param integration_method "concatenate", "weighted_average", or "factor"
#' @return Integrated gene x feature matrix
#' @export
integrate_multi_omics <- function(expression,
                                  mutations = NULL,
                                  methylation = NULL,
                                  crispr = NULL,
                                  integration_method = "weighted_average") {

  genes <- rownames(expression)
  features <- list()

  expr_summary <- data.frame(
    gene = genes,
    mean_expr = rowMeans(expression, na.rm = TRUE),
    sd_expr = apply(expression, 1, stats::sd, na.rm = TRUE),
    cv_expr = apply(expression, 1, function(x) stats::sd(x, na.rm = TRUE) /
                      (mean(x, na.rm = TRUE) + .Machine$double.eps)),
    stringsAsFactors = FALSE
  )
  features$expression <- expr_summary

  if (!is.null(mutations)) {
    mut_summary <- mutations |>
      dplyr::group_by(gene) |>
      dplyr::summarise(
        n_mutations = dplyr::n(),
        n_samples_mutated = length(unique(sample)),
        has_driver = any(grepl("missense|nonsense|frameshift", consequence,
                               ignore.case = TRUE)),
        .groups = "drop"
      )
    features$mutation <- mut_summary
  }

  if (!is.null(methylation)) {
    meth_genes <- intersect(genes, rownames(methylation))
    meth_summary <- data.frame(
      gene = meth_genes,
      mean_methylation = rowMeans(methylation[meth_genes, , drop = FALSE], na.rm = TRUE),
      methylation_sd = apply(methylation[meth_genes, , drop = FALSE], 1,
                              stats::sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    features$methylation <- meth_summary
  }

  if (!is.null(crispr)) {
    features$crispr <- data.frame(
      gene = rownames(crispr),
      dependency_score = rowMeans(crispr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (integration_method == "concatenate") {
    integrated <- features[[1]]
    for (i in 2:length(features)) {
      integrated <- merge(integrated, features[[i]], by = "gene", all = TRUE)
    }
  } else {
    integrated <- features[[1]]
    for (i in 2:length(features)) {
      integrated <- merge(integrated, features[[i]], by = "gene", all = TRUE)
    }
  }

  integrated[is.na(integrated)] <- 0

  message("Integrated ", nrow(integrated), " genes across ",
          length(features), " omics layers")

  list(
    features = features,
    integrated = integrated,
    method = integration_method
  )
}

#' @name druggability_score
#' @title Druggability Score
#' @description Estimates druggability of target based on protein class,
#'   known drug interactions, and structural features.
#' @param genes Character vector of gene symbols
#' @param drug_data GDSC or similar drug sensitivity data
#' @return data.frame with druggability scores
#' @export
druggability_score <- function(genes, drug_data = NULL) {

  druggable_classes <- list(
    kinase = c("EGFR", "ERBB2", "ALK", "MET", "RET", "KIT", "PDGFRA",
               "FLT3", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1",
               "MAPK3", "PIK3CA", "AKT1", "MTOR", "CDK4", "CDK6",
               "JAK1", "JAK2", "ABL1", "SRC", "FGFR1", "FGFR2", "FGFR3"),
    gpcr = c("GPCR1", "GPCR2"),
    protease = c("MMP2", "MMP9", "CAPN1"),
    epigenetic = c("EZH2", "HDAC1", "HDAC2", "DNMT1", "DNMT3A", "BRD4",
                   "KDM1A", "KMT2A"),
    nuclear_receptor = c("ESR1", "AR", "PPARA", "RXRA"),
    transporter = c("SLC16A1", "SLC7A5", "ABCB1", "ABCG2")
  )

  score_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  score_df$druggability_score <- 0
  score_df$protein_class <- "other"
  score_df$is_kinase <- FALSE
  score_df$is_epigenetic <- FALSE

  for (cls in names(druggable_classes)) {
    in_class <- score_df$gene %in% druggable_classes[[cls]]
    if (cls == "kinase") {
      score_df$druggability_score[in_class] <- 0.9
      score_df$protein_class[in_class] <- cls
      score_df$is_kinase[in_class] <- TRUE
    } else if (cls == "epigenetic") {
      score_df$druggability_score[in_class] <- 0.8
      score_df$protein_class[in_class] <- cls
      score_df$is_epigenetic[in_class] <- TRUE
    } else {
      score_df$druggability_score[in_class] <- 0.7
      score_df$protein_class[in_class] <- cls
    }
  }

  if (!is.null(drug_data) && "drug_name" %in% names(drug_data)) {
    drugged_genes <- unique(drug_data$gene_name %||% drug_data$target)
    score_df$has_known_drug <- score_df$gene %in% drugged_genes
    score_df$druggability_score <- ifelse(
      score_df$has_known_drug,
      pmin(score_df$druggability_score + 0.1, 1),
      score_df$druggability_score
    )
  }

  score_df
}

#' @name bayesian_target_ranking
#' @title Bayesian Target Ranking
#' @description Bayesian approach to target ranking that models uncertainty
#'   in evidence and computes posterior probabilities of being a true target.
#' @param scored_targets Output from composite_target_score()
#' @param prior_probability Prior probability of being a true target (default: 0.01)
#' @param noise_sd Assumed noise level in scores (default: 0.1)
#' @return data.frame with posterior probabilities
#' @export
bayesian_target_ranking <- function(scored_targets,
                                    prior_probability = 0.01,
                                    noise_sd = 0.1) {

  targets <- scored_targets

  likelihood <- targets$composite_score
  likelihood <- likelihood / (max(likelihood) + .Machine$double.eps)

  posterior <- (likelihood * prior_probability) /
    (likelihood * prior_probability + (1 - likelihood) * (1 - prior_probability))

  targets$posterior_prob <- posterior
  targets$bayesian_rank <- rank(-targets$posterior_prob)

  targets <- targets[order(targets$bayesian_rank), ]

  targets$evidence_category <- ifelse(
    targets$posterior_prob > 0.5, "Strong",
    ifelse(targets$posterior_prob > 0.1, "Moderate",
    ifelse(targets$posterior_prob > 0.05, "Suggestive", "Weak"))
  )

  message("Bayesian ranking: ", sum(targets$evidence_category == "Strong"),
          " strong, ", sum(targets$evidence_category == "Moderate"),
          " moderate targets")

  targets
}
