#' @title Differential Expression Module
#' @description Unified interface for differential expression analysis using
#'   DESeq2, edgeR, and limma-voom. Includes consensus ranking across methods.

#' @name run_deseq2
#' @title Run DESeq2 Differential Expression
#' @description Performs differential expression analysis using DESeq2 with
#'   shrinkage estimation and multiple testing correction.
#' @param count_matrix Gene x Sample count matrix (integer)
#' @param col_data Sample metadata with 'condition' column
#' @param contrast Character vector of length 3: c("condition", "tumor", "normal")
#' @param lfc_threshold Log2 fold-change threshold for significance (default: 0)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param shrink Logical, apply apeglm shrinkage (default: TRUE)
#' @return data.frame with DE results sorted by adjusted p-value
#' @export
run_deseq2 <- function(count_matrix,
                       col_data,
                       contrast = NULL,
                       lfc_threshold = 0,
                       padj_threshold = 0.05,
                       shrink = TRUE) {
    storage.mode(count_matrix) <- "integer"

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = col_data,
        design = ~condition
    )

    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    if (!is.null(contrast)) {
        res <- DESeq2::results(
            dds,
            contrast = contrast,
            lfcThreshold = lfc_threshold,
            alpha = padj_threshold
        )
    } else {
        res <- DESeq2::results(
            dds,
            lfcThreshold = lfc_threshold,
            alpha = padj_threshold
        )
    }

    if (shrink && requireNamespace("apeglm", quietly = TRUE)) {
        coef_name <- DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))]
        res <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", res = res) # changed from apeglm::lfcShrink
    }

    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL
    res_df <- res_df[order(res_df$padj), ]

    message(
        "DESeq2: ", sum(res_df$padj < padj_threshold & !is.na(res_df$padj)),
        " significant DEGs (padj < ", padj_threshold, ")"
    )

    res_df
}

#' @name run_edger
#' @title Run edgeR Differential Expression
#' @description Performs differential expression using edgeR's quasi-likelihood
#'   pipeline with TMM normalization.
#' @param count_matrix Gene x Sample count matrix (integer)
#' @param col_data Sample metadata with 'condition' column
#' @param contrast_levels Character vector: c("tumor", "normal")
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @return data.frame with DE results
#' @export
run_edger <- function(count_matrix,
                      col_data,
                      contrast_levels = NULL,
                      padj_threshold = 0.05) {
    y <- edgeR::DGEList(counts = count_matrix, group = col_data$condition)
    y <- edgeR::calcNormFactors(y, method = "TMM")

    design <- stats::model.matrix(~ 0 + condition, data = col_data)
    colnames(design) <- levels(col_data$condition)

    if (is.null(contrast_levels)) {
        contrast_levels <- rev(levels(col_data$condition))
    }

    contrast_matrix <- limma::makeContrasts(
        contrasts = paste(contrast_levels, collapse = " - "),
        levels = design
    )

    y <- edgeR::estimateDisp(y, design)

    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix[, 1])

    res <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    res$gene <- rownames(res)
    rownames(res) <- NULL

    names(res)[names(res) == "PValue"] <- "pvalue"
    names(res)[names(res) == "FDR"] <- "padj"
    names(res)[names(res) == "logFC"] <- "log2FoldChange"

    message(
        "edgeR: ", sum(res$padj < padj_threshold & !is.na(res$padj)),
        " significant DEGs (FDR < ", padj_threshold, ")"
    )

    res
}

#' @name run_limma_voom
#' @title Run limma-voom Differential Expression
#' @description Performs differential expression using limma-voom with
#'   sample quality weights.
#' @param count_matrix Gene x Sample count matrix (integer)
#' @param col_data Sample metadata with 'condition' column
#' @param contrast_levels Character vector: c("tumor", "normal")
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @return data.frame with DE results
#' @export
run_limma_voom <- function(count_matrix,
                           col_data,
                           contrast_levels = NULL,
                           padj_threshold = 0.05) {
    y <- edgeR::DGEList(counts = count_matrix, group = col_data$condition)
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~ 0 + condition, data = col_data)
    colnames(design) <- levels(col_data$condition)

    if (is.null(contrast_levels)) {
        contrast_levels <- rev(levels(col_data$condition))
    }

    contrast_matrix <- limma::makeContrasts(
        contrasts = paste(contrast_levels, collapse = " - "),
        levels = design
    )

    v <- limma::voom(y, design, plot = FALSE)
    vfit <- limma::lmFit(v, design)
    vfit <- limma::contrasts.fit(vfit, contrasts = contrast_matrix)
    efit <- limma::eBayes(vfit)

    res <- limma::topTable(efit, sort.by = "P", n = Inf)
    res$gene <- rownames(res)
    rownames(res) <- NULL

    names(res)[names(res) == "P.Value"] <- "pvalue"
    names(res)[names(res) == "adj.P.Val"] <- "padj"
    names(res)[names(res) == "AveExpr"] <- "baseMean"
    names(res)[names(res) == "t"] <- "stat"

    message(
        "limma-voom: ", sum(res$padj < padj_threshold & !is.na(res$padj)),
        " significant DEGs (adj.P.Val < ", padj_threshold, ")"
    )

    res
}

#' @name consensus_de
#' @title Consensus Differential Expression
#' @description Runs multiple DE methods and computes a consensus ranking.
#'   Genes detected by multiple methods are ranked higher.
#' @param count_matrix Gene x Sample count matrix
#' @param col_data Sample metadata with 'condition' column
#' @param methods Character vector of methods: "deseq2", "edger", "limma"
#' @param padj_threshold Adjusted p-value threshold
#' @param min_methods Minimum number of methods that must agree (default: 2)
#' @return data.frame with consensus results
#' @export
consensus_de <- function(count_matrix,
                         col_data,
                         methods = c("deseq2", "edger", "limma"),
                         padj_threshold = 0.05,
                         min_methods = 2) {
    all_results <- list()

    if ("deseq2" %in% methods) {
        all_results$deseq2 <- tryCatch(
            run_deseq2(count_matrix, col_data, padj_threshold = padj_threshold),
            error = function(e) {
                message("DESeq2 failed: ", e$message)
                NULL
            }
        )
    }

    if ("edger" %in% methods) {
        all_results$edger <- tryCatch(
            run_edger(count_matrix, col_data, padj_threshold = padj_threshold),
            error = function(e) {
                message("edgeR failed: ", e$message)
                NULL
            }
        )
    }

    if ("limma" %in% methods) {
        all_results$limma <- tryCatch(
            run_limma_voom(count_matrix, col_data, padj_threshold = padj_threshold),
            error = function(e) {
                message("limma failed: ", e$message)
                NULL
            }
        )
    }

    all_results <- purrr::compact(all_results)

    if (length(all_results) == 0) {
        stop("All DE methods failed")
    }

    sig_lists <- lapply(names(all_results), function(m) {
        res <- all_results[[m]]
        res$gene[res$padj < padj_threshold & !is.na(res$padj)]
    })
    names(sig_lists) <- names(all_results)

    all_genes <- unique(unlist(sig_lists))

    consensus <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

    for (m in names(sig_lists)) {
        consensus[[paste0("sig_", m)]] <- consensus$gene %in% sig_lists[[m]]
    }

    sig_cols <- paste0("sig_", names(sig_lists))
    consensus$n_methods <- rowSums(consensus[, sig_cols, drop = FALSE])

    rank_cols <- lapply(names(all_results), function(m) {
        res <- all_results[[m]]
        res_df <- data.frame(
            gene = res$gene, rank = seq_len(nrow(res)),
            stringsAsFactors = FALSE
        )
        names(res_df)[2] <- paste0("rank_", m)
        res_df
    })

    consensus <- Reduce(
        function(x, y) merge(x, y, by = "gene", all.x = TRUE),
        c(list(consensus), rank_cols)
    )

    rank_only_cols <- grep("^rank_", names(consensus), value = TRUE)
    consensus$mean_rank <- rowMeans(consensus[, rank_only_cols, drop = FALSE],
        na.rm = TRUE
    )

    consensus <- consensus[order(consensus$mean_rank), ]

    consensus$consensus_score <- consensus$n_methods / length(all_results) *
        (1 / (1 + log1p(consensus$mean_rank)))

    consensus <- consensus[order(-consensus$consensus_score), ]

    top_de <- consensus[consensus$n_methods >= min_methods, ]

    message(
        "Consensus DE: ", nrow(top_de), " genes detected by >= ",
        min_methods, " methods"
    )

    consensus
}

#' @name rank_de_genes
#' @title Rank DE Genes by Combined Score
#' @description Ranks DE results by a combination of statistical significance
#'   and effect size.
#' @param de_results data.frame with columns: gene, log2FoldChange, padj
#' @param w_pvalue Weight for p-value component (default: 0.5)
#' @param w_lfc Weight for log-fold-change component (default: 0.5)
#' @return data.frame with added rank and combined_score columns
#' @export
rank_de_genes <- function(de_results, w_pvalue = 0.5, w_lfc = 0.5) {
    df <- de_results

    df$neg_log_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
    df$abs_lfc <- abs(df$log2FoldChange)

    df$score_pvalue <- (df$neg_log_padj - min(df$neg_log_padj)) /
        (max(df$neg_log_padj) - min(df$neg_log_padj) + .Machine$double.eps)
    df$score_lfc <- (df$abs_lfc - min(df$abs_lfc)) /
        (max(df$abs_lfc) - min(df$abs_lfc) + .Machine$double.eps)

    df$combined_score <- w_pvalue * df$score_pvalue + w_lfc * df$score_lfc
    df$rank <- rank(-df$combined_score)

    df[order(df$rank), ]
}
