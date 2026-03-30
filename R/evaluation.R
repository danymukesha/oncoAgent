#' @title Evaluation Module
#' @description Benchmarking functions for comparing predicted targets against
#'   gold standard gene sets. Computes AUROC, precision@K, recall, and
#'   generates comprehensive evaluation reports.

#' @name evaluate_targets
#' @title Evaluate Predicted Targets
#' @description Computes precision of predicted targets against a known set.
#' @param predicted Character vector of predicted gene symbols
#' @param known_targets Character vector of known gold standard genes
#' @return Precision value (0-1)
#' @export
evaluate_targets <- function(predicted, known_targets) {
  tp <- sum(predicted %in% known_targets)
  precision <- tp / length(predicted)

  message("Evaluation: ", tp, " true positives out of ",
          length(predicted), " predicted (precision = ",
          round(precision, 3), ")")

  precision
}

#' @name benchmark_vs_goldstandard
#' @title Benchmark Against Gold Standard
#' @description Comprehensive benchmarking against multiple gold standard
#'   gene sets with stratification.
#' @param scored_targets Output from composite_target_score or rank_targets
#' @param gold_standard Named list with $genes and optionally $per_source
#' @param n_permutations Number of permutation tests for significance
#' @return List with precision, recall, F1, permutation p-values
#' @export
benchmark_vs_goldstandard <- function(scored_targets,
                                      gold_standard,
                                      n_permutations = 1000) {

  predicted <- scored_targets$gene
  known <- gold_standard$genes

  tp <- sum(predicted %in% known)
  fp <- sum(!predicted %in% known)
  fn <- sum(!known %in% predicted)

  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0

  set.seed(42)
  null_precision <- replicate(n_permutations, {
    random_pred <- sample(
      unique(c(predicted, known)),
      length(predicted),
      replace = FALSE
    )
    sum(random_pred %in% known) / length(random_pred)
  })

  perm_pvalue <- mean(null_precision >= precision)

  source_metrics <- NULL
  if (!is.null(gold_standard$per_source)) {
    source_metrics <- lapply(names(gold_standard$per_source), function(src) {
      src_genes <- gold_standard$per_source[[src]]
      src_tp <- sum(predicted %in% src_genes)
      list(
        source = src,
        precision = src_tp / length(predicted),
        recall = src_tp / length(src_genes),
        n_known = length(src_genes),
        n_recovered = src_tp
      )
    })
  }

  result <- list(
    precision = precision,
    recall = recall,
    f1 = f1,
    true_positives = tp,
    false_positives = fp,
    false_negatives = fn,
    permutation_pvalue = perm_pvalue,
    null_precision_mean = mean(null_precision),
    null_precision_sd = stats::sd(null_precision),
    source_metrics = source_metrics,
    n_predicted = length(predicted),
    n_known = length(known)
  )

  message("Benchmark: Precision=", round(precision, 3),
          " Recall=", round(recall, 3),
          " F1=", round(f1, 3),
          " (permutation p=", round(perm_pvalue, 4), ")")

  result
}

#' @name calculate_auroc
#' @title Calculate AUROC
#' @description Computes Area Under the Receiver Operating Characteristic
#'   curve for target prediction.
#' @param scored_targets data.frame with columns: gene, composite_score
#' @param known_targets Character vector of known gold standard genes
#' @return AUROC value (0-1)
#' @export
calculate_auroc <- function(scored_targets, known_targets) {

  df <- data.frame(
    gene = scored_targets$gene,
    score = scored_targets$composite_score,
    is_known = scored_targets$gene %in% known_targets,
    stringsAsFactors = FALSE
  )

  if (length(unique(df$is_known)) < 2) {
    message("AUROC: Cannot compute (only one class present)")
    return(NA_real_)
  }

  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(
      response = df$is_known,
      predictor = df$score,
      quiet = TRUE,
      direction = "<"
    )
    auroc <- as.numeric(pROC::auc(roc_obj))
  } else {
    auroc <- .manual_auroc(df$score, df$is_known)
  }

  message("AUROC: ", round(auroc, 3))
  auroc
}

#' @keywords internal
.manual_auroc <- function(scores, labels) {
  n_pos <- sum(labels)
  n_neg <- sum(!labels)

  if (n_pos == 0 || n_neg == 0) return(0.5)

  scores_pos <- scores[labels]
  scores_neg <- scores[!labels]

  auc <- sum(outer(scores_pos, scores_neg, ">")) + 0.5 * sum(outer(scores_pos, scores_neg, "=="))
  auc / (n_pos * n_neg)
}

#' @name precision_at_k
#' @title Precision at K
#' @description Computes precision for the top K predicted targets.
#' @param scored_targets data.frame with columns: gene, composite_score
#' @param known_targets Character vector of known gold standard genes
#' @param k Number of top predictions to consider
#' @return Precision value (0-1)
#' @export
precision_at_k <- function(scored_targets, known_targets, k = 10) {

  top_k <- utils::head(scored_targets$gene, k)
  tp <- sum(top_k %in% known_targets)
  precision <- tp / k

  message("Precision@", k, ": ", round(precision, 3),
          " (", tp, "/", k, ")")
  precision
}

#' @name recall_of_known_targets
#' @title Recall of Known Targets
#' @description Computes recall (sensitivity) for recovering known targets.
#' @param predicted Character vector of predicted genes
#' @param known_targets Character vector of known gold standard genes
#' @return Recall value (0-1)
#' @export
recall_of_known_targets <- function(predicted, known_targets) {
  recovered <- sum(predicted %in% known_targets)
  recall <- recovered / length(known_targets)

  message("Recall: ", round(recall, 3),
          " (", recovered, "/", length(known_targets), " known targets recovered)")
  recall
}

#' @name generate_evaluation_report
#' @title Generate Evaluation Report
#' @description Creates a comprehensive markdown evaluation report
#'   comparing predicted targets against gold standards.
#' @param scored_targets Output from composite_target_score
#' @param gold_standard Named list with $genes and $per_source
#' @param output_file Optional file path to write report
#' @return Character string with the report
#' @export
generate_evaluation_report <- function(scored_targets,
                                       gold_standard,
                                       output_file = NULL) {

  predicted <- scored_targets$gene

  bench <- benchmark_vs_goldstandard(scored_targets, gold_standard)
  auroc <- calculate_auroc(scored_targets, gold_standard$genes)

  p10 <- precision_at_k(scored_targets, gold_standard$genes, k = 10)
  p20 <- precision_at_k(scored_targets, gold_standard$genes, k = 20)
  p50 <- precision_at_k(scored_targets, gold_standard$genes, k = 50)

  top_recovered <- intersect(head(predicted, 20), gold_standard$genes)
  top_missed <- setdiff(head(gold_standard$genes, 20), predicted)

  report <- c(
    "# oncoAgent Evaluation Report",
    paste0("Generated: ", Sys.time()),
    "",
    "## Summary",
    "",
    paste0("- **Total targets predicted**: ", bench$n_predicted),
    paste0("- **Known targets in gold standard**: ", bench$n_known),
    paste0("- **True positives**: ", bench$true_positives),
    paste0("- **False positives**: ", bench$false_positives),
    paste0("- **False negatives**: ", bench$false_negatives),
    "",
    "## Performance Metrics",
    "",
    paste0("- **Precision**: ", round(bench$precision, 3)),
    paste0("- **Recall**: ", round(bench$recall, 3)),
    paste0("- **F1 Score**: ", round(bench$f1, 3)),
    paste0("- **AUROC**: ", round(auroc, 3)),
    paste0("- **Permutation p-value**: ", round(bench$permutation_pvalue, 4)),
    "",
    "## Precision@K",
    "",
    paste0("- **P@10**: ", round(p10, 3)),
    paste0("- **P@20**: ", round(p20, 3)),
    paste0("- **P@50**: ", round(p50, 3)),
    "",
    "## Top Recovered Targets",
    "",
    paste0("  ", paste(top_recovered, collapse = ", ")),
    "",
    "## Top Missed Targets",
    "",
    paste0("  ", paste(top_missed, collapse = ", "))
  )

  if (!is.null(bench$source_metrics)) {
    report <- c(report, "", "## Per-Source Metrics", "")
    for (sm in bench$source_metrics) {
      report <- c(report,
        paste0("### ", sm$source),
        paste0("- Precision: ", round(sm$precision, 3)),
        paste0("- Recall: ", round(sm$recall, 3)),
        paste0("- Recovered: ", sm$n_recovered, "/", sm$n_known),
        ""
      )
    }
  }

  report_text <- paste(report, collapse = "\n")

  if (!is.null(output_file)) {
    writeLines(report_text, output_file)
    message("Report written to ", output_file)
  }

  invisible(report_text)
}
