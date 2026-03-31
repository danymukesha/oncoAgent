#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @keywords internal
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts vst
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmLRT cpm glmQLFit glmQLFTest topTags
#' @importFrom limma voom lmFit eBayes topTable makeContrasts contrasts.fit removeBatchEffect
#' @importFrom survival Surv coxph survfit
#' @importFrom igraph graph_from_data_frame cluster_louvain degree betweenness page_rank
#' @importFrom dplyr filter select mutate arrange left_join group_by summarise pull rename bind_rows across everything desc
#' @importFrom tidyr pivot_wider pivot_longer drop_na
#' @importFrom purrr map map_dfr map_chr map_lgl keep compact reduce safely possibly
#' @importFrom stringr str_detect str_extract str_replace str_to_lower
#' @importFrom data.table as.data.table rbindlist
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline theme_minimal labs scale_color_manual
#' @importFrom pROC roc auc
## usethis namespace: end
NULL
