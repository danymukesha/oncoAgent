#' @title Preprocessing Module
#' @description Quality control, normalization, batch correction, and
#'   data harmonization functions for multi-omics preprocessing.

#' @name normalize_counts
#' @title Normalize Count Data
#' @description Applies normalization to RNA-seq count matrices. Supports
#'   DESeq2 variance-stabilizing transformation (VST), regularized log (rlog),
#'   and edgeR TMM normalization.
#' @param count_matrix Gene x Sample count matrix (integer)
#' @param method Normalization method: "vst", "rlog", "tmm", "cpm"
#' @param col_data Sample metadata data.frame with condition column
#' @param blind Logical, blind VST/rlog to experimental design (default: TRUE)
#' @return Normalized expression matrix
#' @export
normalize_counts <- function(count_matrix,
                             method = "vst",
                             col_data = NULL,
                             blind = TRUE) {

  if (!is.matrix(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }

  storage.mode(count_matrix) <- "integer"

  if (method %in% c("vst", "rlog")) {
    if (is.null(col_data)) {
      col_data <- data.frame(
        condition = rep("A", ncol(count_matrix)),
        row.names = colnames(count_matrix)
      )
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = col_data,
      design = ~ 1
    )

    if (method == "vst") {
      vsd <- DESeq2::vst(dds, blind = blind)
      norm_mat <- SummarizedExperiment::assay(vsd)
    } else {
      rld <- DESeq2::rlog(dds, blind = blind)
      norm_mat <- SummarizedExperiment::assay(rld)
    }

  } else if (method == "tmm") {
    y <- edgeR::DGEList(counts = count_matrix)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    norm_mat <- edgeR::cpm(y, log = TRUE, prior.count = 1)

  } else if (method == "cpm") {
    y <- edgeR::DGEList(counts = count_matrix)
    norm_mat <- edgeR::cpm(y, log = TRUE, prior.count = 1)

  } else {
    stop("Unknown normalization method: ", method)
  }

  message("Normalized ", nrow(norm_mat), " genes x ", ncol(norm_mat),
          " samples using ", toupper(method))
  norm_mat
}

#' @name filter_low_expression
#' @title Filter Low-Expression Genes
#' @description Removes genes with low expression across samples to reduce
#'   noise and improve statistical power.
#' @param count_matrix Gene x Sample count matrix
#' @param min_count Minimum count threshold (default: 10)
#' @param min_samples Minimum proportion of samples that must exceed
#'   min_count (default: 0.1)
#' @return Filtered count matrix
#' @export
filter_low_expression <- function(count_matrix,
                                  min_count = 10,
                                  min_samples = 0.1) {

  n_samples <- ncol(count_matrix)
  threshold <- ceiling(min_samples * n_samples)

  keep <- rowSums(count_matrix >= min_count) >= threshold

  message("Keeping ", sum(keep), " of ", length(keep), " genes ",
          "(expressed in >= ", min_samples * 100, "% of samples ",
          "with count >= ", min_count, ")")

  count_matrix[keep, , drop = FALSE]
}

#' @name batch_correct
#' @title Batch Correction
#' @description Applies batch correction using ComBat (sva) or limma's
#'   removeBatchEffect.
#' @param expr_matrix Normalized expression matrix (genes x samples)
#' @param batch_factor Factor or vector indicating batch membership
#' @param method "combat" or "limma"
#' @param covariates Optional data.frame of biological covariates to preserve
#' @return Batch-corrected expression matrix
#' @export
batch_correct <- function(expr_matrix,
                          batch_factor,
                          method = "combat",
                          covariates = NULL) {

  batch_factor <- as.factor(batch_factor)

  if (method == "combat") {
    if (!requireNamespace("sva", quietly = TRUE)) {
      message("sva package not available, falling back to limma removeBatchEffect")
      method <- "limma"
    }
  }

  if (method == "combat") {
    mod <- NULL
    if (!is.null(covariates)) {
      mod <- stats::model.matrix(~ ., data = covariates)
    }

    corrected <- sva::ComBat(
      dat = expr_matrix,
      batch = batch_factor,
      mod = mod,
      par.prior = TRUE
    )

  } else if (method == "limma") {
    corrected <- limma::removeBatchEffect(
      x = expr_matrix,
      batch = batch_factor
    )

  } else {
    stop("Unknown batch correction method: ", method)
  }

  message("Applied ", toupper(method), " batch correction across ",
          nlevels(batch_factor), " batches")
  corrected
}

#' @name quality_control_metrics
#' @title Compute Quality Control Metrics
#' @description Calculates per-sample QC metrics for RNA-seq data including
#'   library size, detection rate, and outlier flags.
#' @param count_matrix Gene x Sample count matrix
#' @param col_data Sample metadata data.frame
#' @return data.frame with QC metrics per sample
#' @export
quality_control_metrics <- function(count_matrix, col_data = NULL) {

  if (is.null(col_data)) {
    col_data <- data.frame(
      sample_id = colnames(count_matrix),
      row.names = colnames(count_matrix)
    )
  }

  qc <- data.frame(
    sample_id = colnames(count_matrix),
    library_size = colSums(count_matrix),
    n_genes_detected = colSums(count_matrix > 0),
    pct_genes_detected = colSums(count_matrix > 0) / nrow(count_matrix) * 100,
    median_expression = apply(count_matrix, 2, median),
    mean_expression = colMeans(count_matrix),
    stringsAsFactors = FALSE
  )

  qc$log_library_size <- log10(qc$library_size + 1)

  lib_median <- stats::median(qc$library_size)
  lib_mad <- stats::mad(qc$library_size)
  qc$library_size_outlier <- abs(qc$library_size - lib_median) > 3 * lib_mad

  det_median <- stats::median(qc$n_genes_detected)
  det_mad <- stats::mad(qc$n_genes_detected)
  qc$detection_rate_outlier <- abs(qc$n_genes_detected - det_median) > 3 * det_mad

  qc$is_outlier <- qc$library_size_outlier | qc$detection_rate_outlier

  message("QC summary: ", sum(!qc$is_outlier), "/", nrow(qc),
          " samples pass filters")

  qc
}

#' @name harmonize_clinical
#' @title Harmonize Clinical Data
#' @description Standardizes clinical data from TCGA into a consistent format
#'   for downstream survival and subgroup analyses.
#' @param clinical Raw clinical data.frame (from TCGAbiolinks or similar)
#' @param vital_status_col Column name for vital status
#' @param time_col Column name for follow-up time
#' @param time_unit Unit of time: "days", "months", "years"
#' @return Harmonized data.frame with standardized column names
#' @export
harmonize_clinical <- function(clinical,
                               vital_status_col = NULL,
                               time_col = NULL,
                               time_unit = "days") {

  harmonized <- clinical

  if (!is.null(vital_status_col)) {
    harmonized$vital_status <- clinical[[vital_status_col]]
  } else {
    vital_candidates <- c("vital_status", "days_to_death", "Overall.Survival.Status")
    for (col in vital_candidates) {
      if (col %in% names(clinical)) {
        harmonized$vital_status <- clinical[[col]]
        break
      }
    }
  }

  if (!is.null(harmonized$vital_status)) {
    harmonized$status <- as.integer(
      grepl("dead|1|deceased", harmonized$vital_status, ignore.case = TRUE)
    )
  }

  if (!is.null(time_col)) {
    harmonized$time <- as.numeric(clinical[[time_col]])
  } else {
    time_candidates <- c("days_to_last_follow_up", "Overall.Survival..Months.",
                         "last_followup", "days_to_death")
    for (col in time_candidates) {
      if (col %in% names(clinical)) {
        harmonized$time <- as.numeric(clinical[[col]])
        break
      }
    }
  }

  if (!is.null(harmonized$time) && time_unit == "months") {
    harmonized$time <- harmonized$time * 30.44
  } else if (!is.null(harmonized$time) && time_unit == "years") {
    harmonized$time <- harmonized$time * 365.25
  }

  if (is.null(harmonized$sample_id)) {
    if ("bcr_patient_barcode" %in% names(clinical)) {
      harmonized$sample_id <- clinical$bcr_patient_barcode
    } else if (!is.null(rownames(clinical))) {
      harmonized$sample_id <- rownames(clinical)
    }
  }

  harmonized
}
