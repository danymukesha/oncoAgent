## Agentic Layer
# R6 classes implementing the five core agent types:
#   DataDiscoveryAgent, AnalysisAgent, HypothesisAgent,
#   GapDetectionAgent, and EvaluationAgent.
# Each agent has decision logic, state tracking, and action capabilities.

# ============================================================================
# DataDiscoveryAgent
# ============================================================================

#' @title Data Discovery Agent
#' @description R6 class that checks data availability, identifies missing
#'   modalities, and recommends data sources for a given cancer project.
#' @export
DataDiscoveryAgent <- R6::R6Class(
    "DataDiscoveryAgent",
    public = list(
        #' @field state Current agent state
        state = NULL,

        #' @description Create a new DataDiscoveryAgent
        #' @param project TCGA project identifier
        #' @param required_modalities Character vector of required data types
        #' @param cache_dir Cache directory for data
        initialize = function(project = "TCGA-BRCA",
                              required_modalities = c(
                                  "expression", "clinical",
                                  "mutation", "ppi"
                              ),
                              cache_dir = NULL) {
            self$state <- list(
                project = project,
                required_modalities = required_modalities,
                available_modalities = character(),
                missing_modalities = character(),
                data_quality = list(),
                recommendations = character(),
                action_log = list(),
                cache_dir = cache_dir %||% file.path(tempdir(), "oncoagent_cache")
            )
        },

        #' @description Assess what data is available and what is missing
        #' @param existing_data Named list of already-fetched data
        #' @return Updated state with availability assessment
        assess = function(existing_data = list()) {
            available <- names(existing_data)
            missing <- setdiff(self$state$required_modalities, available)

            self$state$available_modalities <- available
            self$state$missing_modalities <- missing

            for (mod in available) {
                self$state$data_quality[[mod]] <- private$.assess_quality(
                    existing_data[[mod]], mod
                )
            }

            self$state$recommendations <- private$.generate_recommendations(missing)

            private$.log_action("assess", paste(
                "Found", length(available), "of",
                length(self$state$required_modalities),
                "modalities"
            ))

            self$state
        },

        #' @description Fetch missing data automatically
        #' @param existing_data Named list of current data
        #' @return Named list with all fetched data
        fetch_missing = function(existing_data = list()) {
            missing <- self$state$missing_modalities

            for (mod in missing) {
                message("DataDiscoveryAgent: Fetching missing modality: ", mod)
                existing_data[[mod]] <- private$.fetch_modality(mod)
                private$.log_action("fetch", paste("Retrieved", mod))
            }

            self$assess(existing_data)
            existing_data
        },

        #' @description Get agent status report
        #' @return List with current status
        status = function() {
            list(
                project = self$state$project,
                available = self$state$available_modalities,
                missing = self$state$missing_modalities,
                quality = self$state$data_quality,
                recommendations = self$state$recommendations,
                n_actions = length(self$state$action_log)
            )
        }
    ),
    private = list(
        .assess_quality = function(data, modality) {
            quality <- list(completeness = 0, n_samples = 0, n_features = 0)

            if (modality == "expression" && methods::is(data, "SummarizedExperiment")) {
                quality$n_samples <- ncol(data)
                quality$n_features <- nrow(data)
                quality$completeness <- 1 - sum(is.na(SummarizedExperiment::assay(data))) /
                    length(SummarizedExperiment::assay(data))
            } else if (is.data.frame(data)) {
                quality$n_samples <- nrow(data)
                quality$n_features <- ncol(data)
                quality$completeness <- 1 - sum(is.na(data)) / length(data)
            } else if (is.matrix(data)) {
                quality$n_samples <- ncol(data)
                quality$n_features <- nrow(data)
                quality$completeness <- 1 - sum(is.na(data)) / length(data)
            }

            quality
        },
        .generate_recommendations = function(missing) {
            recs <- character()
            if ("expression" %in% missing) {
                recs <- c(recs, "Fetch RNA-seq counts from TCGA via TCGAbiolinks::GDCquery")
            }
            if ("clinical" %in% missing) {
                recs <- c(recs, "Fetch clinical annotations from TCGA via GDCquery_clinic")
            }
            if ("mutation" %in% missing) {
                recs <- c(recs, "Fetch somatic mutations from TCGA (Masked Somatic Mutation)")
            }
            if ("ppi" %in% missing) {
                recs <- c(recs, "Fetch PPI network from STRING-db for DEG list")
            }
            if ("crispr" %in% missing) {
                recs <- c(recs, "Fetch CRISPR dependency scores from DepMap portal")
            }
            if ("drug" %in% missing) {
                recs <- c(recs, "Fetch drug sensitivity data from GDSC")
            }
            recs
        },
        .fetch_modality = function(modality) {
            switch(modality,
                expression = tryCatch(
                    fetch_tcga_data(self$state$project, "expression")$expression,
                    error = function(e) {
                        message("Failed: ", e$message)
                        NULL
                    }
                ),
                clinical = tryCatch(
                    fetch_tcga_data(self$state$project, "clinical")$clinical,
                    error = function(e) {
                        message("Failed: ", e$message)
                        NULL
                    }
                ),
                mutation = tryCatch(
                    fetch_tcga_data(self$state$project, "mutation")$mutation,
                    error = function(e) {
                        message("Failed: ", e$message)
                        NULL
                    }
                ),
                crispr = tryCatch(
                    fetch_depmap_crispr(),
                    error = function(e) {
                        message("Failed: ", e$message)
                        NULL
                    }
                ),
                drug = tryCatch(
                    fetch_gdsc_drug_sensitivity(),
                    error = function(e) {
                        message("Failed: ", e$message)
                        NULL
                    }
                ),
                {
                    message("Unknown modality: ", modality)
                    NULL
                }
            )
        },
        .log_action = function(action, detail) {
            entry <- list(
                action = action,
                detail = detail,
                timestamp = Sys.time()
            )
            self$state$action_log[[length(self$state$action_log) + 1]] <- entry
        }
    )
)

# ============================================================================
# AnalysisAgent
# ============================================================================

#' @title Analysis Agent
#' @description R6 class that executes the analytical pipeline (DE, survival,
#'   network analysis) and manages analysis state.
#'
#' @field state A list storing the internal state of the agent, including:
#'   \itemize{
#'     \item \code{de_methods}: Character vector of differential expression methods.
#'     \item \code{padj_threshold}: Adjusted p-value threshold.
#'     \item \code{lfc_threshold}: Log-fold-change threshold.
#'     \item \code{results}: List of results from different analysis steps.
#'     \item \code{execution_log}: List of log messages with timestamps.
#'     \item \code{completed_steps}: Character vector of completed steps.
#'   }
#'
#' @export
AnalysisAgent <- R6::R6Class(
    "AnalysisAgent",
    public = list(
        state = NULL,

        #' @description Create a new AnalysisAgent
        #' @param de_methods Character vector of DE methods
        #' @param padj_threshold Adjusted p-value threshold
        #' @param lfc_threshold Log-fold-change threshold
        initialize = function(de_methods = c("deseq2", "edger", "limma"),
                              padj_threshold = 0.05,
                              lfc_threshold = 0) {
            self$state <- list(
                de_methods = de_methods,
                padj_threshold = padj_threshold,
                lfc_threshold = lfc_threshold,
                results = list(),
                execution_log = list(),
                completed_steps = character()
            )
        },

        #' @description Run full analysis pipeline
        #' @param data Named list of data (expression, clinical, ppi, etc.)
        #' @return Named list of all analysis results
        run = function(data) {
            results <- list()

            if (!is.null(data$expression)) {
                results$de <- private$.run_de(data$expression, data$col_data)
                self$state$completed_steps <- c(self$state$completed_steps, "de")
            }

            if (!is.null(results$de) && !is.null(data$clinical)) {
                results$survival <- private$.run_survival(results$de, data)
                self$state$completed_steps <- c(self$state$completed_steps, "survival")
            }

            if (!is.null(results$de) && !is.null(data$ppi)) {
                results$network <- private$.run_network(results$de, data$ppi)
                self$state$completed_steps <- c(self$state$completed_steps, "network")
            }

            if (!is.null(results$de)) {
                results$scoring <- private$.run_scoring(results, data)
                self$state$completed_steps <- c(self$state$completed_steps, "scoring")
            }

            self$state$results <- results
            results
        },

        #' @description Get analysis status
        #' @return List with analysis status
        status = function() {
            list(
                completed_steps = self$state$completed_steps,
                methods = self$state$de_methods,
                thresholds = list(
                    padj = self$state$padj_threshold,
                    lfc = self$state$lfc_threshold
                ),
                n_results = length(self$state$results)
            )
        }
    ),
    private = list(
        .run_de = function(expression, col_data) {
            private$.log("Starting differential expression analysis")

            expr_mat <- if (methods::is(expression, "SummarizedExperiment")) {
                SummarizedExperiment::assay(expression)
            } else {
                as.matrix(expression)
            }

            if (is.null(col_data)) {
                col_data <- data.frame(
                    condition = ifelse(grepl("Tumor|Cancer", colnames(expr_mat), ignore.case = TRUE),
                        "tumor", "normal"
                    ),
                    row.names = colnames(expr_mat)
                )
            }
            col_data$condition <- as.factor(col_data$condition)

            result <- tryCatch(
                consensus_de(expr_mat, col_data,
                    methods = self$state$de_methods,
                    padj_threshold = self$state$padj_threshold
                ),
                error = function(e) {
                    message("Consensus DE failed, trying single method: ", e$message)
                    tryCatch(
                        run_deseq2(expr_mat, col_data,
                            padj_threshold = self$state$padj_threshold
                        ),
                        error = function(e2) {
                            message("All DE failed")
                            NULL
                        }
                    )
                }
            )

            result
        },
        .run_survival = function(de_results, data) {
            private$.log("Starting survival analysis")

            clin <- harmonize_clinical(data$clinical)

            if (is.null(clin$time) || is.null(clin$status)) {
                message("No survival data available in clinical data")
                return(NULL)
            }

            expr_mat <- if (methods::is(data$expression, "SummarizedExperiment")) {
                SummarizedExperiment::assay(data$expression)
            } else {
                as.matrix(data$expression)
            }

            top_genes <- if (!is.null(de_results$gene)) {
                head(de_results$gene, 100)
            } else {
                head(rownames(expr_mat), 100)
            }

            common_samples <- intersect(colnames(expr_mat), clin$sample_id)
            if (length(common_samples) < 20) {
                message(
                    "Too few overlapping samples for survival analysis: ",
                    length(common_samples)
                )
                return(NULL)
            }

            expr_sub <- expr_mat[top_genes, common_samples, drop = FALSE]
            clin_sub <- clin[match(common_samples, clin$sample_id), ]

            calculate_hazard_ratios(
                expr_sub,
                time = clin_sub$time,
                status = clin_sub$status,
                genes = top_genes
            )
        },
        .run_network = function(de_results, ppi_data) {
            private$.log("Starting network analysis")

            top_genes <- if (!is.null(de_results$gene)) {
                head(de_results$gene, 200)
            } else {
                character()
            }

            if (length(top_genes) < 5) {
                message("Too few genes for network analysis")
                return(NULL)
            }

            if (is.data.frame(ppi_data) && nrow(ppi_data) > 0) {
                g <- build_ppi_network(ppi_data, score_threshold = 400)
            } else {
                synthetic_ppi <- .generate_synthetic_ppi(top_genes)
                g <- build_ppi_network(synthetic_ppi, score_threshold = 400)
            }

            if (igraph::vcount(g) < 3) {
                return(NULL)
            }

            centrality <- calculate_network_centrality(g)
            modules <- detect_network_modules(g)

            list(
                network = g,
                centrality = centrality,
                modules = modules
            )
        },
        .run_scoring = function(results, data) {
            private$.log("Starting target scoring")

            crispr_scores <- NULL
            if (!is.null(data$crispr)) {
                crispr_df <- data.frame(
                    gene = colnames(data$crispr),
                    dependency_score = colMeans(data$crispr, na.rm = TRUE),
                    stringsAsFactors = FALSE
                )
                crispr_scores <- crispr_df
            }

            composite_target_score(
                de_results = results$de,
                survival_hr = results$survival,
                centrality = results$network$centrality,
                crispr_scores = crispr_scores
            )
        },
        .log = function(msg) {
            entry <- list(message = msg, timestamp = Sys.time())
            self$state$execution_log[[length(self$state$execution_log) + 1]] <- entry
            message("AnalysisAgent: ", msg)
        }
    )
)

# ============================================================================
# HypothesisAgent
# ============================================================================

#' @title Hypothesis Agent
#' @description R6 class that generates and ranks hypotheses about target
#'   mechanisms, combining multi-omics evidence with literature signals.
#'
#' @field state A list storing the internal state of the agent, including:
#'   \itemize{
#'     \item \code{top_n}: Maximum number of hypotheses to generate.
#'     \item \code{min_score}: Minimum evidence score threshold.
#'     \item \code{hypotheses}: List of generated hypotheses.
#'     \item \code{evidence_log}: List of evidence collected for each hypothesis.
#'   }
#'
#' @export
HypothesisAgent <- R6::R6Class(
    "HypothesisAgent",
    public = list(
        state = NULL,

        #' @description Create a new HypothesisAgent
        #' @param top_n_hypotheses Maximum hypotheses to generate
        #' @param min_evidence_score Minimum evidence score threshold
        initialize = function(top_n_hypotheses = 20,
                              min_evidence_score = 0.1) {
            self$state <- list(
                top_n = top_n_hypotheses,
                min_score = min_evidence_score,
                hypotheses = list(),
                evidence_log = list()
            )
        },

        #' @description Generate hypotheses from analysis results
        #' @param scored_targets Output from composite_target_score
        #' @param literature Optional literature data.frame
        #' @param network_results Optional network analysis results
        #' @return data.frame of ranked hypotheses
        generate = function(scored_targets, literature = NULL,
                            network_results = NULL) {
            hypotheses <- list()

            for (i in seq_len(min(self$state$top_n, nrow(scored_targets)))) {
                gene <- scored_targets$gene[i]
                score <- scored_targets$composite_score[i]

                hyp <- list(
                    gene = gene,
                    rank = i,
                    composite_score = score,
                    evidence = private$.collect_evidence(
                        gene, scored_targets,
                        literature, network_results
                    ),
                    mechanism_hypothesis = private$.infer_mechanism(
                        gene, scored_targets,
                        network_results
                    ),
                    confidence = private$.estimate_confidence(gene, scored_targets)
                )

                hypotheses[[length(hypotheses) + 1]] <- hyp
            }

            hyp_df <- purrr::map_dfr(hypotheses, function(h) {
                data.frame(
                    gene = h$gene,
                    rank = h$rank,
                    composite_score = h$composite_score,
                    mechanism = h$mechanism_hypothesis,
                    confidence = h$confidence,
                    n_evidence_types = length(h$evidence$types),
                    evidence_summary = paste(h$evidence$types, collapse = ";"),
                    stringsAsFactors = FALSE
                )
            })

            self$state$hypotheses <- hypotheses
            message("HypothesisAgent: Generated ", nrow(hyp_df), " hypotheses")
            hyp_df
        },

        #' @description Get top hypotheses
        #' @param n Number to return
        #' @return data.frame of top hypotheses
        top = function(n = 10) {
            if (length(self$state$hypotheses) == 0) {
                return(NULL)
            }

            hyp_df <- purrr::map_dfr(self$state$hypotheses, function(h) {
                data.frame(
                    gene = h$gene,
                    rank = h$rank,
                    mechanism = h$mechanism_hypothesis,
                    confidence = h$confidence,
                    stringsAsFactors = FALSE
                )
            })

            utils::head(hyp_df, n)
        }
    ),
    private = list(
        .collect_evidence = function(gene, scored_targets, literature, network) {
            evidence <- list(types = character(), details = list())

            row <- scored_targets[scored_targets$gene == gene, ]

            if ("de_score" %in% names(row) && row$de_score > 0) {
                evidence$types <- c(evidence$types, "differential_expression")
                evidence$details$de_score <- row$de_score
            }

            if ("hr_score" %in% names(row) && row$hr_score > 0) {
                evidence$types <- c(evidence$types, "survival_association")
                evidence$details$survival_score <- row$hr_score
            }

            if ("centrality_score" %in% names(row) && row$centrality_score > 0) {
                evidence$types <- c(evidence$types, "network_hub")
                evidence$details$centrality <- row$centrality_score
            }

            if ("crispr_score" %in% names(row) && row$crispr_score > 0) {
                evidence$types <- c(evidence$types, "crispr_dependency")
                evidence$details$dependency <- row$crispr_score
            }

            if (!is.null(literature) && gene %in% literature$gene) {
                evidence$types <- c(evidence$types, "literature_support")
                evidence$details$n_papers <- sum(literature$gene == gene)
            }

            if (!is.null(network) && gene %in% network$centrality$gene) {
                evidence$types <- c(evidence$types, "network_connectivity")
            }

            evidence
        },
        .infer_mechanism = function(gene, scored_targets, network) {
            row <- scored_targets[scored_targets$gene == gene, ]

            mechanisms <- character()

            if ("de_score" %in% names(row) && row$de_score > 0.5) {
                mechanisms <- c(mechanisms, "altered transcription regulation")
            }

            if ("hr_score" %in% names(row) && row$hr_score > 0.5) {
                mechanisms <- c(mechanisms, "prognostic biomarker / survival driver")
            }

            if ("centrality_score" %in% names(row) && row$centrality_score > 0.5) {
                mechanisms <- c(mechanisms, "network hub / pathway bottleneck")
            }

            if ("crispr_score" %in% names(row) && row$crispr_score > 0.5) {
                mechanisms <- c(mechanisms, "essential dependency / synthetic lethality")
            }

            if (length(mechanisms) == 0) {
                mechanisms <- "insufficient evidence for mechanism inference"
            }

            paste(mechanisms, collapse = "; ")
        },
        .estimate_confidence = function(gene, scored_targets) {
            row <- scored_targets[scored_targets$gene == gene, ]

            score_cols <- c("de_score", "hr_score", "centrality_score", "crispr_score")
            available <- intersect(score_cols, names(row))

            if (length(available) == 0) {
                return(0)
            }

            n_active <- sum(row[, available, drop = FALSE] > 0.1)
            n_active / length(available)
        }
    )
)

# ============================================================================
# GapDetectionAgent
# ============================================================================

#' @title Gap Detection Agent
#' @description R6 class that identifies data gaps, statistical weaknesses,
#'   and conflicting signals across evidence dimensions.
#'
#' @field state A list storing the internal state of the agent, including:
#'   \itemize{
#'     \item \code{thresholds}: List of minimum thresholds for analysis (samples, DEGs, modules)
#'     \item \code{gaps}: List of detected gaps with type, severity, and description
#'     \item \code{recommendations}: Character vector of recommendations based on gaps
#'     \item \code{severity}: Overall severity of detected gaps (critical/warning/info)
#'   }
#'
#' @export
GapDetectionAgent <- R6::R6Class(
    "GapDetectionAgent",
    public = list(
        state = NULL,

        #' @description Create a new GapDetectionAgent
        #' @param min_samples Minimum sample size for adequate power
        #' @param min_degs Minimum DEGs expected
        #' @param min_modules Minimum network modules expected
        initialize = function(min_samples = 50,
                              min_degs = 100,
                              min_modules = 3) {
            self$state <- list(
                thresholds = list(
                    min_samples = min_samples,
                    min_degs = min_degs,
                    min_modules = min_modules
                ),
                gaps = list(),
                recommendations = character(),
                severity = character()
            )
        },

        #' @description Detect gaps in analysis results
        #' @param results Named list of pipeline results
        #' @param data Named list of input data
        #' @return List of detected gaps with severity and recommendations
        detect = function(results, data) {
            gaps <- list()

            data_gaps <- private$.check_data_gaps(data)
            analysis_gaps <- private$.check_analysis_gaps(results)
            power_gaps <- private$.check_statistical_power(results, data)
            conflict_gaps <- private$.check_conflicting_signals(results)

            all_gaps <- c(data_gaps, analysis_gaps, power_gaps, conflict_gaps)

            self$state$gaps <- all_gaps
            self$state$recommendations <- private$.generate_gap_recommendations(all_gaps)
            self$state$severity <- private$.assess_overall_severity(all_gaps)

            all_gaps
        },

        #' @description Get gap summary
        #' @return List with gap statistics
        summary = function() {
            gap_types <- if (length(self$state$gaps) > 0) {
                table(sapply(self$state$gaps, `[[`, "type"))
            } else {
                table(character(0))
            }
            list(
                n_gaps = length(self$state$gaps),
                by_type = gap_types,
                severity = self$state$severity,
                recommendations = self$state$recommendations
            )
        },

        #' @description Check if pipeline should be re-run
        #' @return Logical, TRUE if critical gaps exist
        needs_rerun = function() {
            any(sapply(self$state$gaps, `[[`, "severity") == "critical")
        }
    ),
    private = list(
        .check_data_gaps = function(data) {
            gaps <- list()

            required <- c("expression", "clinical")
            for (req in required) {
                if (is.null(data[[req]])) {
                    gaps[[length(gaps) + 1]] <- list(
                        type = "missing_data",
                        severity = "critical",
                        description = paste("Missing required data:", req),
                        recommendation = paste("Fetch", req, "data before proceeding")
                    )
                }
            }

            if (!is.null(data$expression)) {
                expr <- if (methods::is(data$expression, "SummarizedExperiment")) {
                    SummarizedExperiment::assay(data$expression)
                } else {
                    as.matrix(data$expression)
                }

                if (ncol(expr) < self$state$thresholds$min_samples) {
                    gaps[[length(gaps) + 1]] <- list(
                        type = "small_sample",
                        severity = "warning",
                        description = paste(
                            "Only", ncol(expr), "samples, minimum recommended:",
                            self$state$thresholds$min_samples
                        ),
                        recommendation = "Increase sample size or use regularized methods"
                    )
                }

                missing_pct <- sum(is.na(expr)) / length(expr)
                if (missing_pct > 0.1) {
                    gaps[[length(gaps) + 1]] <- list(
                        type = "missing_values",
                        severity = if (missing_pct > 0.3) "critical" else "warning",
                        description = paste0(
                            round(missing_pct * 100, 1),
                            "% missing values in expression data"
                        ),
                        recommendation = "Apply imputation or filter genes with excessive missingness"
                    )
                }
            }

            if (is.null(data$clinical)) {
                gaps[[length(gaps) + 1]] <- list(
                    type = "missing_clinical",
                    severity = "critical",
                    description = "No clinical data available for survival analysis",
                    recommendation = "Fetch clinical annotations to enable survival analysis"
                )
            }

            if (is.null(data$ppi)) {
                gaps[[length(gaps) + 1]] <- list(
                    type = "missing_ppi",
                    severity = "warning",
                    description = "No PPI data available for network analysis",
                    recommendation = "Fetch PPI data from STRING-db after identifying DEGs"
                )
            }

            gaps
        },
        .check_analysis_gaps = function(results) {
            gaps <- list()

            if (is.null(results$de)) {
                gaps[[length(gaps) + 1]] <- list(
                    type = "no_de_results",
                    severity = "critical",
                    description = "Differential expression analysis produced no results",
                    recommendation = "Check data quality and DE parameters"
                )
            } else {
                n_sig <- if (!is.null(results$de$n_methods)) {
                    sum(results$de$n_methods >= 2)
                } else if (!is.null(results$de$padj)) {
                    sum(results$de$padj < 0.05, na.rm = TRUE)
                } else {
                    nrow(results$de)
                }

                if (n_sig < self$state$thresholds$min_degs) {
                    gaps[[length(gaps) + 1]] <- list(
                        type = "low_deg_signal",
                        severity = "warning",
                        description = paste(
                            "Only", n_sig, "significant DEGs detected (threshold:",
                            self$state$thresholds$min_degs, ")"
                        ),
                        recommendation = "Consider relaxing thresholds or checking data normalization"
                    )
                }
            }

            if (is.null(results$survival)) {
                gaps[[length(gaps) + 1]] <- list(
                    type = "no_survival",
                    severity = "warning",
                    description = "Survival analysis not performed",
                    recommendation = "Ensure clinical data with survival endpoints is available"
                )
            }

            if (is.null(results$network)) {
                gaps[[length(gaps) + 1]] <- list(
                    type = "no_network",
                    severity = "info",
                    description = "Network analysis not performed",
                    recommendation = "Fetch PPI data and rerun network analysis"
                )
            }

            gaps
        },
        .check_statistical_power = function(results, data) {
            gaps <- list()

            if (!is.null(results$de) && !is.null(data$expression)) {
                expr <- if (methods::is(data$expression, "SummarizedExperiment")) {
                    SummarizedExperiment::assay(data$expression)
                } else {
                    as.matrix(data$expression)
                }

                n_samples <- ncol(expr)

                if (n_samples < 30) {
                    gaps[[length(gaps) + 1]] <- list(
                        type = "low_power",
                        severity = "critical",
                        description = paste(
                            "Very low sample size (n =", n_samples,
                            "). Results may be unreliable."
                        ),
                        recommendation = "Increase sample size or use penalized/shrinkage methods"
                    )
                }
            }

            gaps
        },
        .check_conflicting_signals = function(results) {
            gaps <- list()

            if (!is.null(results$scoring) && !is.null(results$de)) {
                top_scored <- utils::head(results$scoring$gene, 20)

                if ("n_methods" %in% names(results$de)) {
                    top_de <- results$de$gene[results$de$n_methods >= 2]
                    overlap <- intersect(top_scored, top_de)

                    if (length(overlap) < length(top_scored) * 0.3) {
                        gaps[[length(gaps) + 1]] <- list(
                            type = "conflicting_signals",
                            severity = "warning",
                            description = paste(
                                "Low overlap between top-scored targets and",
                                "consensus DEGs. Evidence may be conflicting."
                            ),
                            recommendation = "Review individual DE method results and survival data"
                        )
                    }
                }
            }

            gaps
        },
        .generate_gap_recommendations = function(gaps) {
            unique(sapply(gaps, `[[`, "recommendation"))
        },
        .assess_overall_severity = function(gaps) {
            severities <- sapply(gaps, `[[`, "severity")
            if ("critical" %in% severities) {
                return("critical")
            }
            if ("warning" %in% severities) {
                return("warning")
            }
            return("info")
        }
    )
)

# ============================================================================
# EvaluationAgent
# ============================================================================

#' @title Evaluation Agent
#' @description R6 class that benchmarks predicted targets against gold
#'   standard gene sets and computes evaluation metrics.
#'
#' @field state A list storing the internal state of the agent, including:
#'   \itemize{
#'     \item \code{thresholds}: List of minimum thresholds (samples, DEGs, modules).
#'     \item \code{gaps}: List of detected gaps with type, severity, and description.
#'     \item \code{recommendations}: Character vector of recommendations based on gaps.
#'     \item \code{severity}: Overall severity of detected gaps (critical/warning/info).
#'   }
#'
#' @export
EvaluationAgent <- R6::R6Class(
    "EvaluationAgent",
    public = list(
        state = NULL,

        #' @description Create a new EvaluationAgent
        #' @param gold_standard Named list of gene sets for benchmarking
        #' @param k_values K values for Precision@K (default: c(10, 20, 50))
        initialize = function(gold_standard = NULL,
                              k_values = c(10, 20, 50)) {
            self$state <- list(
                gold_standard = gold_standard,
                k_values = k_values,
                metrics = list(),
                evaluation_log = list()
            )
        },

        #' @description Evaluate targets against gold standards
        #' @param scored_targets Output from composite_target_score or rank_targets
        #' @param gold_standard Named list of gene sets (if NULL, uses constructor value)
        #' @return List of evaluation metrics
        evaluate = function(scored_targets, gold_standard = NULL) {
            gs <- gold_standard %||% self$state$gold_standard

            if (is.null(gs)) {
                message("EvaluationAgent: No gold standard provided, building reference panel")
                gs <- build_reference_panel()
                self$state$gold_standard <- gs
            }

            predicted_genes <- scored_targets$gene

            metrics <- list()

            if (!is.null(gs$genes)) {
                metrics$overall <- evaluate_targets(predicted_genes, gs$genes)
                metrics$auroc <- calculate_auroc(scored_targets, gs$genes)

                for (k in self$state$k_values) {
                    metrics[[paste0("precision_at_", k)]] <- precision_at_k(
                        scored_targets, gs$genes, k
                    )
                }

                metrics$recall <- recall_of_known_targets(predicted_genes, gs$genes)
            }

            if (!is.null(gs$per_source)) {
                metrics$per_source <- lapply(names(gs$per_source), function(src) {
                    known <- gs$per_source[[src]]
                    list(
                        precision = evaluate_targets(predicted_genes, known),
                        recall = recall_of_known_targets(predicted_genes, known),
                        n_known = length(known),
                        n_recovered = length(intersect(predicted_genes, known))
                    )
                })
                names(metrics$per_source) <- names(gs$per_source)
            }

            self$state$metrics <- metrics
            metrics
        },

        #' @description Generate evaluation report
        #' @return Formatted evaluation report
        report = function() {
            if (length(self$state$metrics) == 0) {
                return("No evaluation has been performed yet. Run $evaluate() first.")
            }

            m <- self$state$metrics

            report_lines <- c(
                "=== oncoAgent Target Evaluation Report ===",
                "",
                paste("Overall Precision:", round(m$overall, 3)),
                paste("Recall:", round(m$recall, 3)),
                "",
                "Precision@K:"
            )

            for (k in self$state$k_values) {
                key <- paste0("precision_at_", k)
                if (!is.null(m[[key]])) {
                    report_lines <- c(
                        report_lines,
                        paste0("  K=", k, ": ", round(m[[key]], 3))
                    )
                }
            }

            if (!is.null(m$auroc)) {
                report_lines <- c(
                    report_lines, "",
                    paste("AUROC:", round(m$auroc, 3))
                )
            }

            paste(report_lines, collapse = "\n")
        }
    )
)
