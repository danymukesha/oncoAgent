#' @title Orchestration Module
#' @description Pipeline orchestration with adaptive agentic loop.
#'   Manages the full workflow from data ingestion through target scoring
#'   with gap-driven re-execution.

#' @name create_pipeline_config
#' @title Create Pipeline Configuration
#' @description Creates a standardized configuration list for the oncoAgent
#'   pipeline.
#' @param project TCGA project identifier
#' @param cancer_type Cancer type label for reference panel
#' @param de_methods DE methods to use
#' @param padj_threshold Adjusted p-value threshold
#' @param max_iterations Maximum agentic loop iterations
#' @param output_dir Directory for results
#' @return Configuration list
#' @export
create_pipeline_config <- function(project = "TCGA-BRCA",
                                   cancer_type = "breast cancer",
                                   de_methods = c("deseq2", "edger", "limma"),
                                   padj_threshold = 0.05,
                                   max_iterations = 5,
                                   output_dir = "results") {

  list(
    project = project,
    cancer_type = cancer_type,
    de_methods = de_methods,
    padj_threshold = padj_threshold,
    lfc_threshold = 0,
    required_modalities = c("expression", "clinical"),
    optional_modalities = c("mutation", "ppi", "crispr", "drug"),
    top_genes_for_network = 200,
    top_genes_for_survival = 100,
    max_iterations = max_iterations,
    output_dir = output_dir,
    cache_dir = file.path(tempdir(), "oncoagent_cache"),
    verbose = TRUE,
    timestamp = Sys.time()
  )
}

#' @name run_oncopipeline
#' @title Run oncoAgent Pipeline
#' @description Executes the full agentic pipeline for oncology target discovery.
#'   This is the main entry point for running the system.
#' @param config Pipeline configuration (from create_pipeline_config)
#' @param existing_data Optional pre-loaded data to skip ingestion
#' @param gold_standard Optional gold standard gene set for evaluation
#' @return List with: targets, gaps, evaluation, agent_states, config
#' @export
run_oncopipeline <- function(config,
                             existing_data = NULL,
                             gold_standard = NULL) {

  message("========================================")
  message("oncoAgent Pipeline: ", config$project)
  message("========================================")

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize agents
  data_agent <- DataDiscoveryAgent$new(
    project = config$project,
    required_modalities = config$required_modalities,
    cache_dir = config$cache_dir
  )

  analysis_agent <- AnalysisAgent$new(
    de_methods = config$de_methods,
    padj_threshold = config$padj_threshold,
    lfc_threshold = config$lfc_threshold
  )

  hypothesis_agent <- HypothesisAgent$new(
    top_n_hypotheses = 20
  )

  gap_agent <- GapDetectionAgent$new()

  eval_agent <- EvaluationAgent$new()

  # Phase 1: Data Discovery
  message("\n--- Phase 1: Data Discovery ---")
  data_state <- data_agent$assess(existing_data %||% list())

  if (length(data_state$missing_modalities) > 0) {
    message("Missing modalities: ", paste(data_state$missing_modalities, collapse = ", "))
    for (rec in data_state$recommendations) {
      message("  Recommendation: ", rec)
    }

    data <- tryCatch(
      data_agent$fetch_missing(existing_data %||% list()),
      error = function(e) {
        message("Data fetch failed: ", e$message)
        existing_data %||% list()
      }
    )
  } else {
    data <- existing_data %||% list()
  }

  # Phase 2: Analysis
  message("\n--- Phase 2: Analysis ---")
  results <- analysis_agent$run(data)

  # Phase 3: Gap Detection
  message("\n--- Phase 3: Gap Detection ---")
  gaps <- gap_agent$detect(results, data)
  gap_summary <- gap_agent$summary()

  message("Gaps detected: ", gap_summary$n_gaps)
  message("Overall severity: ", gap_summary$severity)

  for (gap in gaps) {
    message("  [", toupper(gap$severity), "] ", gap$description)
  }

  # Phase 4: Adaptive Re-analysis (if needed)
  iteration <- 1
  while (gap_agent$needs_rerun() && iteration < config$max_iterations) {
    message("\n--- Re-analysis iteration ", iteration + 1, " ---")

    if (any(sapply(gaps, function(g) grepl("missing", g$type)))) {
      message("Attempting to fill data gaps...")
      data <- tryCatch(
        data_agent$fetch_missing(data),
        error = function(e) { message("Re-fetch failed: ", e$message); data }
      )
    }

    results <- analysis_agent$run(data)
    gaps <- gap_agent$detect(results, data)

    iteration <- iteration + 1
  }

  # Phase 5: Hypothesis Generation
  message("\n--- Phase 5: Hypothesis Generation ---")
  if (!is.null(results$scoring)) {
    hypotheses <- hypothesis_agent$generate(
      scored_targets = results$scoring,
      literature = data$literature,
      network_results = results$network
    )
    message("Generated ", nrow(hypotheses), " hypotheses")
  } else {
    hypotheses <- NULL
  }

  # Phase 6: Target Ranking
  message("\n--- Phase 6: Target Ranking ---")
  ranked_targets <- NULL
  if (!is.null(results$scoring)) {
    ranked_targets <- rank_targets(results$scoring, top_n = 50)
    message("Top 10 targets:")
    print(utils::head(ranked_targets[, c("gene", "composite_score", "target_tier")], 10))
  }

  # Phase 7: Evaluation
  message("\n--- Phase 7: Evaluation ---")
  evaluation <- NULL
  if (!is.null(ranked_targets)) {
    evaluation <- eval_agent$evaluate(ranked_targets, gold_standard)
    message(eval_agent$report())
  }

  # Assemble output
  output <- list(
    targets = ranked_targets,
    hypotheses = hypotheses,
    gaps = gap_summary,
    evaluation = evaluation,
    results = results,
    agent_states = list(
      data_discovery = data_agent$status(),
      analysis = analysis_agent$status(),
      gap_detection = gap_agent$summary()
    ),
    config = config,
    timestamp = Sys.time()
  )

  # Export
  export_results(output, config$output_dir)

  message("\n========================================")
  message("Pipeline complete. Results in: ", config$output_dir)
  message("========================================")

  output
}

#' @name adaptive_pipeline_loop
#' @title Adaptive Pipeline Loop
#' @description Runs the pipeline iteratively, adjusting configuration
#'   based on gap detection results. This is the "agentic iteration" loop.
#' @param config Initial pipeline configuration
#' @param max_iterations Maximum iterations (default: 5)
#' @param improvement_threshold Minimum improvement to continue (default: 0.01)
#' @return List with final results and iteration history
#' @export
adaptive_pipeline_loop <- function(config,
                                   max_iterations = 5,
                                   improvement_threshold = 0.01) {

  history <- list()
  best_result <- NULL
  best_score <- -Inf

  for (i in seq_len(max_iterations)) {
    message("\n=== Adaptive Loop Iteration ", i, "/", max_iterations, " ===")

    result <- tryCatch(
      run_oncopipeline(config),
      error = function(e) {
        message("Pipeline failed on iteration ", i, ": ", e$message)
        NULL
      }
    )

    if (is.null(result)) {
      message("Pipeline failed, stopping loop")
      break
    }

    current_score <- if (!is.null(result$targets)) {
      mean(result$targets$composite_score, na.rm = TRUE)
    } else {
      0
    }

    history[[i]] <- list(
      iteration = i,
      score = current_score,
      n_targets = if (!is.null(result$targets)) nrow(result$targets) else 0,
      n_gaps = result$gaps$n_gaps,
      severity = result$gaps$severity
    )

    if (current_score > best_score + improvement_threshold) {
      best_score <- current_score
      best_result <- result
      message("New best score: ", round(current_score, 4))
    }

    if (result$gaps$severity == "info") {
      message("No critical gaps remaining. Converged at iteration ", i)
      break
    }

    config <- .adjust_config(config, result$gaps)
  }

  list(
    best_result = best_result,
    history = history,
    n_iterations = length(history),
    converged = length(history) < max_iterations
  )
}

#' @keywords internal
.adjust_config <- function(config, gaps) {
  for (gap in gaps$gaps) {
    if (gap$type == "low_deg_signal") {
      config$padj_threshold <- min(config$padj_threshold * 2, 0.25)
      message("Adjusted padj_threshold to ", config$padj_threshold)
    }
    if (gap$type == "small_sample") {
      config$de_methods <- c("limma")
      message("Switched to limma-only for small sample size")
    }
  }
  config
}

#' @name pipeline_status
#' @title Pipeline Status
#' @description Returns current status of a pipeline run
#' @param result Output from run_oncopipeline()
#' @return Formatted status string
#' @export
pipeline_status <- function(result) {
  lines <- c(
    "=== Pipeline Status ===",
    paste("Project:", result$config$project),
    paste("Timestamp:", result$timestamp),
    paste("Targets identified:", if (!is.null(result$targets)) nrow(result$targets) else 0),
    paste("Gaps detected:", result$gaps$n_gaps),
    paste("Overall severity:", result$gaps$severity),
    "",
    "Agent States:"
  )

  for (agent_name in names(result$agent_states)) {
    lines <- c(lines, paste0("  ", agent_name, ": ",
                              result$agent_states[[agent_name]]$n_actions %||% "N/A",
                              " actions"))
  }

  paste(lines, collapse = "\n")
}

#' @name export_results
#' @title Export Pipeline Results
#' @description Exports pipeline results to CSV and JSON formats
#' @param result Output from run_oncopipeline()
#' @param output_dir Output directory
#' @return Invisibly returns file paths
#' @export
export_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- character()

  if (!is.null(result$targets)) {
    target_file <- file.path(output_dir, "ranked_targets.csv")
    utils::write.csv(result$targets, target_file, row.names = FALSE)
    files <- c(files, target_file)
  }

  if (!is.null(result$hypotheses)) {
    hyp_file <- file.path(output_dir, "hypotheses.csv")
    utils::write.csv(result$hypotheses, hyp_file, row.names = FALSE)
    files <- c(files, hyp_file)
  }

  if (!is.null(result$evaluation)) {
    eval_file <- file.path(output_dir, "evaluation.json")
    jsonlite::write_json(result$evaluation, eval_file, pretty = TRUE, auto_unbox = TRUE)
    files <- c(files, eval_file)
  }

  config_file <- file.path(output_dir, "config.json")
  jsonlite::write_json(result$config, config_file, pretty = TRUE, auto_unbox = TRUE)
  files <- c(files, config_file)

  summary_file <- file.path(output_dir, "pipeline_summary.txt")
  writeLines(pipeline_status(result), summary_file)
  files <- c(files, summary_file)

  message("Exported ", length(files), " files to ", output_dir)
  invisible(files)
}
