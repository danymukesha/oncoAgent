#' @title Utility Functions
#' @description Package-level utility functions and configuration.

#' @name oncoAgent_version
#' @title Get oncoAgent Version
#' @description Returns the current version of the oncoAgent package.
#' @return Version string
#' @export
oncoAgent_version <- function() {
  utils::packageVersion("oncoAgent")
}

#' @name set_oncoagent_config
#' @title Set oncoAgent Configuration
#' @description Sets global configuration options for the oncoAgent package.
#' @param cache_dir Default cache directory
#' @param n_cores Number of parallel cores to use
#' @param verbose Logical, enable verbose logging
#' @param log_level Logging level: "INFO", "DEBUG", "WARN", "ERROR"
#' @return Invisibly returns the configuration list
#' @export
set_oncoagent_config <- function(cache_dir = NULL,
                                 n_cores = 1,
                                 verbose = TRUE,
                                 log_level = "INFO") {

  config <- list(
    cache_dir = cache_dir %||% file.path(tempdir(), "oncoagent_cache"),
    n_cores = n_cores,
    verbose = verbose,
    log_level = log_level
  )

  if (is.null(cache_dir)) {
    dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  options(oncoagent.config = config)

  if (verbose) {
    message("oncoAgent configuration set:")
    message("  Cache dir: ", config$cache_dir)
    message("  Cores: ", config$n_cores)
    message("  Log level: ", config$log_level)
  }

  invisible(config)
}

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  set_oncoagent_config(verbose = FALSE)
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "oncoAgent v", utils::packageVersion("oncoAgent"),
    " - Agentic AI for Oncology Target Discovery"
  )
}

#' @keywords internal
utils::globalVariables(c(
  "IC50", "coexpression_score", "combined_score", "consequence",
  "database_score", "experimental_score", "from", "gene", "pmid",
  "preferredName_A", "preferredName_B", "private", "textmining_score",
  "to", "composite_centrality", "dependency_score", "n",
  "druggability_score", "is_kinase", "has_known_drug"
))
