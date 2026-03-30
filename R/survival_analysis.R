#' @title Survival Analysis Module
#' @description Cox regression, Kaplan-Meier, and multivariate survival
#'   analysis for linking gene expression to patient outcomes.

#' @name run_cox_regression
#' @title Cox Proportional Hazards Regression
#' @description Fits univariate Cox PH model for a single gene's expression
#'   against survival outcomes.
#' @param gene_expression Numeric vector of expression values for one gene
#' @param time Numeric vector of survival times
#' @param status Integer vector of event status (1 = event, 0 = censored)
#' @param gene_name Gene symbol for labeling
#' @param continuous Logical, treat expression as continuous (default: TRUE)
#' @param median_split Logical, dichotomize at median (default: FALSE)
#' @return List with: model, summary, hazard_ratio, pvalue, concordance
#' @export
run_cox_regression <- function(gene_expression,
                               time,
                               status,
                               gene_name = NULL,
                               continuous = TRUE,
                               median_split = FALSE) {

  if (median_split) {
    med <- stats::median(gene_expression, na.rm = TRUE)
    gene_expression <- ifelse(gene_expression > med, "High", "Low")
    gene_expression <- factor(gene_expression, levels = c("Low", "High"))
  }

  df <- data.frame(
    expr = gene_expression,
    time = as.numeric(time),
    status = as.integer(status),
    stringsAsFactors = !continuous && !is.factor(gene_expression)
  )

  complete <- stats::complete.cases(df)
  df <- df[complete, ]

  if (nrow(df) < 10) {
    warning("Too few observations (", nrow(df), ") for Cox regression")
    return(NULL)
  }

  if (length(unique(df$status)) < 2) {
    warning("No events observed, cannot fit Cox model")
    return(NULL)
  }

  fit <- survival::coxph(survival::Surv(time, status) ~ expr, data = df)
  sfit <- summary(fit)

  hr <- exp(stats::coef(fit))
  pval <- sfit$coefficients[, "Pr(>|z|)"]
  ci <- exp(stats::confint(fit))

  result <- list(
    gene = gene_name,
    model = fit,
    summary = sfit,
    n_samples = nrow(df),
    n_events = sum(df$status),
    hazard_ratio = as.numeric(hr),
    hr_ci_lower = as.numeric(ci[1]),
    hr_ci_upper = as.numeric(ci[2]),
    pvalue = as.numeric(pval),
    concordance = sfit$concordance["C"],
    call = match.call()
  )

  result
}

#' @name run_kaplan_meier
#' @title Kaplan-Meier Survival Analysis
#' @description Fits Kaplan-Meier survival curves and performs log-rank test.
#' @param gene_expression Numeric vector of expression values
#' @param time Numeric vector of survival times
#' @param status Integer vector of event status
#' @param gene_name Gene symbol
#' @param split_method Dichotomization method: "median", "quartile", "optimal"
#' @return List with: fit, pvalue, group_survival, plot_data
#' @export
run_kaplan_meier <- function(gene_expression,
                             time,
                             status,
                             gene_name = NULL,
                             split_method = "median") {

  df <- data.frame(
    expr = gene_expression,
    time = as.numeric(time),
    status = as.integer(status)
  )

  df <- df[stats::complete.cases(df), ]

  if (split_method == "median") {
    cutoff <- stats::median(df$expr)
  } else if (split_method == "quartile") {
    q1 <- stats::quantile(df$expr, 0.25)
    q3 <- stats::quantile(df$expr, 0.75)
    df <- df[df$expr <= q1 | df$expr >= q3, ]
    cutoff <- (q1 + q3) / 2
  } else if (split_method == "optimal") {
    cutoff <- .find_optimal_cutoff(df$expr, df$time, df$status)
  } else {
    stop("Unknown split method: ", split_method)
  }

  df$group <- ifelse(df$expr > cutoff, "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))

  km_fit <- survival::survfit(survival::Surv(time, status) ~ group, data = df)

  surv_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = df)
  pval <- 1 - stats::pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

  group_summary <- data.frame(
    group = names(km_fit$strata),
    n = km_fit$n,
    median_survival = km_fit$time[which.min(abs(km_fit$surv - 0.5))],
    stringsAsFactors = FALSE
  )

  plot_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper,
    strata = rep(names(km_fit$strata), km_fit$strata),
    stringsAsFactors = FALSE
  )

  list(
    gene = gene_name,
    fit = km_fit,
    surv_diff = surv_diff,
    pvalue = pval,
    cutoff = cutoff,
    split_method = split_method,
    group_survival = group_summary,
    plot_data = plot_data
  )
}

#' @keywords internal
.find_optimal_cutoff <- function(expr, time, status, min_group = 10) {
  candidates <- unique(stats::quantile(expr, probs = seq(0.25, 0.75, by = 0.05)))

  best_pval <- 1
  best_cut <- stats::median(expr)

  for (cut in candidates) {
    group <- ifelse(expr > cut, "High", "Low")
    if (sum(group == "High") < min_group || sum(group == "Low") < min_group) next

    surv_diff <- survival::survdiff(
      survival::Surv(time, status) ~ group
    )
    pval <- 1 - stats::pchisq(surv_diff$chisq, 1)

    if (pval < best_pval) {
      best_pval <- pval
      best_cut <- cut
    }
  }

  best_cut
}

#' @name survival_forest
#' @title Survival Random Forest
#' @description Runs Cox regression across multiple genes simultaneously using
#'   penalized regression (LASSO) for variable selection.
#' @param expr_matrix Gene x Sample expression matrix
#' @param time Numeric vector of survival times
#' @param status Integer vector of event status
#' @param alpha Elastic net mixing parameter (1 = LASSO, 0 = Ridge)
#' @param n_folds Cross-validation folds (default: 10)
#' @return List with: selected_genes, coefficients, lambda, cv_error
#' @export
survival_forest <- function(expr_matrix,
                            time,
                            status,
                            alpha = 1,
                            n_folds = 10) {

  df <- data.frame(t(expr_matrix), time = time, status = status)
  df <- df[stats::complete.cases(df), ]

  x <- as.matrix(df[, !(names(df) %in% c("time", "status"))])
  y <- survival::Surv(df$time, df$status)

  set.seed(42)
  cv_fit <- glmnet::cv.glmnet(
    x, y,
    family = "cox",
    alpha = alpha,
    nfolds = n_folds,
    type.measure = "C"
  )

  coefs <- stats::coef(cv_fit, s = "lambda.min")
  selected <- which(coefs[, 1] != 0)
  selected_genes <- rownames(coefs)[selected]
  selected_coefs <- coefs[selected, 1]

  message("Survival LASSO: ", length(selected_genes),
          " genes selected at lambda.min = ", round(cv_fit$lambda.min, 4))

  list(
    selected_genes = selected_genes,
    coefficients = data.frame(
      gene = selected_genes,
      coefficient = selected_coefs,
      hazard_ratio = exp(selected_coefs),
      stringsAsFactors = FALSE
    ),
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cv_error = min(cv_fit$cvm),
    model = cv_fit
  )
}

#' @name multivariate_survival
#' @title Multivariate Survival Analysis
#' @description Fits a multivariate Cox model with gene expression and
#'   clinical covariates.
#' @param expr_matrix Gene x Sample expression matrix (top genes)
#' @param clinical_df data.frame with time, status, and clinical covariates
#' @param clinical_vars Character vector of clinical variable names to include
#' @return List with: model, coefficients, forest_plot_data
#' @export
multivariate_survival <- function(expr_matrix,
                                  clinical_df,
                                  clinical_vars = NULL) {

  merged <- merge(
    data.frame(t(expr_matrix), check.names = FALSE),
    clinical_df,
    by = 0
  )
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL

  gene_vars <- intersect(colnames(expr_matrix), colnames(merged))
  gene_vars <- gene_vars[1:min(20, length(gene_vars))]

  all_vars <- c(gene_vars, intersect(clinical_vars, colnames(merged)))
  all_vars <- intersect(all_vars, colnames(merged))
  all_vars <- all_vars[all_vars != "time" & all_vars != "status"]

  formula_str <- paste("survival::Surv(time, status) ~",
                       paste(all_vars, collapse = " + "))
  formula_obj <- stats::as.formula(formula_str)

  fit <- survival::coxph(formula_obj, data = merged)
  sfit <- summary(fit)

  coefs <- data.frame(
    variable = rownames(sfit$coefficients),
    coefficient = sfit$coefficients[, "coef"],
    hr = exp(sfit$coefficients[, "coef"]),
    hr_lower = exp(sfit$conf.int[, "lower .95"]),
    hr_upper = exp(sfit$conf.int[, "upper .95"]),
    pvalue = sfit$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  coefs$significance <- ifelse(coefs$pvalue < 0.001, "***",
                        ifelse(coefs$pvalue < 0.01, "**",
                        ifelse(coefs$pvalue < 0.05, "*", "ns")))

  list(
    model = fit,
    summary = sfit,
    coefficients = coefs,
    concordance = sfit$concordance["C"],
    formula = formula_str
  )
}

#' @name calculate_hazard_ratios
#' @title Calculate Hazard Ratios for Multiple Genes
#' @description Runs univariate Cox regression for each gene and returns
#'   a summary table of hazard ratios.
#' @param expr_matrix Gene x Sample expression matrix
#' @param time Numeric vector of survival times
#' @param status Integer vector of event status
#' @param genes Character vector of genes to test (NULL = all)
#' @param n_cores Number of parallel cores (default: 1)
#' @return data.frame with hazard ratios, CIs, and p-values per gene
#' @export
calculate_hazard_ratios <- function(expr_matrix,
                                    time,
                                    status,
                                    genes = NULL,
                                    n_cores = 1) {

  if (is.null(genes)) {
    genes <- rownames(expr_matrix)
  }

  genes <- intersect(genes, rownames(expr_matrix))

  run_one <- function(gene) {
    res <- tryCatch(
      run_cox_regression(
        expr_matrix[gene, ], time, status, gene_name = gene
      ),
      error = function(e) NULL
    )

    if (is.null(res)) return(NULL)

    data.frame(
      gene = gene,
      hr = res$hazard_ratio,
      hr_lower = res$hr_ci_lower,
      hr_upper = res$hr_ci_upper,
      pvalue = res$pvalue,
      concordance = res$concordance,
      n_samples = res$n_samples,
      n_events = res$n_events,
      stringsAsFactors = FALSE
    )
  }

  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(genes, run_one, mc.cores = n_cores)
  } else {
    results <- lapply(genes, run_one)
  }

  hr_table <- dplyr::bind_rows(purrr::compact(results))

  if (nrow(hr_table) > 0) {
    hr_table$padj <- stats::p.adjust(hr_table$pvalue, method = "BH")
    hr_table <- hr_table[order(hr_table$pvalue), ]
  }

  message("Computed HRs for ", nrow(hr_table), "/", length(genes), " genes")
  hr_table
}
