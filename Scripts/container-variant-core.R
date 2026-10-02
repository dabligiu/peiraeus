# Shared implementation for the two container-dashboard sensitivity variants.
#
# This file is sourced by:
#   container-dashboard-base-2016q2.R
#   container-dashboard-average-pre.R
#
# The calling script must define `variant_config`. Each wrapper is therefore a
# complete, reproducible entry point while the common estimation code remains
# in one auditable location.

suppressPackageStartupMessages({
  library(fixest)
  library(readxl)
  library(jsonlite)
})

if (!exists("variant_config")) stop("The wrapper must define variant_config before sourcing this file")
required_config <- c("id", "title", "output_html", "normalization", "drop_2016_q3")
if (!all(required_config %in% names(variant_config))) stop("variant_config is incomplete")
if (!variant_config$normalization %in% c("q2", "average_pre")) stop("Unknown normalization")

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
wrapper_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(wrapper_path))
source(file.path(project_dir, "Scripts", "event-study-inference-helpers.R"))

template_path <- "/Users/nikos/Desktop/peiraeus/piraeus-port-traffic-containers.html"
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
output_path <- file.path(project_dir, variant_config$output_html)
result_path <- file.path(project_dir, ".traffic-work", paste0(variant_config$id, "-results.json"))
if (!file.exists(template_path)) stop("Container HTML template not found: ", template_path)
reuse_honest <- identical(Sys.getenv("CONTAINER_REUSE_HONEST", unset = "false"), "true") && file.exists(result_path)
restore_json_nulls <- function(value) {
  if (is.list(value) && length(value) == 0L) return(NA_real_)
  if (is.list(value)) return(lapply(value, restore_json_nulls))
  value
}
cached_bundle <- if (reuse_honest) restore_json_nulls(read_json(result_path, simplifyVector = FALSE)) else NULL

TREATMENT_INDEX <- 11L
Q3_INDEX <- 10L
Q2_EVENT_TIME <- -2L
PRE_END_YEAR <- 2019L

donor_pools <- list(
  "Core container" = c("Volos", "Thessaloniki"),
  "Core + Heraklio" = c("Volos", "Thessaloniki", "Heraklio"),
  "All container-active" = c("Volos", "Thessaloniki", "Heraklio", "Lavrio")
)

raw <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
outcome_config <- analysis_outcome_config()
outcome_config$key <- "containers"
raw <- set_analysis_outcome(raw, outcome_config)
raw$quarter_index <- (raw$Year - 2014L) * 4L + raw$Quarter_Number - 1L
raw$event_time <- raw$quarter_index - TREATMENT_INDEX
raw$treated <- as.integer(raw$Port == "Piraeus")
raw$post <- as.integer(raw$event_time >= 0L)
raw$log_total <- ifelse(raw$total_cargo > 0, log(raw$total_cargo), NA_real_)

# Match the original dashboard's full-sample balanced, positive container panel.
eligible <- ave(
  raw$total_cargo > 0,
  raw$Port,
  FUN = function(value) length(value) == 48L && all(value)
)
balanced <- raw[eligible & !grepl("^GR -", raw$Port),]

make_panel <- function(end_year, ports = NULL, force_drop_q3 = variant_config$drop_2016_q3) {
  panel <- balanced[balanced$Year <= end_year,]
  if (!is.null(ports)) panel <- panel[panel$Port %in% ports,]
  if (force_drop_q3) panel <- panel[panel$quarter_index != Q3_INDEX,]
  panel$Quarter <- factor(panel$Quarter, levels = unique(panel$Quarter[order(panel$quarter_index)]))
  panel
}

event_term_local <- function(k) paste0("event_time::", k, ":treated")

joint_test <- function(estimates, covariance, indexes, cluster_count) {
  b <- estimates[indexes]
  v <- covariance[indexes, indexes, drop = FALSE]
  decomposition <- eigen((v + t(v)) / 2, symmetric = TRUE)
  tolerance <- max(c(decomposition$values, 0)) * 1e-8
  keep <- decomposition$values > tolerance
  if (!any(keep)) return(list(f = NA_real_, p = NA_real_, rank = 0L))
  projected <- crossprod(decomposition$vectors[, keep, drop = FALSE], b)
  statistic <- sum(projected^2 / decomposition$values[keep]) / sum(keep)
  list(f = statistic, p = pf(statistic, sum(keep), cluster_count - 1L, lower.tail = FALSE), rank = sum(keep))
}

linear_test <- function(estimates, covariance, indexes, cluster_count) {
  weights <- rep(1 / length(indexes), length(indexes))
  estimate <- sum(weights * estimates[indexes])
  standard_error <- sqrt(max(as.numeric(t(weights) %*% covariance[indexes, indexes, drop = FALSE] %*% weights), 0))
  statistic <- estimate / standard_error
  list(estimate = estimate, se = standard_error, p = 2 * pt(-abs(statistic), cluster_count - 1L))
}

fit_event <- function(panel, full = FALSE, normalization = variant_config$normalization) {
  model <- feols(
    log_total ~ i(event_time, treated, ref = Q2_EVENT_TIME) | Port + Quarter,
    data = panel,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  all_times <- if (full) -11L:36L else -11L:12L
  available_times <- all_times[all_times != -1L | !variant_config$drop_2016_q3]
  coefficient_names <- names(coef(model))
  selector <- matrix(0, length(available_times), length(coefficient_names),
                     dimnames = list(as.character(available_times), coefficient_names))
  for (i in seq_along(available_times)) {
    term <- event_term_local(available_times[i])
    if (term %in% coefficient_names) selector[i, term] <- 1
  }
  contrast <- selector
  omitted <- available_times == Q2_EVENT_TIME
  if (normalization == "average_pre") {
    pre_rows <- which(available_times < 0L)
    contrast <- sweep(selector, 2, colMeans(selector[pre_rows,, drop = FALSE]), "-")
    omitted[] <- FALSE
  }
  beta <- coef(model)[coefficient_names]
  covariance_raw <- vcov(model)[coefficient_names, coefficient_names, drop = FALSE]
  estimates <- as.numeric(contrast %*% beta)
  covariance <- contrast %*% covariance_raw %*% t(contrast)
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  cluster_count <- length(unique(panel$Port))
  post_rows <- which(available_times >= 0L)
  pre_rows <- which(available_times < 0L & !omitted)
  average_post <- linear_test(estimates, covariance, post_rows, cluster_count)
  post_joint <- joint_test(estimates, covariance, post_rows, cluster_count)
  pre_joint <- joint_test(estimates, covariance, pre_rows, cluster_count)
  rows <- lapply(seq_along(available_times), function(i) list(
    quarter = sprintf("%d Q%d", 2014L + (available_times[i] + TREATMENT_INDEX) %/% 4L,
                      (available_times[i] + TREATMENT_INDEX) %% 4L + 1L),
    eventTime = available_times[i],
    estimate = estimates[i],
    se = if (omitted[i]) NA_real_ else standard_errors[i],
    low = if (omitted[i]) NA_real_ else estimates[i] - qnorm(.975) * standard_errors[i],
    high = if (omitted[i]) NA_real_ else estimates[i] + qnorm(.975) * standard_errors[i],
    p = if (omitted[i]) NA_real_ else 2 * pnorm(-abs(estimates[i] / standard_errors[i])),
    omitted = omitted[i]
  ))
  list(model = model, rows = rows, times = available_times, estimates = estimates,
       covariance = covariance, average_post = average_post, post_joint = post_joint,
       pre_joint = pre_joint, cluster_count = cluster_count)
}

fit_did <- function(panel, seasonal = FALSE) {
  formula <- if (seasonal) {
    log_total ~ treated:post | Port^Quarter_Number + Quarter
  } else {
    log_total ~ treated:post | Port + Quarter
  }
  model <- feols(formula, data = panel, cluster = ~Port,
                 ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE))
  table <- coeftable(model)
  list(estimate = unname(table["treated:post", "Estimate"]),
       se = unname(table["treated:post", "Std. Error"]),
       p = unname(table["treated:post", "Pr(>|t|)"]))
}

fit_seasonal_event <- function(panel, full = FALSE) {
  all_times <- if (full) -11L:36L else -11L:12L
  available_times <- all_times[all_times != -1L | !variant_config$drop_2016_q3]
  base_times <- available_times[1:4]
  model <- feols(
    log_total ~ i(event_time, treated, ref = base_times) | Port^Quarter_Number + Quarter,
    data = panel,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  quarter_of_year <- ((available_times + TREATMENT_INDEX) %% 4L) + 1L
  contrasts <- seasonal_average_event_study(
    model,
    event_times = available_times,
    quarter_of_year = quarter_of_year,
    pre_event_times = available_times[available_times < 0L],
    cluster_count = length(unique(panel$Port))
  )
  rows <- lapply(seq_along(available_times), function(i) list(
    quarter = sprintf("%d Q%d", 2014L + (available_times[i] + TREATMENT_INDEX) %/% 4L,
                      (available_times[i] + TREATMENT_INDEX) %% 4L + 1L),
    eventTime = available_times[i], estimate = contrasts$estimates[i],
    se = contrasts$standard_errors[i], low = contrasts$ci_low[i], high = contrasts$ci_high[i],
    p = contrasts$p_values[i], omitted = FALSE
  ))
  list(rows = rows, contrasts = contrasts)
}

tests_object <- function(panel, event_fit, did_fit, seasonal = FALSE) {
  out <- list(
    didPercent = 100 * (exp(did_fit$estimate) - 1),
    averagePostPercent = 100 * (exp(event_fit$average_post$estimate) - 1),
    averagePostP = event_fit$average_post$p,
    jointPostP = event_fit$post_joint$p,
    pretrendP = event_fit$pre_joint$p
  )
  if (seasonal) {
    out$didClusteredP <- did_fit$p
  } else {
    out$wildBootstrapP <- as.numeric(wild_cluster_did(panel)["p_value"])
  }
  out
}

study_object <- function(panel, event_fit, did_fit, full = FALSE) {
  quarters <- length(unique(panel$quarter_index))
  base_label <- if (variant_config$normalization == "q2") "2016 Q2" else "Average of fully untreated pre-treatment quarters"
  list(
    meta = list(
      sample = if (full) "2014 Q1–2025 Q4" else "2014 Q1–2019 Q4",
      treatmentQuarter = "2016 Q4", baseQuarter = base_label,
      outcome = "ln(container cargo tonnes)", ports = event_fit$cluster_count,
      controls = event_fit$cluster_count - 1L, quarters = quarters,
      observations = nrow(panel), standardErrors = "Port-clustered CR1",
      sampleRule = if (variant_config$drop_2016_q3)
        paste0("Balanced positive-outcome ports; 2016 Q3 excluded (", quarters, " retained quarters)")
      else paste0("Balanced ports with positive outcome in all ", quarters, " quarters"),
      source = "data/port-cargo-master-2014-2025.xlsx"
    ),
    estimates = event_fit$rows,
    tests = tests_object(panel, event_fit, did_fit)
  )
}

seasonal_study_object <- function(panel, fit, did_fit, baseline_did) {
  cfit <- fit$contrasts
  proxy <- list(average_post = list(estimate = as.numeric(cfit$average_post["estimate"]),
                                    p = as.numeric(cfit$average_post["p_value"])),
                post_joint = list(p = as.numeric(cfit$joint_post["p_value"])),
                pre_joint = list(p = as.numeric(cfit$joint_pre["p_value"])))
  seasonal_did_pct <- 100 * (exp(did_fit$estimate) - 1)
  list(
    meta = list(sample = "2014 Q1–2019 Q4", treatmentQuarter = "2016 Q4",
                reference = if (variant_config$drop_2016_q3)
                  "same-quarter fully untreated pre-treatment average; 2016 Q3 excluded"
                else "same-quarter pre-treatment average",
                outcome = "ln(container cargo tonnes)", ports = length(unique(panel$Port)),
                controls = length(unique(panel$Port)) - 1L,
                quarters = length(unique(panel$quarter_index)), observations = nrow(panel),
                fixedEffects = "port × quarter-of-year + calendar quarter",
                standardErrors = "Port-clustered CR1 linear-contrast intervals"),
    estimates = fit$rows,
    tests = c(tests_object(panel, proxy, did_fit, seasonal = TRUE)),
    comparison = list(baselineDidPercent = baseline_did,
                      seasonalDidPercent = seasonal_did_pct,
                      differencePoints = seasonal_did_pct - baseline_did)
  )
}

donor_object <- function(end_year, seasonal = FALSE) {
  pools <- lapply(names(donor_pools), function(pool_name) {
    donors <- donor_pools[[pool_name]]
    panel <- make_panel(end_year, c("Piraeus", donors))
    if (seasonal) {
      event_fit <- fit_seasonal_event(panel, full = FALSE)
      did_fit <- fit_did(panel, seasonal = TRUE)
      cfit <- event_fit$contrasts
      proxy <- list(average_post = list(estimate = as.numeric(cfit$average_post["estimate"]),
                                        p = as.numeric(cfit$average_post["p_value"])),
                    post_joint = list(p = as.numeric(cfit$joint_post["p_value"])),
                    pre_joint = list(p = as.numeric(cfit$joint_pre["p_value"])))
      test <- tests_object(panel, proxy, did_fit, seasonal = TRUE)
      baseline <- fit_did(panel, seasonal = FALSE)
      omitted <- rep(FALSE, length(event_fit$rows))
    } else {
      event_fit <- fit_event(panel, full = end_year == 2025L)
      did_fit <- fit_did(panel)
      test <- tests_object(panel, event_fit, did_fit)
      baseline <- NULL
      omitted <- vapply(event_fit$rows, function(x) x$omitted, logical(1))
    }
    list(
      name = pool_name, donors = donors,
      did = 100 * (exp(did_fit$estimate) - 1),
      b = vapply(event_fit$rows, function(x) x$estimate, numeric(1)),
      lo = vapply(event_fit$rows, function(x) x$low, numeric(1)),
      hi = vapply(event_fit$rows, function(x) x$high, numeric(1)),
      tests = test,
      baselineDid = if (seasonal) 100 * (exp(baseline$estimate) - 1) else NULL,
      differencePoints = if (seasonal) 100 * (exp(did_fit$estimate) - exp(baseline$estimate)) else NULL,
      omitted = omitted
    )
  })
  list(pools = pools)
}

# HonestDiD requires the omitted period to be the final clean pre-period.
# Therefore both variants use a dedicated sample with 2016 Q3 removed and Q2
# as the computational reference. This is stated explicitly in the dashboard.
run_honest <- function(panel, pool_name = NULL, donors = NULL) {
  honest_panel <- panel[panel$quarter_index != Q3_INDEX,]
  honest_panel$Quarter <- droplevels(honest_panel$Quarter)
  model <- feols(
    log_total ~ i(event_time, treated, ref = Q2_EVENT_TIME) | Port + Quarter,
    data = honest_panel, cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  pre_times <- -11L:-3L
  post_times <- 0L:12L
  times <- c(pre_times, post_times)
  terms <- vapply(times, event_term_local, character(1))
  beta <- coef(model)[terms]
  sigma <- vcov(model)[terms, terms, drop = FALSE]
  l_vec <- rep(1 / length(post_times), length(post_times))
  rm <- suppressWarnings(HonestDiD::createSensitivityResults_relativeMagnitudes(
    betahat = beta, sigma = sigma, numPrePeriods = length(pre_times),
    numPostPeriods = length(post_times), method = "C-LF",
    Mbarvec = c(0, .025, .05, .1, .25, .5), l_vec = l_vec,
    gridPoints = 1201, grid.lb = -10, grid.ub = 10, seed = 20260822
  ))
  rm <- as.data.frame(rm)
  sd <- do.call(rbind, lapply(c(0, .025, .05, .1), function(m_value) {
    cs <- suppressWarnings(HonestDiD::computeConditionalCS_DeltaSD(
      betahat = beta, sigma = sigma, numPrePeriods = length(pre_times),
      numPostPeriods = length(post_times), l_vec = l_vec, M = m_value,
      alpha = .05, hybrid_flag = "LF", gridPoints = 4001,
      grid.lb = -20, grid.ub = 20, seed = 20260822
    ))
    accepted <- cs$grid[cs$accept == 1]
    data.frame(M = m_value, lb = min(accepted), ub = max(accepted))
  }))
  post_terms <- vapply(post_times, event_term_local, character(1))
  post_beta <- coef(model)[post_terms]
  post_v <- vcov(model)[post_terms, post_terms, drop = FALSE]
  w <- rep(1 / length(post_terms), length(post_terms))
  estimate <- sum(w * post_beta)
  se <- sqrt(as.numeric(t(w) %*% post_v %*% w))
  rm_rows <- lapply(seq_len(nrow(rm)), function(i) list(
    value = rm$Mbar[i], low = rm$lb[i], high = rm$ub[i],
    includesZero = rm$lb[i] <= 0 && rm$ub[i] >= 0
  ))
  sd_rows <- lapply(seq_len(nrow(sd)), function(i) list(
    value = sd$M[i], low = sd$lb[i], high = sd$ub[i],
    includesZero = sd$lb[i] <= 0 && sd$ub[i] >= 0
  ))
  summary <- list(
    estimate = estimate, percentEffect = 100 * (exp(estimate) - 1),
    conventionalLow = estimate - qnorm(.975) * se,
    conventionalHigh = estimate + qnorm(.975) * se,
    firstRelativeMagnitudeIncludingZero = {
      z <- which(vapply(rm_rows, function(x) x$includesZero, logical(1))); if (length(z)) rm_rows[[z[1]]]$value else NA_real_
    },
    lastRelativeMagnitudeExcludingZero = NA_real_,
    firstSmoothnessIncludingZero = {
      z <- which(vapply(sd_rows, function(x) x$includesZero, logical(1))); if (length(z)) sd_rows[[z[1]]]$value else NA_real_
    },
    lastSmoothnessExcludingZero = NA_real_
  )
  if (is.null(pool_name)) {
    list(
      meta = list(method = "Rambachan–Roth HonestDiD", packageVersion = as.character(packageVersion("HonestDiD")),
                  sample = "2014 Q1–2019 Q4; 2016 Q3 excluded", estimand = "Mean of 13 post-treatment event coefficients relative to 2016 Q2",
                  observedPreQuarters = 10L, estimatedPreCoefficients = 9L, postCoefficients = 13L,
                  inference = "C-LF robust 95% confidence sets; port-clustered covariance"),
      summary = summary, relativeMagnitude = rm_rows, smoothness = sd_rows
    )
  } else {
    c(list(name = pool_name, donors = donors), summary,
      list(relativeMagnitude = rm_rows, smoothness = sd_rows))
  }
}

fit_staggered <- function() {
  treatments <- data.frame(
    port = c("Piraeus", "Thessaloniki", "Heraklio"),
    treatment_quarter = c("2016 Q4", "2018 Q2", "2024 Q4"),
    treatment_index = c(11L, 17L, 43L),
    event = c("Transfer of 51% of PPA shares", "Transfer of 67% of ThPA shares", "Transfer of 67% of HPA shares"),
    source = c("https://hradf.com/en/piraeus-port-authority-s-a-ppa/",
               "https://hradf.com/wp-content/uploads/2021/11/hradf-thpa-share-transfer.pdf",
               "https://growthfund.gr/en/hradf-acquisition-of-a-majority-stake-in-the-share-capital-of-heraklion-port-authority-hpa-s-a-by-the-consortium-grimaldi-euromed-s-p-a-minoan-lines-s-a-for-80-million-euros/"),
    stringsAsFactors = FALSE
  )
  panel <- make_panel(2025L)
  cohort_lookup <- setNames(treatments$treatment_index, treatments$port)
  panel$cohort <- ifelse(panel$Port %in% treatments$port, unname(cohort_lookup[panel$Port]), 1000L)

  fit_one <- function(data) feols(
    log_total ~ sunab(cohort, quarter_index, ref.p = Q2_EVENT_TIME) | Port + Quarter,
    data = data, cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )

  summarize_model <- function(model, included_treatments, normalization, model_cluster_count) {
    beta <- model$coefficients
    covariance <- vcov(model)
    coefficient_names <- names(beta)
    cells <- list()
    for (g in included_treatments$treatment_index) {
      calendar_times <- sort(unique(panel$quarter_index))
      rel_times <- calendar_times - g
      rel_times <- rel_times[rel_times >= -12L & rel_times <= 36L]
      selector <- matrix(0, length(rel_times), length(beta), dimnames = list(as.character(rel_times), coefficient_names))
      for (i in seq_along(rel_times)) {
        term <- paste0("quarter_index::", rel_times[i], ":cohort::", g)
        if (term %in% coefficient_names) selector[i, term] <- 1
      }
      contrast <- selector
      if (normalization == "average_pre") {
        pre <- which(rel_times < 0L)
        contrast <- sweep(selector, 2, colMeans(selector[pre,, drop = FALSE]), "-")
      }
      cells[[as.character(g)]] <- list(times = rel_times, contrast = contrast)
    }
    cell_rows <- list()
    for (g in names(cells)) {
      item <- cells[[g]]
      for (i in seq_along(item$times)) {
        k <- item$times[i]
        if (normalization == "q2" && k == Q2_EVENT_TIME) next
        cell_rows[[length(cell_rows) + 1L]] <- list(g = as.integer(g), k = k, c = item$contrast[i,])
      }
    }
    aggregate_contrast <- function(rows) Reduce("+", lapply(rows, `[[`, "c")) / length(rows)
    estimate_contrast <- function(cvec) {
      estimate <- sum(cvec * beta)
      se <- sqrt(max(as.numeric(t(cvec) %*% covariance %*% cvec), 0))
      list(estimate = estimate, se = se, p = 2 * pt(-abs(estimate / se), model_cluster_count - 1L))
    }
    period_rows <- lapply(-12L:12L, function(k) {
      if (normalization == "q2" && k == Q2_EVENT_TIME) {
        return(list(event_time = k, estimate = 0, std_error = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                    p_value = NA_real_, contributing_cohorts = nrow(included_treatments), omitted = TRUE))
      }
      matching <- Filter(function(x) x$k == k, cell_rows)
      if (!length(matching)) return(NULL)
      value <- estimate_contrast(aggregate_contrast(matching))
      list(event_time = k, estimate = value$estimate, std_error = value$se,
           ci_low = value$estimate - qt(.975, model_cluster_count - 1L) * value$se,
           ci_high = value$estimate + qt(.975, model_cluster_count - 1L) * value$se,
           p_value = value$p, contributing_cohorts = length(matching), omitted = FALSE)
    })
    period_rows <- Filter(Negate(is.null), period_rows)
    post_cells <- Filter(function(x) x$k >= 0L, cell_rows)
    att <- estimate_contrast(aggregate_contrast(post_cells))
    cohort_rows <- lapply(seq_len(nrow(included_treatments)), function(i) {
      g <- included_treatments$treatment_index[i]
      matching <- Filter(function(x) x$g == g && x$k >= 0L, cell_rows)
      value <- estimate_contrast(aggregate_contrast(matching))
      c(as.list(included_treatments[i,]), list(estimate = value$estimate, std_error = value$se,
        ci_low = value$estimate - qt(.975, model_cluster_count - 1L) * value$se,
        ci_high = value$estimate + qt(.975, model_cluster_count - 1L) * value$se,
        p_value = value$p, post_quarters = length(matching), percent_effect = 100 * (exp(value$estimate) - 1)))
    })
    pre_period <- Filter(function(x) x$event_time < 0L && !x$omitted, period_rows)
    pre_p <- NA_real_
    if (length(pre_period)) {
      # A compact Wald diagnostic on the displayed aggregate lead estimates.
      pre_cs <- lapply(pre_period, function(row) {
        matching <- Filter(function(x) x$k == row$event_time, cell_rows)
        aggregate_contrast(matching)
      })
      C <- do.call(rbind, pre_cs)
      b <- as.numeric(C %*% beta)
      V <- C %*% covariance %*% t(C)
      decomp <- eigen((V + t(V)) / 2, symmetric = TRUE)
      keep <- decomp$values > max(c(decomp$values, 0)) * 1e-8
      if (any(keep)) {
        projected <- crossprod(decomp$vectors[, keep, drop = FALSE], b)
        f <- sum(projected^2 / decomp$values[keep]) / sum(keep)
        pre_p <- pf(f, sum(keep), model_cluster_count - 1L, lower.tail = FALSE)
      }
    }
    list(att = att, periods = period_rows, cohorts = cohort_rows, pre_p = pre_p)
  }

  model <- fit_one(panel)
  summary <- summarize_model(model, treatments, variant_config$normalization, length(unique(panel$Port)))
  leave_one_out <- lapply(treatments$port, function(port) {
    subset_treatments <- treatments[treatments$port != port,]
    subset_panel <- panel[panel$Port != port,]
    subset_model <- fit_one(subset_panel)
    value <- summarize_model(subset_model, subset_treatments, variant_config$normalization,
                             length(unique(subset_panel$Port)))$att
    list(omitted_port = port, estimate = value$estimate, std_error = value$se,
         p_value = value$p, percent_effect = 100 * (exp(value$estimate) - 1))
  })
  cluster_count <- length(unique(panel$Port))
  list(
    meta = list(estimator = "Sun-Abraham interaction-weighted; cohort paths re-normalized",
                sample = if (variant_config$drop_2016_q3) "2014 Q1-2025 Q4 container cargo; 2016 Q3 excluded" else "2014 Q1-2025 Q4 container cargo",
                outcome = "ln(container cargo tonnes)", ports = cluster_count,
                treatedPorts = nrow(treatments), neverTreatedPorts = cluster_count - nrow(treatments),
                quarters = length(unique(panel$quarter_index)), observations = nrow(panel),
                standardErrors = "Port-clustered CR1",
                sampleRule = "Balanced positive-container-cargo ports in retained quarters",
                eventWindow = if (variant_config$normalization == "q2") "k = -12 to +12; k = -2 omitted" else "k = -12 to +12; cohort-specific fully untreated pre-period average = 0",
                comparisonGroup = "Never-treated ports", source = "data/port-cargo-master-2014-2025.xlsx"),
    summary = list(estimate = summary$att$estimate, percentEffect = 100 * (exp(summary$att$estimate) - 1),
                   standardError = summary$att$se,
                   ciLow = summary$att$estimate - qt(.975, cluster_count - 1L) * summary$att$se,
                   ciHigh = summary$att$estimate + qt(.975, cluster_count - 1L) * summary$att$se,
                   pValue = summary$att$p, pretrendF = NA_real_, pretrendDf1 = NA_real_,
                   pretrendDf2 = cluster_count - 1L, pretrendP = summary$pre_p),
    treatments = lapply(seq_len(nrow(treatments)), function(i) as.list(treatments[i,])),
    cohorts = summary$cohorts, events = summary$periods, leaveOneOut = leave_one_out
  )
}

pre_panel <- make_panel(PRE_END_YEAR)
full_panel <- make_panel(2025L)
pre_event <- fit_event(pre_panel)
pre_did <- fit_did(pre_panel)
full_event <- fit_event(full_panel, full = TRUE)
full_did <- fit_did(full_panel)
seasonal_event <- fit_seasonal_event(pre_panel)
seasonal_did <- fit_did(pre_panel, seasonal = TRUE)

eventStudy <- study_object(pre_panel, pre_event, pre_did)
fullSampleEventStudy <- study_object(full_panel, full_event, full_did, full = TRUE)
seasonalEventStudy <- seasonal_study_object(
  pre_panel, seasonal_event, seasonal_did, eventStudy$tests$didPercent
)
donorPoolStudy <- donor_object(PRE_END_YEAR)
fullSampleDonorPoolStudy <- donor_object(2025L)
seasonalDonorPoolStudy <- donor_object(PRE_END_YEAR, seasonal = TRUE)
if (reuse_honest) {
  honestDidStudy <- cached_bundle$honestDidStudy
  donorPoolHonestDidStudy <- cached_bundle$donorPoolHonestDidStudy
} else {
  honestDidStudy <- run_honest(make_panel(PRE_END_YEAR, force_drop_q3 = TRUE))
  donor_honest_pools <- lapply(names(donor_pools), function(pool_name) {
    donors <- donor_pools[[pool_name]]
    run_honest(make_panel(PRE_END_YEAR, c("Piraeus", donors), force_drop_q3 = TRUE), pool_name, donors)
  })
  donorPoolHonestDidStudy <- list(
    method = "Rambachan–Roth HonestDiD",
    estimand = "Mean of 13 post-treatment coefficients relative to 2016 Q2",
    sample = "2014 Q1–2019 Q4; 2016 Q3 excluded",
    pools = donor_honest_pools
  )
}
staggeredDidStudy <- fit_staggered()

bundle <- list(
  eventStudy = eventStudy,
  honestDidStudy = honestDidStudy,
  donorPoolStudy = donorPoolStudy,
  donorPoolHonestDidStudy = donorPoolHonestDidStudy,
  staggeredDidStudy = staggeredDidStudy,
  fullSampleEventStudy = fullSampleEventStudy,
  fullSampleDonorPoolStudy = fullSampleDonorPoolStudy,
  seasonalEventStudy = seasonalEventStudy,
  seasonalDonorPoolStudy = seasonalDonorPoolStudy
)

dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
write_json(bundle, result_path, auto_unbox = TRUE, pretty = TRUE, digits = 15, na = "null")

# The HTML assembly and static-text refresh are kept in a small deterministic
# builder so the statistical work remains entirely visible in these R files.
builder <- file.path(project_dir, ".traffic-work", "build-container-variant.py")
status <- system2("python3", c(builder, shQuote(template_path), shQuote(result_path),
                               shQuote(output_path), shQuote(variant_config$title),
                               shQuote(variant_config$normalization),
                               if (variant_config$drop_2016_q3) "true" else "false"))
if (status != 0) stop("HTML builder failed with status ", status)

cat("Saved analysis bundle:", result_path, "\n")
cat("Saved dashboard:", output_path, "\n")
cat("Primary DiD:", sprintf("%+.1f%%", eventStudy$tests$didPercent), "\n")
cat("Average post event-study effect:", sprintf("%+.1f%%", eventStudy$tests$averagePostPercent), "\n")
