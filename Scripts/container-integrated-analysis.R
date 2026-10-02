# Integrated container-cargo report.
#
# This single script contains every estimation used in the accompanying HTML:
#   1. Seasonally adjusted Piraeus DiD and event study.
#   2. Within-Piraeus container versus non-container comparison.
#   3. Short- and full-sample seasonally adjusted all-port event studies.
#   4. Donor-pool sensitivity for both windows.
#   5. Seasonally adjusted stacked staggered DiD for Piraeus, Thessaloniki,
#      and Heraklio, using Volos and Lavrio as never-treated controls.
#
# Timing convention: the 10 August 2016 transfer occurred in 2016 Q3. That
# transition quarter is excluded. 2016 Q2 is the last clean pre-treatment
# quarter and 2016 Q4 is the first post-treatment quarter.

suppressPackageStartupMessages({
  library(readxl)
  library(fixest)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
template_path <- file.path(project_dir, "Scripts", "container-integrated-template.html")
output_path <- file.path(project_dir, "piraeus-container-integrated-analysis.html")
json_path <- file.path(project_dir, ".traffic-work", "container-integrated-results.json")
csv_path <- file.path(project_dir, ".traffic-work", "container-integrated-estimates.csv")

PIRAEUS_TREATMENT <- 11L       # 2016 Q4, first full post-treatment quarter
TRANSITION_QUARTER <- 10L      # 2016 Q3, excluded
LAST_CLEAN_PRE <- 9L           # 2016 Q2
HAC_LAG <- 4L

total_columns <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)
container_columns <- c(
  "External_Unloaded_Containers", "External_Loaded_Containers",
  "Internal_Unloaded_Containers", "Internal_Loaded_Containers"
)

raw <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
required <- unique(c(total_columns, container_columns))
if (length(setdiff(required, names(raw)))) stop("Required cargo columns are missing")
raw[required] <- lapply(raw[required], function(x) { x[is.na(x)] <- 0; as.numeric(x) })
raw$total_cargo <- rowSums(raw[total_columns])
raw$container <- rowSums(raw[container_columns])
raw$noncontainer <- raw$total_cargo - raw$container
raw$quarter_index <- (raw$Year - 2014L) * 4L + raw$Quarter_Number - 1L
raw$event_time <- raw$quarter_index - PIRAEUS_TREATMENT
raw$treated <- as.integer(raw$Port == "Piraeus")
raw$post <- as.integer(raw$quarter_index >= PIRAEUS_TREATMENT)
raw$log_total <- ifelse(raw$container > 0, log(raw$container), NA_real_)

eligible <- ave(
  is.finite(raw$log_total), raw$Port,
  FUN = function(x) length(x) == 48L && all(x)
)
balanced <- raw[eligible & !grepl("^GR -", raw$Port),]
balanced <- balanced[balanced$quarter_index != TRANSITION_QUARTER,]
balanced <- balanced[order(balanced$Port, balanced$quarter_index),]
balanced$Quarter <- factor(
  balanced$Quarter,
  levels = unique(balanced$Quarter[order(balanced$quarter_index)])
)

all_ports <- sort(unique(balanced$Port))
expected_ports <- c("Heraklio", "Lavrio", "Piraeus", "Thessaloniki", "Volos")
if (!setequal(all_ports, expected_ports)) {
  stop("Unexpected balanced positive-container port set: ", paste(all_ports, collapse = ", "))
}

donor_pools <- list(
  "Core container" = c("Volos", "Thessaloniki"),
  "Core + Heraklio" = c("Volos", "Thessaloniki", "Heraklio"),
  "All container-active" = c("Volos", "Thessaloniki", "Heraklio", "Lavrio")
)

event_term <- function(k) paste0("event_time::", k, ":treated")

contrast_result <- function(cvec, beta, covariance, df) {
  estimate <- sum(cvec * beta)
  variance <- as.numeric(t(cvec) %*% covariance %*% cvec)
  se <- sqrt(max(variance, 0))
  statistic <- estimate / se
  critical <- qt(.975, df)
  list(
    estimate = estimate,
    percent = 100 * (exp(estimate) - 1),
    standardError = se,
    ciLow = estimate - critical * se,
    ciHigh = estimate + critical * se,
    percentCiLow = 100 * (exp(estimate - critical * se) - 1),
    percentCiHigh = 100 * (exp(estimate + critical * se) - 1),
    pValue = 2 * pt(-abs(statistic), df)
  )
}

joint_test <- function(C, beta, covariance, df) {
  b <- as.numeric(C %*% beta)
  V <- C %*% covariance %*% t(C)
  eig <- eigen((V + t(V)) / 2, symmetric = TRUE)
  tolerance <- max(c(eig$values, 0)) * 1e-8
  keep <- eig$values > tolerance
  if (!any(keep)) return(list(fStatistic = NA_real_, rank = 0L, pValue = NA_real_))
  projected <- crossprod(eig$vectors[, keep, drop = FALSE], b)
  f <- sum(projected^2 / eig$values[keep]) / sum(keep)
  list(fStatistic = f, rank = sum(keep), pValue = pf(f, sum(keep), df, lower.tail = FALSE))
}

make_panel <- function(end_year, donors = NULL) {
  ports <- if (is.null(donors)) all_ports else c("Piraeus", donors)
  panel <- balanced[balanced$Year <= end_year & balanced$Port %in% ports,]
  panel$Quarter <- droplevels(panel$Quarter)
  panel
}

# Seasonality-adjusted event-study coefficients are normalized against the
# average clean pre-treatment observation from the same quarter of the year.
# This is the identified normalization when port-specific quarter-of-year FE
# are combined with a fully flexible event path. A second contrast shifts the
# complete path to Piraeus's 2016 Q2 estimate.
fit_seasonal_study <- function(panel, label) {
  event_times <- sort(unique(panel$event_time))
  pre_times <- event_times[event_times < 0L]
  post_times <- event_times[event_times >= 0L]
  qnum_for_time <- ((event_times + PIRAEUS_TREATMENT) %% 4L) + 1L
  base_times <- vapply(1:4, function(q) as.integer(min(event_times[qnum_for_time == q])), integer(1))

  event_model <- feols(
    log_total ~ i(event_time, treated, ref = base_times) | Port^Quarter_Number + Quarter,
    data = panel,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  beta <- coef(event_model)
  covariance <- vcov(event_model)
  coefficient_names <- names(beta)
  selector <- matrix(0, nrow = length(event_times), ncol = length(beta),
                     dimnames = list(as.character(event_times), coefficient_names))
  for (i in seq_along(event_times)) {
    term <- event_term(event_times[i])
    if (term %in% coefficient_names) selector[i, term] <- 1
  }
  seasonal_contrast <- selector
  for (i in seq_along(event_times)) {
    same_season_pre <- which(event_times < 0L & qnum_for_time == qnum_for_time[i])
    seasonal_contrast[i,] <- selector[i,] - colMeans(selector[same_season_pre,, drop = FALSE])
  }
  q2_row <- match(LAST_CLEAN_PRE - PIRAEUS_TREATMENT, event_times)
  q2_contrast <- sweep(seasonal_contrast, 2, seasonal_contrast[q2_row,], "-")
  cluster_count <- length(unique(panel$Port))
  df <- cluster_count - 1L
  estimates <- as.numeric(seasonal_contrast %*% beta)
  V_event <- seasonal_contrast %*% covariance %*% t(seasonal_contrast)
  q2_estimates <- as.numeric(q2_contrast %*% beta)
  V_q2 <- q2_contrast %*% covariance %*% t(q2_contrast)
  critical <- qt(.975, df)

  post_rows <- which(event_times >= 0L)
  pre_rows <- which(event_times < 0L)
  avg_pre_c <- colMeans(seasonal_contrast[post_rows,, drop = FALSE])
  q2_avg_c <- colMeans(q2_contrast[post_rows,, drop = FALSE])

  did_model <- feols(
    log_total ~ treated:post | Port^Quarter_Number + Quarter,
    data = panel,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  did_table <- coeftable(did_model)
  did_est <- unname(did_table["treated:post", "Estimate"])
  did_se <- unname(did_table["treated:post", "Std. Error"])
  did_p <- unname(did_table["treated:post", "Pr(>|t|)"])
  did_critical <- qt(.975, df)
  did <- list(
    estimate = did_est, percent = 100 * (exp(did_est) - 1), standardError = did_se,
    percentCiLow = 100 * (exp(did_est - did_critical * did_se) - 1),
    percentCiHigh = 100 * (exp(did_est + did_critical * did_se) - 1), pValue = did_p
  )

  event_rows <- lapply(seq_along(event_times), function(i) list(
    quarter = sprintf("%d Q%d", 2014L + (event_times[i] + PIRAEUS_TREATMENT) %/% 4L,
                      (event_times[i] + PIRAEUS_TREATMENT) %% 4L + 1L),
    eventTime = event_times[i],
    estimate = estimates[i],
    standardError = sqrt(max(V_event[i, i], 0)),
    ciLow = estimates[i] - critical * sqrt(max(V_event[i, i], 0)),
    ciHigh = estimates[i] + critical * sqrt(max(V_event[i, i], 0)),
    q2Estimate = q2_estimates[i],
    q2StandardError = sqrt(max(V_q2[i, i], 0)),
    q2CiLow = q2_estimates[i] - critical * sqrt(max(V_q2[i, i], 0)),
    q2CiHigh = q2_estimates[i] + critical * sqrt(max(V_q2[i, i], 0))
  ))

  list(
    label = label,
    meta = list(ports = cluster_count, controls = cluster_count - 1L,
                quarters = length(unique(panel$quarter_index)), observations = nrow(panel),
                fixedEffects = "port × quarter-of-year + calendar quarter",
                standardErrors = "port-clustered CR1", excludedQuarter = "2016 Q3"),
    did = did,
    averagePre = contrast_result(avg_pre_c, beta, covariance, df),
    q2 = contrast_result(q2_avg_c, beta, covariance, df),
    pretrend = joint_test(seasonal_contrast[pre_rows,, drop = FALSE], beta, covariance, df),
    events = event_rows
  )
}

fit_sensitivity <- function(end_year) {
  pools <- donor_pools
  if (end_year == 2019L) {
    pools[["Untreated through 2019"]] <- c("Volos", "Heraklio", "Lavrio")
  } else {
    pools[["Never-treated only"]] <- c("Volos", "Lavrio")
  }
  lapply(names(pools), function(name) {
    study <- fit_seasonal_study(make_panel(end_year, pools[[name]]), name)
    list(name = name, donors = pools[[name]], meta = study$meta, did = study$did,
         averagePre = study$averagePre, q2 = study$q2, pretrend = study$pretrend,
         events = study$events)
  })
}

# ---------- Within-Piraeus segment comparison ----------
hac_fit <- function(y, X, lag = 4L) {
  X <- as.matrix(X)
  fit <- lm.fit(X, y)
  u <- as.numeric(fit$residuals)
  n <- nrow(X); k <- ncol(X); xu <- X * u
  meat <- crossprod(xu)
  for (ell in seq_len(min(lag, n - 1L))) {
    weight <- 1 - ell / (lag + 1)
    gamma <- crossprod(xu[(ell + 1L):n,, drop = FALSE], xu[1L:(n - ell),, drop = FALSE])
    meat <- meat + weight * (gamma + t(gamma))
  }
  bread <- solve(crossprod(X))
  list(coefficients = setNames(as.numeric(fit$coefficients), colnames(X)),
       vcov = (n / (n - k)) * bread %*% meat %*% bread, df = n - k)
}

hac_result <- function(fit, term) {
  j <- match(term, names(fit$coefficients)); estimate <- fit$coefficients[j]
  se <- sqrt(fit$vcov[j,j]); critical <- qt(.975, fit$df)
  list(estimate = estimate, percent = 100 * (exp(estimate) - 1), standardError = se,
       percentCiLow = 100 * (exp(estimate - critical * se) - 1),
       percentCiHigh = 100 * (exp(estimate + critical * se) - 1),
       pValue = 2 * pt(-abs(estimate / se), fit$df))
}

piraeus <- raw[raw$Port == "Piraeus" & raw$quarter_index != TRANSITION_QUARTER,]
piraeus <- piraeus[order(piraeus$quarter_index),]
if (any(piraeus$container <= 0) || any(piraeus$noncontainer <= 0)) stop("Piraeus segment outcomes must be positive")
piraeus$log_ratio <- log(piraeus$container) - log(piraeus$noncontainer)
piraeus$post <- as.integer(piraeus$quarter_index >= PIRAEUS_TREATMENT)

within_window <- function(end_year) {
  d <- piraeus[piraeus$Year <= end_year,]
  pre <- d[d$post == 0L,]; post <- d[d$post == 1L,]
  seasonal_X <- model.matrix(~ post + factor(Quarter_Number), data = d)
  seasonal <- hac_result(hac_fit(d$log_ratio, seasonal_X, HAC_LAG), "post")
  list(
    label = if (end_year == 2019L) "2014 Q1–2019 Q4" else "2014 Q1–2025 Q4",
    quarters = nrow(d), preQuarters = nrow(pre), postQuarters = nrow(post),
    seasonalDid = seasonal
  )
}

within_pre <- piraeus[piraeus$post == 0L,]
within_pre_time <- within_pre$quarter_index - min(within_pre$quarter_index)
within_pretrend_X <- model.matrix(~ within_pre_time + factor(Quarter_Number), data = within_pre)
within_pretrend <- hac_result(
  hac_fit(within_pre$log_ratio, within_pretrend_X, HAC_LAG), "within_pre_time"
)
pre_container_by_quarter <- tapply(log(within_pre$container), within_pre$Quarter_Number, mean)
pre_noncontainer_by_quarter <- tapply(log(within_pre$noncontainer), within_pre$Quarter_Number, mean)
pre_ratio_by_quarter <- tapply(within_pre$log_ratio, within_pre$Quarter_Number, mean)
seasonal_ratio_gap <- piraeus$log_ratio - pre_ratio_by_quarter[as.character(piraeus$Quarter_Number)]
within_series <- lapply(seq_len(nrow(piraeus)), function(i) list(
  quarter = piraeus$Quarter[i], year = piraeus$Year[i], quarterNumber = piraeus$Quarter_Number[i],
  eventTime = piraeus$quarter_index[i] - PIRAEUS_TREATMENT,
  containerIndex = 100 * exp(log(piraeus$container[i]) - pre_container_by_quarter[as.character(piraeus$Quarter_Number[i])]),
  noncontainerIndex = 100 * exp(log(piraeus$noncontainer[i]) - pre_noncontainer_by_quarter[as.character(piraeus$Quarter_Number[i])]),
  ratioAveragePre = 100 * (exp(seasonal_ratio_gap[i]) - 1)
))
within <- list(
  windows = list(within_window(2019L), within_window(2025L)),
  pretrend = within_pretrend,
  series = within_series
)

# ---------- Seasonally adjusted staggered stacked DiD ----------
treatments <- data.frame(
  port = c("Piraeus", "Thessaloniki", "Heraklio"),
  treatmentQuarter = c("2016 Q4", "2018 Q2", "2024 Q4"),
  treatmentIndex = c(11L, 17L, 43L),
  stringsAsFactors = FALSE
)
never_treated <- c("Volos", "Lavrio")

stack_list <- lapply(seq_len(nrow(treatments)), function(i) {
  tr <- treatments[i,]
  d <- balanced[balanced$Port %in% c(tr$port, never_treated),]
  d$stack <- tr$port
  d$stackTreated <- as.integer(d$Port == tr$port)
  d$stackPost <- as.integer(d$quarter_index >= tr$treatmentIndex)
  d$stackPort <- interaction(d$stack, d$Port, drop = TRUE)
  d$stackQuarter <- interaction(d$stack, d$Quarter, drop = TRUE)
  d
})
stacked <- do.call(rbind, stack_list)
stacked_model <- feols(
  log_total ~ stackTreated:stackPost | stackPort^Quarter_Number + stackQuarter,
  data = stacked,
  cluster = ~Port,
  ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
)
stacked_table <- coeftable(stacked_model)
stacked_est <- unname(stacked_table["stackTreated:stackPost", "Estimate"])
stacked_se <- unname(stacked_table["stackTreated:stackPost", "Std. Error"])
stacked_df <- length(unique(stacked$Port)) - 1L
stacked_critical <- qt(.975, stacked_df)
stacked_summary <- list(
  estimate = stacked_est, percent = 100 * (exp(stacked_est) - 1), standardError = stacked_se,
  percentCiLow = 100 * (exp(stacked_est - stacked_critical * stacked_se) - 1),
  percentCiHigh = 100 * (exp(stacked_est + stacked_critical * stacked_se) - 1),
  pValue = unname(stacked_table["stackTreated:stackPost", "Pr(>|t|)"])
)

cohort_results <- lapply(seq_len(nrow(treatments)), function(i) {
  tr <- treatments[i,]
  d <- balanced[balanced$Port %in% c(tr$port, never_treated),]
  d$cohortTreated <- as.integer(d$Port == tr$port)
  d$cohortPost <- as.integer(d$quarter_index >= tr$treatmentIndex)
  model <- feols(
    log_total ~ cohortTreated:cohortPost | Port^Quarter_Number + Quarter,
    data = d, cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  tab <- coeftable(model); est <- unname(tab["cohortTreated:cohortPost", "Estimate"])
  se <- unname(tab["cohortTreated:cohortPost", "Std. Error"]); crit <- qt(.975, 2L)
  list(port = tr$port, treatmentQuarter = tr$treatmentQuarter,
       postQuarters = length(unique(d$quarter_index[d$cohortPost == 1L])),
       estimate = est, percent = 100 * (exp(est) - 1), standardError = se,
       percentCiLow = 100 * (exp(est - crit * se) - 1),
       percentCiHigh = 100 * (exp(est + crit * se) - 1),
       pValue = unname(tab["cohortTreated:cohortPost", "Pr(>|t|)"]))
})

staggered <- list(
  meta = list(method = "Seasonally adjusted stacked cohort DiD",
              controls = never_treated,
              fixedEffects = "stack × port × quarter-of-year + stack × calendar quarter",
              standardErrors = "clustered by original port (5 clusters)",
              observations = nrow(stacked)),
  summary = stacked_summary,
  cohorts = cohort_results
)

# ---------- Assemble output ----------
short_study <- fit_seasonal_study(make_panel(2019L), "All container-active ports · short sample")
full_study <- fit_seasonal_study(make_panel(2025L), "All container-active ports · full sample")
short_sensitivity <- fit_sensitivity(2019L)
full_sensitivity <- fit_sensitivity(2025L)

payload <- list(
  meta = list(
    title = "Piraeus container cargo: integrated causal-design dashboard",
    treatmentTransition = "2016 Q3 excluded",
    firstPostQuarter = "2016 Q4",
    lastCleanPreQuarter = "2016 Q2",
    source = "data/port-cargo-master-2014-2025.xlsx",
    generated = format(Sys.Date(), "%Y-%m-%d")
  ),
  headline = short_study,
  within = within,
  short = list(study = short_study, sensitivity = short_sensitivity),
  full = list(study = full_study, sensitivity = full_sensitivity),
  staggered = staggered
)

dir.create(dirname(json_path), showWarnings = FALSE, recursive = TRUE)
write_json(payload, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 15)

rows <- list()
add_row <- function(section, model, result) {
  rows[[length(rows) + 1L]] <<- data.frame(
    section = section, model = model, estimate = result$estimate,
    percent = result$percent, standard_error = result$standardError,
    percent_ci_low = result$percentCiLow, percent_ci_high = result$percentCiHigh,
    p_value = result$pValue
  )
}
add_row("short", "seasonal static DiD", short_study$did)
add_row("short", "average post vs clean pre", short_study$averagePre)
add_row("short", "average post vs 2016 Q2", short_study$q2)
add_row("full", "seasonal static DiD", full_study$did)
add_row("full", "average post vs clean pre", full_study$averagePre)
add_row("full", "average post vs 2016 Q2", full_study$q2)
for (x in within$windows) add_row(paste0("within ", x$label), "seasonally adjusted container vs non-container", x$seasonalDid)
add_row("staggered", "pooled stacked cohort DiD", staggered$summary)
write.csv(do.call(rbind, rows), csv_path, row.names = FALSE)

template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
json_text <- toJSON(payload, auto_unbox = TRUE, na = "null", digits = 15)
html <- sub("__INTEGRATED_CONTAINER_DATA__", json_text, template, fixed = TRUE)
writeLines(html, output_path, useBytes = TRUE)

cat("Created:", output_path, "\n")
cat(sprintf("Headline seasonal DiD: %+.1f%% (p = %.3f)\n", short_study$did$percent, short_study$did$pValue))
cat(sprintf("Short average post: clean pre %+.1f%%; 2016 Q2 %+.1f%%\n",
            short_study$averagePre$percent, short_study$q2$percent))
cat(sprintf("Within-Piraeus: short %+.1f%%; full %+.1f%%\n",
            within$windows[[1]]$seasonalDid$percent, within$windows[[2]]$seasonalDid$percent))
cat(sprintf("Full all-port seasonal DiD: %+.1f%% (p = %.3f)\n",
            full_study$did$percent, full_study$did$pValue))
cat(sprintf("Stacked staggered seasonal DiD: %+.1f%% (p = %.3f)\n",
            staggered$summary$percent, staggered$summary$pValue))
