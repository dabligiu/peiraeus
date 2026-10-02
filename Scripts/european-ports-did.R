# European-port difference-in-differences analysis of Piraeus container cargo.
#
# Source: data/Eurpean Ports Container Data.xlsx
# Units in the workbook match thousand tonnes when compared with the Greek
# source panel. The regression outcome is log container cargo, so rescaling the
# level unit would not change any treatment coefficient.
#
# Timing: the August 2016 transfer makes 2016 Q3 a transition quarter. It is
# excluded. 2016 Q2 is the last clean pre-treatment quarter and 2016 Q4 is the
# first full post-treatment quarter.

suppressPackageStartupMessages({
  library(readxl)
  library(fixest)
  library(jsonlite)
  library(quadprog)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "Eurpean Ports Container Data.xlsx")
template_path <- file.path(project_dir, "Scripts", "european-ports-did-template.html")
output_path <- file.path(project_dir, "piraeus-european-ports-did.html")
json_path <- file.path(project_dir, ".traffic-work", "european-ports-did-results.json")
csv_path <- file.path(project_dir, ".traffic-work", "european-ports-did-estimates.csv")

TREATMENT_YEAR <- 2016L
TREATMENT_QUARTER <- 4L
TRANSITION_LABEL <- "2016-Q3"
Q2_LABEL <- "2016-Q2"

wide <- as.data.frame(read_excel(data_path, sheet = "Sheet1", .name_repair = "minimal"))
quarter_columns <- names(wide)[grepl("^[0-9]{4}-Q[1-4]$", names(wide))]
if (length(quarter_columns) != 84L) stop("Expected 84 quarterly columns from 2005 Q1 to 2025 Q4")
if (!"Piraeus" %in% wide[[1]]) stop("Piraeus row not found")

panel <- do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
  labels <- quarter_columns
  data.frame(
    Port = as.character(wide[i, 1]),
    Quarter = labels,
    Year = as.integer(substr(labels, 1, 4)),
    Quarter_Number = as.integer(substr(labels, 7, 7)),
    container = as.numeric(unlist(wide[i, quarter_columns], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}))
if (any(!is.finite(panel$container)) || any(panel$container <= 0)) stop("All port-quarter outcomes must be positive and observed")
panel$quarter_index <- (panel$Year - 2005L) * 4L + panel$Quarter_Number - 1L
treatment_index <- (TREATMENT_YEAR - 2005L) * 4L + TREATMENT_QUARTER - 1L
panel$event_time <- panel$quarter_index - treatment_index
panel$treated <- as.integer(panel$Port == "Piraeus")
panel$post <- as.integer(panel$event_time >= 0L)
panel$log_total <- log(panel$container)
panel <- panel[panel$Quarter != TRANSITION_LABEL,]
panel <- panel[order(panel$Port, panel$quarter_index),]
panel$Quarter <- factor(panel$Quarter, levels = quarter_columns[quarter_columns != TRANSITION_LABEL])

all_ports <- sort(unique(panel$Port))
controls <- setdiff(all_ports, "Piraeus")
if (length(all_ports) != 8L || length(controls) != 7L) stop("Expected Piraeus and seven control ports")

donor_pools <- list(
  "All European controls" = controls,
  "Mediterranean controls" = c("Algeciras", "Valencia", "Gioia Tauro"),
  "Northern-range controls" = c("Antwerpen", "Bremerhaven", "Hamburg", "Rotterdam")
)

event_term <- function(k) paste0("event_time::", k, ":treated")

contrast_result <- function(cvec, beta, covariance, df) {
  estimate <- sum(cvec * beta)
  se <- sqrt(max(as.numeric(t(cvec) %*% covariance %*% cvec), 0))
  critical <- qt(.975, df)
  list(
    estimate = estimate, percent = 100 * (exp(estimate) - 1), standardError = se,
    percentCiLow = 100 * (exp(estimate - critical * se) - 1),
    percentCiHigh = 100 * (exp(estimate + critical * se) - 1),
    pValue = 2 * pt(-abs(estimate / se), df)
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

make_sample <- function(end_year, donors = controls) {
  d <- panel[panel$Year <= end_year & panel$Port %in% c("Piraeus", donors),]
  d$Quarter <- droplevels(d$Quarter)
  d
}

make_window_sample <- function(start_year, end_year, donors = controls) {
  d <- panel[
    panel$Year >= start_year & panel$Year <= end_year &
      panel$Port %in% c("Piraeus", donors),
  ]
  d$Quarter <- droplevels(d$Quarter)
  d
}

fit_study <- function(d, label) {
  event_times <- sort(unique(d$event_time))
  qnum <- ((event_times + treatment_index) %% 4L) + 1L
  base_times <- vapply(1:4, function(q) as.integer(min(event_times[qnum == q])), integer(1))
  clusters <- length(unique(d$Port)); df <- clusters - 1L

  model <- feols(
    log_total ~ i(event_time, treated, ref = base_times) | Port^Quarter_Number + Quarter,
    data = d, cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  beta <- coef(model); covariance <- vcov(model); names_beta <- names(beta)
  selector <- matrix(0, length(event_times), length(beta), dimnames = list(as.character(event_times), names_beta))
  for (i in seq_along(event_times)) {
    term <- event_term(event_times[i])
    if (term %in% names_beta) selector[i, term] <- 1
  }
  seasonal <- selector
  for (i in seq_along(event_times)) {
    pre_same_season <- which(event_times < 0L & qnum == qnum[i])
    seasonal[i,] <- selector[i,] - colMeans(selector[pre_same_season,, drop = FALSE])
  }
  q2_event <- (2016L - 2005L) * 4L + 2L - 1L - treatment_index
  q2_row <- match(q2_event, event_times)
  q2 <- sweep(seasonal, 2, seasonal[q2_row,], "-")
  estimates <- as.numeric(seasonal %*% beta); V_event <- seasonal %*% covariance %*% t(seasonal)
  q2_estimates <- as.numeric(q2 %*% beta); V_q2 <- q2 %*% covariance %*% t(q2)
  critical <- qt(.975, df)
  post_rows <- which(event_times >= 0L); pre_rows <- which(event_times < 0L)

  static_model <- feols(
    log_total ~ treated:post | Port^Quarter_Number + Quarter,
    data = d, cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  static_table <- coeftable(static_model)
  est <- unname(static_table["treated:post", "Estimate"])
  se <- unname(static_table["treated:post", "Std. Error"])
  static <- list(
    estimate = est, percent = 100 * (exp(est) - 1), standardError = se,
    percentCiLow = 100 * (exp(est - critical * se) - 1),
    percentCiHigh = 100 * (exp(est + critical * se) - 1),
    pValue = unname(static_table["treated:post", "Pr(>|t|)"])
  )

  rows <- lapply(seq_along(event_times), function(i) list(
    quarter = sprintf("%d Q%d", 2005L + (event_times[i] + treatment_index) %/% 4L,
                      (event_times[i] + treatment_index) %% 4L + 1L),
    eventTime = event_times[i], estimate = estimates[i],
    standardError = sqrt(max(V_event[i,i], 0)),
    ciLow = estimates[i] - critical * sqrt(max(V_event[i,i], 0)),
    ciHigh = estimates[i] + critical * sqrt(max(V_event[i,i], 0)),
    q2Estimate = q2_estimates[i], q2StandardError = sqrt(max(V_q2[i,i], 0)),
    q2CiLow = q2_estimates[i] - critical * sqrt(max(V_q2[i,i], 0)),
    q2CiHigh = q2_estimates[i] + critical * sqrt(max(V_q2[i,i], 0))
  ))

  list(
    label = label,
    meta = list(ports = clusters, controls = clusters - 1L,
                quarters = length(unique(d$quarter_index)),
                preQuarters = length(unique(d$quarter_index[d$event_time < 0L])),
                observations = nrow(d),
                fixedEffects = "port × quarter-of-year + calendar quarter",
                standardErrors = "port-clustered CR1"),
    did = static,
    averagePre = contrast_result(colMeans(seasonal[post_rows,, drop = FALSE]), beta, covariance, df),
    q2 = contrast_result(colMeans(q2[post_rows,, drop = FALSE]), beta, covariance, df),
    pretrend = joint_test(seasonal[pre_rows,, drop = FALSE], beta, covariance, df),
    events = rows
  )
}

fit_pool_sensitivity <- function(end_year) {
  lapply(names(donor_pools), function(name) {
    study <- fit_study(make_sample(end_year, donor_pools[[name]]), name)
    list(name = name, donors = donor_pools[[name]], did = study$did,
         averagePre = study$averagePre, q2 = study$q2, pretrend = study$pretrend)
  })
}

fit_leave_one_out <- function(end_year) {
  lapply(controls, function(omitted) {
    donors <- setdiff(controls, omitted)
    study <- fit_study(make_sample(end_year, donors), paste("Without", omitted))
    list(omitted = omitted, controls = length(donors), did = study$did,
         averagePre = study$averagePre, q2 = study$q2, pretrend = study$pretrend)
  })
}

fit_start_year_sensitivity <- function(end_year) {
  lapply(c(2005L, 2010L, 2012L, 2014L), function(start_year) {
    study <- fit_study(
      make_window_sample(start_year, end_year),
      sprintf("%d Q1–%d Q4", start_year, end_year)
    )
    list(startYear = start_year, label = study$label, meta = study$meta,
         did = study$did, averagePre = study$averagePre, q2 = study$q2,
         pretrend = study$pretrend)
  })
}

# Indexed level series for the descriptive chart. Every port is indexed to its
# own geometric mean over the clean pre-treatment period.
pre_mask <- panel$event_time < 0L
pre_means <- tapply(panel$log_total[pre_mask], panel$Port[pre_mask], mean)
panel$index <- 100 * exp(panel$log_total - pre_means[panel$Port])
level_series <- lapply(sort(unique(panel$quarter_index)), function(qi) {
  d <- panel[panel$quarter_index == qi,]
  list(
    quarter = as.character(d$Quarter[1]), year = d$Year[1], quarterNumber = d$Quarter_Number[1],
    eventTime = d$event_time[1],
    piraeus = d$index[d$Port == "Piraeus"],
    allControls = exp(mean(log(d$index[d$Port != "Piraeus"]))),
    mediterranean = exp(mean(log(d$index[d$Port %in% donor_pools[["Mediterranean controls"]]]))),
    northern = exp(mean(log(d$index[d$Port %in% donor_pools[["Northern-range controls"]]])))
  )
})

# ---------- Normalized synthetic control ----------
# Port scales differ substantially, and Piraeus lies below the convex hull of
# the largest European ports in levels for much of the pre-period. We therefore
# match seasonally adjusted proportional paths: each port's log series is
# stripped of its pre-treatment quarter-of-year effect and centered on its own
# clean pre-period geometric mean. The resulting counterfactual concerns growth
# relative to the port's own pre-treatment level, consistent with the log DiD.
synthetic_panel <- panel[panel$Year >= 2014L & panel$Year <= 2019L,]
synthetic_panel$seasonal_log <- NA_real_
for (port in all_ports) {
  rows_port <- synthetic_panel$Port == port
  pre_port <- rows_port & synthetic_panel$event_time < 0L
  overall_pre <- mean(synthetic_panel$log_total[pre_port])
  seasonal_effect <- vapply(1:4, function(q) {
    mean(synthetic_panel$log_total[pre_port & synthetic_panel$Quarter_Number == q]) - overall_pre
  }, numeric(1))
  adjusted <- synthetic_panel$log_total[rows_port] - seasonal_effect[synthetic_panel$Quarter_Number[rows_port]]
  synthetic_panel$seasonal_log[rows_port] <- adjusted - mean(adjusted[synthetic_panel$event_time[rows_port] < 0L])
}

quarter_ids <- sort(unique(synthetic_panel$quarter_index))
synthetic_ports <- all_ports
Y <- sapply(synthetic_ports, function(port) {
  d <- synthetic_panel[synthetic_panel$Port == port,]
  d$seasonal_log[match(quarter_ids, d$quarter_index)]
})
colnames(Y) <- synthetic_ports
pre_rows_sc <- quarter_ids < treatment_index

solve_synthetic_weights <- function(target, donors, pre_rows = pre_rows_sc) {
  X <- Y[pre_rows, donors, drop = FALSE]
  y <- Y[pre_rows, target]
  donor_count <- ncol(X)
  Dmat <- 2 * crossprod(X) + diag(1e-8, donor_count)
  dvec <- 2 * as.numeric(crossprod(X, y))
  Amat <- cbind(rep(1, donor_count), diag(donor_count))
  solution <- solve.QP(Dmat, dvec, Amat, c(1, rep(0, donor_count)), meq = 1)
  weights <- pmax(solution$solution, 0)
  weights <- weights / sum(weights)
  setNames(weights, donors)
}

synthetic_summary <- function(gap, pre_rmspe, end_year) {
  years <- 2005L + quarter_ids %/% 4L
  post <- quarter_ids >= treatment_index & years <= end_year
  average_gap <- mean(gap[post])
  post_rmspe <- sqrt(mean(gap[post]^2))
  list(
    endYear = end_year,
    postQuarters = sum(post),
    averageLogGap = average_gap,
    averagePercentGap = 100 * (exp(average_gap) - 1),
    preRmspe = pre_rmspe,
    postRmspe = post_rmspe,
    rmspeRatio = post_rmspe / pre_rmspe
  )
}

sc_weights <- solve_synthetic_weights("Piraeus", controls)
synthetic_path <- as.numeric(Y[, names(sc_weights), drop = FALSE] %*% sc_weights)
piraeus_path <- Y[, "Piraeus"]
sc_gap <- piraeus_path - synthetic_path
sc_pre_rmspe <- sqrt(mean(sc_gap[pre_rows_sc]^2))
sc_estimate <- synthetic_summary(sc_gap, sc_pre_rmspe, 2019L)

# Standard in-space reassignments: each control port is treated in turn and
# matched to all other ports. With only seven placebo ports, rank-based p-values
# cannot be smaller than 1/8 = 0.125.
placebos <- lapply(controls, function(placebo_port) {
  placebo_donors <- setdiff(all_ports, placebo_port)
  weights <- solve_synthetic_weights(placebo_port, placebo_donors)
  path <- as.numeric(Y[, placebo_donors, drop = FALSE] %*% weights)
  gap <- Y[, placebo_port] - path
  pre_rmspe <- sqrt(mean(gap[pre_rows_sc]^2))
  list(
    port = placebo_port,
    preRmspe = pre_rmspe,
    estimate = synthetic_summary(gap, pre_rmspe, 2019L)
  )
})

placebo_rank <- function(treated_summary) {
  placebo_gaps <- vapply(placebos, function(x) x$estimate$averagePercentGap, numeric(1))
  placebo_ratios <- vapply(placebos, function(x) x$estimate$rmspeRatio, numeric(1))
  list(
    gapRankP = (1 + sum(abs(placebo_gaps) >= abs(treated_summary$averagePercentGap))) /
      (1 + length(placebo_gaps)),
    rmspeRankP = (1 + sum(placebo_ratios >= treated_summary$rmspeRatio)) /
      (1 + length(placebo_ratios))
  )
}

sc_series <- lapply(seq_along(quarter_ids), function(i) {
  d <- synthetic_panel[synthetic_panel$quarter_index == quarter_ids[i],]
  list(
    quarter = as.character(d$Quarter[1]), year = d$Year[1],
    quarterNumber = d$Quarter_Number[1], eventTime = d$event_time[1],
    piraeusIndex = 100 * exp(piraeus_path[i]),
    syntheticIndex = 100 * exp(synthetic_path[i]),
    logGap = sc_gap[i], percentGap = 100 * (exp(sc_gap[i]) - 1)
  )
})

synthetic_control <- list(
  meta = list(
    method = "Convex synthetic control on seasonally adjusted, pre-period-centered log cargo",
    donors = controls,
    estimationWindow = "2014 Q1–2019 Q4",
    prePeriod = "2014 Q1–2016 Q2; 2016 Q3 excluded",
    constraints = "Non-negative weights summing to one",
    normalization = "Each port's clean pre-treatment geometric mean = 100"
  ),
  weights = lapply(names(sc_weights), function(port) list(port = port, weight = unname(sc_weights[port]))),
  series = sc_series,
  estimate = c(sc_estimate, placebo_rank(sc_estimate)),
  placebos = placebos
)

short <- fit_study(make_sample(2019L), "Pre-COVID sample: 2005 Q1–2019 Q4")
full <- fit_study(make_sample(2025L), "Full sample: 2005 Q1–2025 Q4")

payload <- list(
  meta = list(
    title = "European ports DiD: Piraeus container cargo",
    source = "data/Eurpean Ports Container Data.xlsx",
    units = "thousand tonnes (inferred by exact scale match to the Greek source panel)",
    treatmentTransition = "2016 Q3 excluded",
    lastCleanPre = "2016 Q2", firstPost = "2016 Q4",
    ports = all_ports
  ),
  levels = level_series,
  short = list(study = short, pools = fit_pool_sensitivity(2019L),
               leaveOneOut = fit_leave_one_out(2019L),
               startYears = fit_start_year_sensitivity(2019L)),
  full = list(study = full, pools = fit_pool_sensitivity(2025L),
              leaveOneOut = fit_leave_one_out(2025L),
              startYears = fit_start_year_sensitivity(2025L)),
  syntheticControl = synthetic_control
)

dir.create(dirname(json_path), showWarnings = FALSE, recursive = TRUE)
write_json(payload, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 15)

rows <- list()
add_result <- function(window, model, result) {
  rows[[length(rows)+1L]] <<- data.frame(
    window = window, model = model, estimate = result$estimate, percent = result$percent,
    standard_error = result$standardError, percent_ci_low = result$percentCiLow,
    percent_ci_high = result$percentCiHigh, p_value = result$pValue
  )
}
add_result("short", "Seasonal DiD", short$did)
add_result("short", "Average post vs clean pre", short$averagePre)
add_result("short", "Average post vs 2016 Q2", short$q2)
add_result("full", "Seasonal DiD", full$did)
add_result("full", "Average post vs clean pre", full$averagePre)
add_result("full", "Average post vs 2016 Q2", full$q2)
rows[[length(rows)+1L]] <- data.frame(
  window = "synthetic 2014 Q1–2019 Q4", model = "Normalized synthetic control average gap",
  estimate = sc_estimate$averageLogGap, percent = sc_estimate$averagePercentGap,
  standard_error = NA_real_, percent_ci_low = NA_real_, percent_ci_high = NA_real_,
  p_value = synthetic_control$estimate$gapRankP
)
write.csv(do.call(rbind, rows), csv_path, row.names = FALSE)

template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
html <- sub("__EUROPEAN_PORT_DATA__", toJSON(payload, auto_unbox = TRUE, na = "null", digits = 15), template, fixed = TRUE)
writeLines(html, output_path, useBytes = TRUE)

cat("Created:", output_path, "\n")
cat(sprintf("Short seasonal DiD: %+.1f%% (p = %.3f); event average vs clean pre %+.1f%%; vs Q2 %+.1f%%\n",
            short$did$percent, short$did$pValue, short$averagePre$percent, short$q2$percent))
cat(sprintf("Full seasonal DiD: %+.1f%% (p = %.3f); event average vs clean pre %+.1f%%; vs Q2 %+.1f%%\n",
            full$did$percent, full$did$pValue, full$averagePre$percent, full$q2$percent))
cat(sprintf("Synthetic control 2014 Q1–2019 Q4 average gap: %+.1f%% (rank p = %.3f); pre RMSPE %.3f; post/pre RMSPE %.3f\n",
            sc_estimate$averagePercentGap, synthetic_control$estimate$gapRankP,
            sc_pre_rmspe, sc_estimate$rmspeRatio))
