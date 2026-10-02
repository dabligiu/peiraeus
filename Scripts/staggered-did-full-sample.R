# Staggered-adoption DiD for the full 2014 Q1-2025 Q4 port panel.
#
# Estimator: Sun and Abraham (2021) interaction-weighted event study, with
# never-treated ports as controls, port and quarter fixed effects, and
# port-clustered standard errors. Treatment begins in the quarter in which the
# ownership transfer or terminal operating concession became effective.

suppressPackageStartupMessages({
  library(fixest)
  library(jsonlite)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
output_dir <- file.path(project_dir, ".traffic-work")
source(file.path(dirname(script_path), "event-study-inference-helpers.R"))
outcome_config <- analysis_outcome_config()

# Quarter indices use 2014 Q1 = 0. These dates are deliberately visible and
# easy to edit because treatment definitions are substantive research choices.
treatments <- data.frame(
  port = c("Piraeus", "Thessaloniki", "Igoumenitsa", "Heraklio"),
  treatment_quarter = c("2016 Q4", "2018 Q2", "2024 Q1", "2024 Q4"),
  treatment_index = c(11L, 17L, 40L, 43L),
  event = c(
    "Transfer of 51% of PPA shares",
    "Transfer of 67% of ThPA shares",
    "Transfer of 67% of IPA shares",
    "Transfer of 67% of HPA shares"
  ),
  source = c(
    "https://hradf.com/en/piraeus-port-authority-s-a-ppa/",
    "https://hradf.com/wp-content/uploads/2021/11/hradf-thpa-share-transfer.pdf",
    "https://hradf.com/wp-content/uploads/2024/05/Asset-Development-Plan-%CE%91DP-of-HRADF-S.A.-29-Dec-2023.pdf",
    "https://growthfund.gr/en/hradf-acquisition-of-a-majority-stake-in-the-share-capital-of-heraklion-port-authority-hpa-s-a-by-the-consortium-grimaldi-euromed-s-p-a-minoan-lines-s-a-for-80-million-euros/"
  ),
  stringsAsFactors = FALSE
)

panel <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
panel <- panel[!grepl("^GR -", panel$Port),]
panel <- set_analysis_outcome(panel, outcome_config)

# Natural logs require positive outcomes. Retain the balanced set of ports with
# positive total cargo in every one of the 48 quarters.
positive_balanced <- ave(
  panel$total_cargo > 0,
  panel$Port,
  FUN = function(value) length(value) == 48L && all(value)
)
panel <- panel[positive_balanced,]
panel$log_total <- log(panel$total_cargo)
panel$quarter_index <- (panel$Year - 2014L) * 4L + panel$Quarter_Number - 1L
panel$Quarter <- factor(
  panel$Quarter,
  levels = sprintf("%d Q%d", 2014L + (0:47) %/% 4L, (0:47) %% 4L + 1L)
)

treatments <- treatments[treatments$port %in% unique(panel$Port),, drop = FALSE]
if (nrow(treatments) < 2) stop("Fewer than two treated ports remain in the positive balanced sample")
if (nrow(panel) != length(unique(panel$Port)) * 48L) {
  stop("The full-sample analysis panel is not balanced")
}

cohort_lookup <- setNames(treatments$treatment_index, treatments$port)
# fixest's sunab() convention uses a cohort outside the observed period for
# never-treated units.
panel$cohort <- ifelse(
  panel$Port %in% treatments$port,
  unname(cohort_lookup[panel$Port]),
  1000L
)

fit_model <- function(data) {
  feols(
    log_total ~ sunab(cohort, quarter_index, ref.p = -1) | Port + Quarter,
    data = data,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
}

model <- fit_model(panel)
cluster_count <- length(unique(panel$Port))
treated_count <- nrow(treatments)
control_count <- cluster_count - treated_count

att_table <- aggregate(model, "att")
att_estimate <- unname(att_table[1, "Estimate"])
att_se <- unname(att_table[1, "Std. Error"])
att_p <- unname(att_table[1, "Pr(>|t|)"])
critical_value <- qt(0.975, df = cluster_count - 1L)

period_table <- aggregate(model, "period")
period_values <- as.integer(sub("^.*::", "", rownames(period_table)))
event_window <- -12:12
event_rows <- lapply(event_window, function(k) {
  if (k == -1L) {
    return(data.frame(
      event_time = k, estimate = 0, std_error = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, p_value = NA_real_, contributing_cohorts = treated_count,
      omitted = TRUE
    ))
  }
  position <- match(k, period_values)
  if (is.na(position)) stop("Missing relative-period estimate: ", k)
  estimate <- period_table[position, "Estimate"]
  standard_error <- period_table[position, "Std. Error"]
  available <- sum(
    treatments$treatment_index + k >= 0L &
      treatments$treatment_index + k <= 47L
  )
  data.frame(
    event_time = k,
    estimate = estimate,
    std_error = standard_error,
    ci_low = estimate - critical_value * standard_error,
    ci_high = estimate + critical_value * standard_error,
    p_value = period_table[position, "Pr(>|t|)"],
    contributing_cohorts = available,
    omitted = FALSE
  )
})
event_results <- do.call(rbind, event_rows)

cohort_table <- aggregate(model, "cohort")
cohort_indices <- as.integer(sub("^cohort::", "", rownames(cohort_table)))
cohort_results <- treatments
cohort_results$estimate <- NA_real_
cohort_results$std_error <- NA_real_
cohort_results$ci_low <- NA_real_
cohort_results$ci_high <- NA_real_
cohort_results$p_value <- NA_real_
cohort_results$post_quarters <- 48L - cohort_results$treatment_index
for (row_index in seq_len(nrow(cohort_results))) {
  position <- match(cohort_results$treatment_index[row_index], cohort_indices)
  estimate <- cohort_table[position, "Estimate"]
  standard_error <- cohort_table[position, "Std. Error"]
  cohort_results$estimate[row_index] <- estimate
  cohort_results$std_error[row_index] <- standard_error
  cohort_results$ci_low[row_index] <- estimate - critical_value * standard_error
  cohort_results$ci_high[row_index] <- estimate + critical_value * standard_error
  cohort_results$p_value[row_index] <- cohort_table[position, "Pr(>|t|)"]
}
cohort_results$percent_effect <- 100 * (exp(cohort_results$estimate) - 1)

# Generalized Wald test for every cohort-specific coefficient in the displayed
# pre-treatment window. A generalized inverse is used because the clustered
# covariance can be rank-deficient when there are only five treated ports.
pre_pattern <- paste0(
  "^quarter_index::-(",
  paste(12:2, collapse = "|"),
  "):cohort::"
)
pre_terms <- names(model$coefficients)[grepl(pre_pattern, names(model$coefficients))]
pre_beta <- model$coefficients[pre_terms]
pre_vcov <- vcov(model)[pre_terms, pre_terms, drop = FALSE]
decomposition <- eigen((pre_vcov + t(pre_vcov)) / 2, symmetric = TRUE)
tolerance <- max(decomposition$values) * 1e-8
keep <- decomposition$values > tolerance
pre_rank <- sum(keep)
projected <- crossprod(decomposition$vectors[, keep, drop = FALSE], pre_beta)
pre_f <- sum((projected^2) / decomposition$values[keep]) / pre_rank
pretrend_p <- pf(pre_f, df1 = pre_rank, df2 = cluster_count - 1L, lower.tail = FALSE)

# Leave-one-treated-port-out estimates show whether the pooled ATT is driven by
# any single treated port. These are sensitivity estimates, not independent
# hypothesis tests.
leave_one_out <- lapply(treatments$port, function(omitted_port) {
  sample <- panel[panel$Port != omitted_port,]
  sensitivity_model <- fit_model(sample)
  sensitivity_att <- aggregate(sensitivity_model, "att")
  data.frame(
    omitted_port = omitted_port,
    estimate = sensitivity_att[1, "Estimate"],
    std_error = sensitivity_att[1, "Std. Error"],
    p_value = sensitivity_att[1, "Pr(>|t|)"],
    percent_effect = 100 * (exp(sensitivity_att[1, "Estimate"]) - 1)
  )
})
leave_one_out <- do.call(rbind, leave_one_out)

summary_results <- data.frame(
  estimator = "Sun-Abraham interaction-weighted",
  sample = paste("2014 Q1-2025 Q4", outcome_config$short_title),
  ports = cluster_count,
  treated_ports = treated_count,
  never_treated_ports = control_count,
  quarters = 48L,
  observations = nrow(panel),
  att_log_points = att_estimate,
  att_percent = 100 * (exp(att_estimate) - 1),
  std_error = att_se,
  ci_low = att_estimate - critical_value * att_se,
  ci_high = att_estimate + critical_value * att_se,
  p_value = att_p,
  displayed_pretrend_f = pre_f,
  displayed_pretrend_df1 = pre_rank,
  displayed_pretrend_df2 = cluster_count - 1L,
  displayed_pretrend_p = pretrend_p,
  stringsAsFactors = FALSE
)

write.csv(event_results, analysis_output_path(project_dir, "staggered-did-full-sample-event", outcome_config), row.names = FALSE, na = "")
write.csv(cohort_results, analysis_output_path(project_dir, "staggered-did-full-sample-cohorts", outcome_config), row.names = FALSE, na = "")
write.csv(summary_results, analysis_output_path(project_dir, "staggered-did-full-sample-summary", outcome_config), row.names = FALSE, na = "")
write.csv(leave_one_out, analysis_output_path(project_dir, "staggered-did-full-sample-leave-one-out", outcome_config), row.names = FALSE, na = "")

json_result <- list(
  meta = list(
    estimator = summary_results$estimator,
    sample = summary_results$sample,
    outcome = outcome_config$log_label,
    ports = cluster_count,
    treatedPorts = treated_count,
    neverTreatedPorts = control_count,
    quarters = 48L,
    observations = nrow(panel),
    standardErrors = "Port-clustered CR1",
    sampleRule = paste("Balanced ports with", outcome_config$sample_label, "in all 48 quarters"),
    eventWindow = "k = -12 to +12; k = -1 omitted",
    comparisonGroup = "Never-treated ports",
    source = "data/port-cargo-master-2014-2025.xlsx"
  ),
  summary = list(
    estimate = att_estimate,
    percentEffect = summary_results$att_percent,
    standardError = att_se,
    ciLow = summary_results$ci_low,
    ciHigh = summary_results$ci_high,
    pValue = att_p,
    pretrendF = pre_f,
    pretrendDf1 = pre_rank,
    pretrendDf2 = cluster_count - 1L,
    pretrendP = pretrend_p
  ),
  treatments = lapply(seq_len(nrow(treatments)), function(i) as.list(treatments[i,])),
  cohorts = lapply(seq_len(nrow(cohort_results)), function(i) as.list(cohort_results[i,])),
  events = lapply(seq_len(nrow(event_results)), function(i) as.list(event_results[i,])),
  leaveOneOut = lapply(seq_len(nrow(leave_one_out)), function(i) as.list(leave_one_out[i,]))
)
write_json(
  json_result,
  analysis_output_path(project_dir, "staggered-did-full-sample", outcome_config, extension = "json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null",
  digits = 15
)

cat("Saved full-sample staggered-DiD outputs in", output_dir, "\n")
print(summary_results, row.names = FALSE)
cat("\nCohort effects:\n")
print(cohort_results[, c("port", "treatment_quarter", "post_quarters", "percent_effect", "p_value")], row.names = FALSE)
cat("\nLeave-one-treated-port-out ATT sensitivity:\n")
print(leave_one_out[, c("omitted_port", "percent_effect", "p_value")], row.names = FALSE)
