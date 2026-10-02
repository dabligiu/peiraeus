# Piraeus event study, 2014 Q1-2019 Q4
# Port and quarter fixed effects; 2016 Q3 (k = -1) is the omitted base.

suppressPackageStartupMessages({
  library(fixest)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-r.csv")
test_output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-tests-r.csv")
seasonal_output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-seasonal-r.csv")
seasonal_test_output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-seasonal-tests-r.csv")
honest_rm_output_path <- file.path(project_dir, ".traffic-work", "honest-did-relative-magnitude-r.csv")
honest_sd_output_path <- file.path(project_dir, ".traffic-work", "honest-did-smoothness-r.csv")
honest_summary_output_path <- file.path(project_dir, ".traffic-work", "honest-did-summary-r.csv")
source(file.path(dirname(script_path), "event-study-inference-helpers.R"))
outcome_config <- analysis_outcome_config()
treatment_config <- analysis_treatment_config()
output_path <- analysis_output_path(project_dir, "event-study-piraeus-r", outcome_config)
test_output_path <- analysis_output_path(project_dir, "event-study-piraeus-tests-r", outcome_config)
seasonal_output_path <- analysis_output_path(project_dir, "event-study-piraeus-seasonal-r", outcome_config)
seasonal_test_output_path <- analysis_output_path(project_dir, "event-study-piraeus-seasonal-tests-r", outcome_config)
honest_rm_output_path <- analysis_output_path(project_dir, "honest-did-relative-magnitude-r", outcome_config)
honest_sd_output_path <- analysis_output_path(project_dir, "honest-did-smoothness-r", outcome_config)
honest_summary_output_path <- analysis_output_path(project_dir, "honest-did-summary-r", outcome_config)

flow_totals <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)

panel <- as.data.frame(read_excel(data_path))
panel <- panel[
  panel$Year >= 2014 & panel$Year <= 2019 &
    panel$Port != "GR - aggregates extraction areas",
]
panel <- set_analysis_outcome(panel, outcome_config)

# Natural logs require positive cargo. Keep a balanced group of ports with a
# positive recorded total in every one of the 24 sample quarters.
positive_port <- ave(panel$total_cargo > 0, panel$Port, FUN = all)
panel <- panel[positive_port,]
panel$log_total <- log(panel$total_cargo)
panel$quarter_index <- (panel$Year - 2014) * 4 + panel$Quarter_Number - 1
panel$event_time <- panel$quarter_index - treatment_config$treatment_index
panel$treated <- as.integer(panel$Port == "Piraeus")
panel$post <- as.integer(panel$event_time >= 0)
panel$Quarter <- factor(panel$Quarter, levels = unique(panel$Quarter[order(panel$quarter_index)]))

model <- feols(
  log_total ~ i(event_time, treated, ref = -1) | Port + Quarter,
  data = panel,
  cluster = ~Port,
  ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
)

coef_table <- coeftable(model)
event_times <- treatment_config$pre_sample_event_times
results <- data.frame(
  quarter = sprintf("%d Q%d", 2014 + (0:23) %/% 4, (0:23) %% 4 + 1),
  event_time = event_times,
  estimate = NA_real_,
  std_error = NA_real_,
  ci_low = NA_real_,
  ci_high = NA_real_,
  p_value = NA_real_,
  omitted = event_times == -1
)

for (k in event_times[event_times != -1]) {
  term <- paste0("event_time::", k, ":treated")
  row <- coef_table[term,]
  position <- which(results$event_time == k)
  results$estimate[position] <- row["Estimate"]
  results$std_error[position] <- row["Std. Error"]
  results$ci_low[position] <- row["Estimate"] - qnorm(0.975) * row["Std. Error"]
  results$ci_high[position] <- row["Estimate"] + qnorm(0.975) * row["Std. Error"]
  results$p_value[position] <- 2 * pnorm(-abs(row["Estimate"] / row["Std. Error"]))
}
results$estimate[results$omitted] <- 0

# Robustness specification: allow every port to have a distinct recurring Q1,
# Q2, Q3, and Q4 pattern while retaining common calendar-quarter effects. The
# regression uses 2014 Q1-Q4 as arbitrary computational bases. Reported event
# coefficients are subsequently transformed into deviations from the average
# of all available untreated observations in the same quarter of the year.
seasonal_base_times <- event_times[1:4]
seasonal_model <- feols(
  log_total ~ i(event_time, treated, ref = seasonal_base_times) | Port^Quarter_Number + Quarter,
  data = panel,
  cluster = ~Port,
  ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
)
cluster_count <- length(unique(panel$Port))
seasonal_contrasts <- seasonal_average_event_study(
  seasonal_model,
  event_times = event_times,
  quarter_of_year = (0:23) %% 4 + 1,
  pre_event_times = treatment_config$all_pre_event_times,
  cluster_count = cluster_count
)
seasonal_results <- results
seasonal_results$estimate <- seasonal_contrasts$estimates
seasonal_results$std_error <- seasonal_contrasts$standard_errors
seasonal_results$ci_low <- seasonal_contrasts$ci_low
seasonal_results$ci_high <- seasonal_contrasts$ci_high
seasonal_results$p_value <- seasonal_contrasts$p_values
seasonal_results$omitted <- FALSE

test_statistics <- event_study_tests(model, cluster_count)
wild_did <- wild_cluster_did(panel)
tests <- data.frame(
  model = paste("Balanced all-port", outcome_config$short_title),
  donor_count = cluster_count - 1,
  observations = nrow(panel),
  did_log_points = as.numeric(wild_did["estimate"]),
  did_percent = 100 * (exp(as.numeric(wild_did["estimate"])) - 1),
  wild_bootstrap_p = as.numeric(wild_did["p_value"]),
  wild_bootstrap_replications = as.integer(wild_did["replications"]),
  wild_bootstrap_method = unname(wild_did["method"]),
  average_post_log_points = as.numeric(test_statistics$average_post["estimate"]),
  average_post_percent = 100 * (exp(as.numeric(test_statistics$average_post["estimate"])) - 1),
  average_post_p = as.numeric(test_statistics$average_post["p_value"]),
  joint_post_f = as.numeric(test_statistics$joint_post["f_statistic"]),
  joint_post_df1 = as.integer(test_statistics$joint_post["df1"]),
  joint_post_df2 = as.integer(test_statistics$joint_post["df2"]),
  joint_post_p = as.numeric(test_statistics$joint_post["p_value"]),
  pretrend_f = as.numeric(test_statistics$joint_pre["f_statistic"]),
  pretrend_df1 = as.integer(test_statistics$joint_pre["df1"]),
  pretrend_df2 = as.integer(test_statistics$joint_pre["df2"]),
  pretrend_p = as.numeric(test_statistics$joint_pre["p_value"])
)

# Rambachan-Roth HonestDiD sensitivity analysis for the mean of all 13
# post-treatment event coefficients. There are eleven observed pre-treatment
# quarters, but k = -1 is the normalization, leaving ten estimated leads.
# The relative-magnitude analysis bounds post-treatment violations by Mbar
# times the largest pre-treatment violation. The smoothness analysis bounds
# quarter-to-quarter changes in the slope of the counterfactual difference.
if (!requireNamespace("HonestDiD", quietly = TRUE)) {
  stop("Package 'HonestDiD' is required. Install it with install.packages('HonestDiD').")
}
honest_pre_times <- treatment_config$estimated_pre_event_times
honest_post_times <- treatment_config$pre_sample_post_event_times
honest_event_times <- c(honest_pre_times, honest_post_times)
honest_terms <- vapply(honest_event_times, event_term, character(1))
honest_beta <- coef(model)[honest_terms]
honest_sigma <- vcov(model)[honest_terms, honest_terms, drop = FALSE]
honest_l_vec <- rep(1 / length(honest_post_times), length(honest_post_times))
honest_mbar_grid <- c(0, 0.025, 0.05, 0.1, 0.25, 0.5)
honest_rm <- suppressWarnings(HonestDiD::createSensitivityResults_relativeMagnitudes(
  betahat = honest_beta,
  sigma = honest_sigma,
  numPrePeriods = length(honest_pre_times),
  numPostPeriods = length(honest_post_times),
  method = "C-LF",
  Mbarvec = honest_mbar_grid,
  l_vec = honest_l_vec,
  gridPoints = 1201,
  grid.lb = -5,
  grid.ub = 5,
  seed = 20260822
))
honest_rm <- as.data.frame(honest_rm)
honest_rm$includes_zero <- honest_rm$lb <= 0 & honest_rm$ub >= 0

honest_m_grid <- c(0, 0.025, 0.05, 0.1)
honest_sd <- do.call(rbind, lapply(honest_m_grid, function(m_value) {
  confidence_set <- suppressWarnings(HonestDiD::computeConditionalCS_DeltaSD(
    betahat = honest_beta,
    sigma = honest_sigma,
    numPrePeriods = length(honest_pre_times),
    numPostPeriods = length(honest_post_times),
    l_vec = honest_l_vec,
    M = m_value,
    alpha = 0.05,
    hybrid_flag = "LF",
    gridPoints = 4001,
    grid.lb = -10,
    grid.ub = 10,
    seed = 20260822
  ))
  accepted <- confidence_set$grid[confidence_set$accept == 1]
  if (length(accepted) == 0) stop("HonestDiD returned an empty smoothness confidence set")
  if (min(accepted) <= -9.999 || max(accepted) >= 9.999) {
    stop("HonestDiD smoothness confidence set reached the search-grid boundary")
  }
  data.frame(
    lb = min(accepted),
    ub = max(accepted),
    method = "C-LF",
    Delta = "DeltaSD",
    M = m_value,
    includes_zero = min(accepted) <= 0 & max(accepted) >= 0
  )
}))

honest_average <- test_statistics$average_post
first_rm_zero <- which(honest_rm$includes_zero)[1]
first_sd_zero <- which(honest_sd$includes_zero)[1]
honest_summary <- data.frame(
  estimand = paste("Mean of", length(honest_post_times), "post-treatment event coefficients"),
  estimate_log_points = as.numeric(honest_average["estimate"]),
  estimate_percent = 100 * (exp(as.numeric(honest_average["estimate"])) - 1),
  conventional_ci_low = as.numeric(honest_average["estimate"]) - qnorm(0.975) * as.numeric(honest_average["std_error"]),
  conventional_ci_high = as.numeric(honest_average["estimate"]) + qnorm(0.975) * as.numeric(honest_average["std_error"]),
  observed_pre_quarters = length(treatment_config$all_pre_event_times),
  estimated_pre_coefficients = length(honest_pre_times),
  post_coefficients = length(honest_post_times),
  first_rm_grid_including_zero = if (is.na(first_rm_zero)) NA_real_ else honest_rm$Mbar[first_rm_zero],
  last_rm_grid_excluding_zero = if (is.na(first_rm_zero) || first_rm_zero == 1) NA_real_ else honest_rm$Mbar[first_rm_zero - 1],
  first_sd_grid_including_zero = if (is.na(first_sd_zero)) NA_real_ else honest_sd$M[first_sd_zero],
  last_sd_grid_excluding_zero = if (is.na(first_sd_zero) || first_sd_zero == 1) NA_real_ else honest_sd$M[first_sd_zero - 1],
  honestdid_version = as.character(utils::packageVersion("HonestDiD"))
)

seasonal_did_model <- feols(
  log_total ~ treated:post | Port^Quarter_Number + Quarter,
  data = panel,
  cluster = ~Port,
  ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
)
seasonal_did_table <- coeftable(seasonal_did_model)
seasonal_did_estimate <- seasonal_did_table["treated:post", "Estimate"]
seasonal_did_se <- seasonal_did_table["treated:post", "Std. Error"]
seasonal_tests <- data.frame(
  model = paste("Port-by-quarter-of-year FE", outcome_config$short_title),
  donor_count = cluster_count - 1,
  observations = nrow(panel),
  did_log_points = seasonal_did_estimate,
  did_percent = 100 * (exp(seasonal_did_estimate) - 1),
  did_std_error = seasonal_did_se,
  did_clustered_p = seasonal_did_table["treated:post", "Pr(>|t|)"],
  average_post_log_points = as.numeric(seasonal_contrasts$average_post["estimate"]),
  average_post_percent = 100 * (exp(as.numeric(seasonal_contrasts$average_post["estimate"])) - 1),
  average_post_p = as.numeric(seasonal_contrasts$average_post["p_value"]),
  joint_post_p = as.numeric(seasonal_contrasts$joint_post["p_value"]),
  pretrend_p = as.numeric(seasonal_contrasts$joint_pre["p_value"])
)

write.csv(results, output_path, row.names = FALSE, na = "")
write.csv(tests, test_output_path, row.names = FALSE, na = "")
write.csv(seasonal_results, seasonal_output_path, row.names = FALSE, na = "")
write.csv(seasonal_tests, seasonal_test_output_path, row.names = FALSE, na = "")
write.csv(honest_rm, honest_rm_output_path, row.names = FALSE, na = "")
write.csv(honest_sd, honest_sd_output_path, row.names = FALSE, na = "")
write.csv(honest_summary, honest_summary_output_path, row.names = FALSE, na = "")
cat("Saved:", output_path, "\n")
cat("Saved:", test_output_path, "\n")
cat("Saved:", seasonal_output_path, "\n")
cat("Saved:", seasonal_test_output_path, "\n")
cat("Saved:", honest_rm_output_path, "\n")
cat("Saved:", honest_sd_output_path, "\n")
cat("Saved:", honest_summary_output_path, "\n")
cat("Ports:", length(unique(panel$Port)), " Controls:", length(unique(panel$Port)) - 1,
    " Observations:", nrow(panel), "\n")
print(results, row.names = FALSE)
cat("\nFour-test summary:\n")
print(tests, row.names = FALSE)
cat("\nPort-by-quarter-of-year robustness summary:\n")
print(seasonal_tests, row.names = FALSE)
cat("\nHonestDiD relative-magnitude sensitivity:\n")
print(honest_rm, row.names = FALSE)
cat("\nHonestDiD smoothness sensitivity:\n")
print(honest_sd, row.names = FALSE)
