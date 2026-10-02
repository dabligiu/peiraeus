# Piraeus event study, full sample 2014 Q1-2025 Q4.
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
output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-full-sample-r.csv")
test_output_path <- file.path(project_dir, ".traffic-work", "event-study-piraeus-full-sample-tests-r.csv")
source(file.path(dirname(script_path), "event-study-inference-helpers.R"))
outcome_config <- analysis_outcome_config()
treatment_config <- analysis_treatment_config()
output_path <- analysis_output_path(project_dir, "event-study-piraeus-full-sample-r", outcome_config)
test_output_path <- analysis_output_path(project_dir, "event-study-piraeus-full-sample-tests-r", outcome_config)

flow_totals <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)

panel <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
panel <- panel[
  panel$Year >= 2014 & panel$Year <= 2025 &
    panel$Port != "GR - aggregates extraction areas",
]
panel <- set_analysis_outcome(panel, outcome_config)

# Natural logs require positive cargo. Keep a balanced group of ports with a
# positive recorded total in every one of the 48 sample quarters.
positive_port <- ave(
  panel$total_cargo > 0,
  panel$Port,
  FUN = function(value) length(value) == 48L && all(value)
)
panel <- panel[positive_port,]
panel$log_total <- log(panel$total_cargo)
panel$quarter_index <- (panel$Year - 2014) * 4 + panel$Quarter_Number - 1
panel$event_time <- panel$quarter_index - treatment_config$treatment_index
panel$treated <- as.integer(panel$Port == "Piraeus")
panel$post <- as.integer(panel$event_time >= 0)
panel$Quarter <- factor(
  panel$Quarter,
  levels = sprintf("%d Q%d", 2014L + (0:47) %/% 4L, (0:47) %% 4L + 1L)
)

model <- feols(
  log_total ~ i(event_time, treated, ref = -1) | Port + Quarter,
  data = panel,
  cluster = ~Port,
  ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
)

coef_table <- coeftable(model)
event_times <- treatment_config$full_sample_event_times
results <- data.frame(
  quarter = sprintf("%d Q%d", 2014L + (0:47) %/% 4L, (0:47) %% 4L + 1L),
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

cluster_count <- length(unique(panel$Port))
test_statistics <- event_study_tests(
  model,
  cluster_count,
  pre_event_times = treatment_config$estimated_pre_event_times,
  post_event_times = treatment_config$full_sample_post_event_times
)
wild_did <- wild_cluster_did(panel)
tests <- data.frame(
  model = paste("Full-sample balanced all-port", outcome_config$short_title),
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

write.csv(results, output_path, row.names = FALSE, na = "")
write.csv(tests, test_output_path, row.names = FALSE, na = "")
cat("Saved:", output_path, "\n")
cat("Saved:", test_output_path, "\n")
cat(
  "Ports:", cluster_count,
  " Controls:", cluster_count - 1,
  " Observations:", nrow(panel), "\n"
)
print(results, row.names = FALSE)
cat("\nFour-test summary:\n")
print(tests, row.names = FALSE)
