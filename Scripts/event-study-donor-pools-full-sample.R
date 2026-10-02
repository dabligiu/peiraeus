# Piraeus donor-pool sensitivity analysis, full sample 2014 Q1-2025 Q4.
# Runs identical event studies for three explicit comparison-port pools and a
# simple pre/post difference-in-differences summary for each pool.

suppressPackageStartupMessages({
  library(fixest)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
event_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pools-full-sample-r.csv")
did_output <- file.path(project_dir, ".traffic-work", "did-donor-pools-full-sample-r.csv")
test_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pool-full-sample-tests-r.csv")
source(file.path(dirname(script_path), "event-study-inference-helpers.R"))
outcome_config <- analysis_outcome_config()
treatment_config <- analysis_treatment_config()
event_output <- analysis_output_path(project_dir, "event-study-donor-pools-full-sample-r", outcome_config)
did_output <- analysis_output_path(project_dir, "did-donor-pools-full-sample-r", outcome_config)
test_output <- analysis_output_path(project_dir, "event-study-donor-pool-full-sample-tests-r", outcome_config)

donor_pools <- analysis_donor_pools(outcome_config)

flow_totals <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)

panel <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
panel <- panel[panel$Year >= 2014 & panel$Year <= 2025,]
panel <- set_analysis_outcome(panel, outcome_config)
panel$log_total <- ifelse(panel$total_cargo > 0, log(panel$total_cargo), NA_real_)
panel$quarter_index <- (panel$Year - 2014) * 4 + panel$Quarter_Number - 1
panel$event_time <- panel$quarter_index - treatment_config$treatment_index
panel$treated <- as.integer(panel$Port == "Piraeus")
panel$post <- as.integer(panel$event_time >= 0)
panel$Quarter <- factor(
  panel$Quarter,
  levels = sprintf("%d Q%d", 2014L + (0:47) %/% 4L, (0:47) %% 4L + 1L)
)

event_rows <- list()
did_rows <- list()
test_rows <- list()

for (pool_name in names(donor_pools)) {
  donors <- donor_pools[[pool_name]]
  required_ports <- c("Piraeus", donors)
  sample <- panel[panel$Port %in% required_ports,]

  counts <- table(sample$Port)
  if (!all(required_ports %in% names(counts)) || any(counts[required_ports] != 48)) {
    stop("The pool is not a balanced 48-quarter panel: ", pool_name)
  }
  if (any(!is.finite(sample$log_total))) {
    stop("A pool contains non-positive total cargo: ", pool_name)
  }

  event_model <- feols(
    log_total ~ i(event_time, treated, ref = -1) | Port + Quarter,
    data = sample,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  event_table <- coeftable(event_model)

  event_times <- treatment_config$full_sample_event_times
  pool_result <- data.frame(
    pool = pool_name,
    donors = paste(donors, collapse = "; "),
    donor_count = length(donors),
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
    row <- event_table[term,]
    position <- which(pool_result$event_time == k)
    pool_result$estimate[position] <- row["Estimate"]
    pool_result$std_error[position] <- row["Std. Error"]
    pool_result$ci_low[position] <- row["Estimate"] - qnorm(0.975) * row["Std. Error"]
    pool_result$ci_high[position] <- row["Estimate"] + qnorm(0.975) * row["Std. Error"]
    pool_result$p_value[position] <- 2 * pnorm(-abs(row["Estimate"] / row["Std. Error"]))
  }
  pool_result$estimate[pool_result$omitted] <- 0
  event_rows[[pool_name]] <- pool_result

  did_model <- feols(
    log_total ~ treated:post | Port + Quarter,
    data = sample,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  did_table <- coeftable(did_model)
  did_term <- "treated:post"
  did_estimate <- did_table[did_term, "Estimate"]
  did_se <- did_table[did_term, "Std. Error"]
  did_rows[[pool_name]] <- data.frame(
    pool = pool_name,
    donors = paste(donors, collapse = "; "),
    donor_count = length(donors),
    observations = nrow(sample),
    log_points = did_estimate,
    percent_effect = 100 * (exp(did_estimate) - 1),
    std_error = did_se,
    ci_low = did_estimate - qnorm(0.975) * did_se,
    ci_high = did_estimate + qnorm(0.975) * did_se
  )

  cluster_count <- length(unique(sample$Port))
  test_statistics <- event_study_tests(
    event_model,
    cluster_count,
    pre_event_times = treatment_config$estimated_pre_event_times,
    post_event_times = treatment_config$full_sample_post_event_times
  )
  wild_did <- wild_cluster_did(sample)
  test_rows[[pool_name]] <- data.frame(
    pool = pool_name,
    donor_count = length(donors),
    observations = nrow(sample),
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
}

event_results <- do.call(rbind, event_rows)
did_results <- do.call(rbind, did_rows)
test_results <- do.call(rbind, test_rows)
rownames(event_results) <- NULL
rownames(did_results) <- NULL
rownames(test_results) <- NULL
write.csv(event_results, event_output, row.names = FALSE, na = "")
write.csv(did_results, did_output, row.names = FALSE, na = "")
write.csv(test_results, test_output, row.names = FALSE, na = "")

cat("Saved:", event_output, "\n")
cat("Saved:", did_output, "\n")
cat("Saved:", test_output, "\n\n")
print(
  did_results[, c("pool", "donor_count", "observations", "log_points", "percent_effect")],
  row.names = FALSE
)
cat("\nFour-test summary:\n")
print(
  test_results[, c(
    "pool", "did_percent", "wild_bootstrap_p", "average_post_percent",
    "average_post_p", "joint_post_p", "pretrend_p"
  )],
  row.names = FALSE
)
