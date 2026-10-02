# Piraeus donor-pool sensitivity analysis, 2014 Q1-2019 Q4
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
event_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pools-r.csv")
did_output <- file.path(project_dir, ".traffic-work", "did-donor-pools-r.csv")
test_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pool-tests-r.csv")
seasonal_event_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pools-seasonal-r.csv")
seasonal_did_output <- file.path(project_dir, ".traffic-work", "did-donor-pools-seasonal-r.csv")
seasonal_test_output <- file.path(project_dir, ".traffic-work", "event-study-donor-pool-seasonal-tests-r.csv")
source(file.path(dirname(script_path), "event-study-inference-helpers.R"))
outcome_config <- analysis_outcome_config()
treatment_config <- analysis_treatment_config()
event_output <- analysis_output_path(project_dir, "event-study-donor-pools-r", outcome_config)
did_output <- analysis_output_path(project_dir, "did-donor-pools-r", outcome_config)
test_output <- analysis_output_path(project_dir, "event-study-donor-pool-tests-r", outcome_config)
seasonal_event_output <- analysis_output_path(project_dir, "event-study-donor-pools-seasonal-r", outcome_config)
seasonal_did_output <- analysis_output_path(project_dir, "did-donor-pools-seasonal-r", outcome_config)
seasonal_test_output <- analysis_output_path(project_dir, "event-study-donor-pool-seasonal-tests-r", outcome_config)
honest_rm_output <- analysis_output_path(project_dir, "honest-did-donor-pools-relative-magnitude-r", outcome_config)
honest_sd_output <- analysis_output_path(project_dir, "honest-did-donor-pools-smoothness-r", outcome_config)
honest_summary_output <- analysis_output_path(project_dir, "honest-did-donor-pools-summary-r", outcome_config)

donor_pools <- analysis_donor_pools(outcome_config)

flow_totals <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)

panel <- as.data.frame(read_excel(data_path))
panel <- panel[panel$Year >= 2014 & panel$Year <= 2019,]
panel <- set_analysis_outcome(panel, outcome_config)
panel$log_total <- ifelse(panel$total_cargo > 0, log(panel$total_cargo), NA_real_)
panel$quarter_index <- (panel$Year - 2014) * 4 + panel$Quarter_Number - 1
panel$event_time <- panel$quarter_index - treatment_config$treatment_index
panel$treated <- as.integer(panel$Port == "Piraeus")
panel$post <- as.integer(panel$event_time >= 0)
panel$Quarter <- factor(panel$Quarter, levels = unique(panel$Quarter[order(panel$quarter_index)]))

event_rows <- list()
did_rows <- list()
test_rows <- list()
seasonal_event_rows <- list()
seasonal_did_rows <- list()
seasonal_test_rows <- list()
honest_rm_rows <- list()
honest_sd_rows <- list()
honest_summary_rows <- list()
event_models <- list()
event_test_statistics <- list()

if (!requireNamespace("HonestDiD", quietly = TRUE)) {
  stop("Package 'HonestDiD' is required. Install it with install.packages('HonestDiD').")
}

run_pool_honest_did <- function(event_model, test_statistics, pool_name, donors) {
  honest_pre_times <- treatment_config$estimated_pre_event_times
  honest_post_times <- treatment_config$pre_sample_post_event_times
  honest_event_times <- c(honest_pre_times, honest_post_times)
  honest_terms <- vapply(honest_event_times, event_term, character(1))
  honest_beta <- coef(event_model)[honest_terms]
  honest_sigma <- vcov(event_model)[honest_terms, honest_terms, drop = FALSE]
  honest_l_vec <- rep(1 / length(honest_post_times), length(honest_post_times))

  rm <- suppressWarnings(HonestDiD::createSensitivityResults_relativeMagnitudes(
    betahat = honest_beta,
    sigma = honest_sigma,
    numPrePeriods = length(honest_pre_times),
    numPostPeriods = length(honest_post_times),
    method = "C-LF",
    Mbarvec = c(0, 0.025, 0.05, 0.1, 0.25, 0.5),
    l_vec = honest_l_vec,
    gridPoints = 2401,
    grid.lb = -20,
    grid.ub = 20,
    seed = 20260822
  ))
  rm <- as.data.frame(rm)
  rm$includes_zero <- rm$lb <= 0 & rm$ub >= 0
  rm$pool <- pool_name
  rm$donors <- paste(donors, collapse = "; ")

  sd <- do.call(rbind, lapply(c(0, 0.025, 0.05, 0.1), function(m_value) {
    used_bound <- 40
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
      grid.lb = -used_bound,
      grid.ub = used_bound,
      seed = 20260822
    ))
    accepted <- confidence_set$grid[confidence_set$accept == 1]
    if (length(accepted) == 0) stop("HonestDiD returned an empty smoothness confidence set for ", pool_name)
    if (min(accepted) <= -used_bound + 1e-8 || max(accepted) >= used_bound - 1e-8) {
      stop("HonestDiD smoothness confidence set reached the search-grid boundary for ", pool_name)
    }
    data.frame(
      lb = min(accepted),
      ub = max(accepted),
      method = "C-LF",
      Delta = "DeltaSD",
      M = m_value,
      includes_zero = min(accepted) <= 0 & max(accepted) >= 0,
      pool = pool_name,
      donors = paste(donors, collapse = "; ")
    )
  }))

  average_post <- test_statistics$average_post
  first_rm_zero <- which(rm$includes_zero)[1]
  first_sd_zero <- which(sd$includes_zero)[1]
  summary <- data.frame(
    pool = pool_name,
    donors = paste(donors, collapse = "; "),
    donor_count = length(donors),
    estimand = paste("Mean of", length(honest_post_times), "post-treatment event coefficients"),
    estimate_log_points = as.numeric(average_post["estimate"]),
    estimate_percent = 100 * (exp(as.numeric(average_post["estimate"])) - 1),
    conventional_ci_low = as.numeric(average_post["estimate"]) - qnorm(0.975) * as.numeric(average_post["std_error"]),
    conventional_ci_high = as.numeric(average_post["estimate"]) + qnorm(0.975) * as.numeric(average_post["std_error"]),
    first_rm_grid_including_zero = if (is.na(first_rm_zero)) NA_real_ else rm$Mbar[first_rm_zero],
    last_rm_grid_excluding_zero = if (is.na(first_rm_zero) || first_rm_zero == 1) NA_real_ else rm$Mbar[first_rm_zero - 1],
    first_sd_grid_including_zero = if (is.na(first_sd_zero)) NA_real_ else sd$M[first_sd_zero],
    last_sd_grid_excluding_zero = if (is.na(first_sd_zero) || first_sd_zero == 1) NA_real_ else sd$M[first_sd_zero - 1],
    honestdid_version = as.character(utils::packageVersion("HonestDiD"))
  )
  list(relative_magnitude = rm, smoothness = sd, summary = summary)
}

for (pool_name in names(donor_pools)) {
  donors <- donor_pools[[pool_name]]
  required_ports <- c("Piraeus", donors)
  sample <- panel[panel$Port %in% required_ports,]

  counts <- table(sample$Port)
  if (!all(required_ports %in% names(counts)) || any(counts[required_ports] != 24)) {
    stop("The pool is not a balanced 24-quarter panel: ", pool_name)
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

  pool_result <- data.frame(
    pool = pool_name,
    donors = paste(donors, collapse = "; "),
    donor_count = length(donors),
    quarter = sprintf("%d Q%d", 2014 + (0:23) %/% 4, (0:23) %% 4 + 1),
    event_time = treatment_config$pre_sample_event_times,
    estimate = NA_real_,
    std_error = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = NA_real_,
    omitted = treatment_config$pre_sample_event_times == -1
  )
  for (k in treatment_config$pre_sample_event_times[treatment_config$pre_sample_event_times != -1]) {
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
  test_statistics <- event_study_tests(event_model, cluster_count)
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

  event_models[[pool_name]] <- event_model
  event_test_statistics[[pool_name]] <- test_statistics

  # Robustness specification with a separate recurring seasonal profile for
  # every port, in addition to common calendar-quarter fixed effects. The four
  # 2014 quarters are arbitrary computational bases; reported coefficients are
  # transformed to deviations from the same-season untreated-period average.
  seasonal_base_times <- treatment_config$pre_sample_event_times[1:4]
  seasonal_event_model <- feols(
    log_total ~ i(event_time, treated, ref = seasonal_base_times) | Port^Quarter_Number + Quarter,
    data = sample,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  seasonal_contrasts <- seasonal_average_event_study(
    seasonal_event_model,
    event_times = treatment_config$pre_sample_event_times,
    quarter_of_year = (0:23) %% 4 + 1,
    pre_event_times = treatment_config$all_pre_event_times,
    cluster_count = cluster_count
  )
  seasonal_pool_result <- pool_result
  seasonal_pool_result$estimate <- seasonal_contrasts$estimates
  seasonal_pool_result$std_error <- seasonal_contrasts$standard_errors
  seasonal_pool_result$ci_low <- seasonal_contrasts$ci_low
  seasonal_pool_result$ci_high <- seasonal_contrasts$ci_high
  seasonal_pool_result$p_value <- seasonal_contrasts$p_values
  seasonal_pool_result$omitted <- FALSE
  seasonal_event_rows[[pool_name]] <- seasonal_pool_result

  seasonal_did_model <- feols(
    log_total ~ treated:post | Port^Quarter_Number + Quarter,
    data = sample,
    cluster = ~Port,
    ssc = ssc(K.adj = TRUE, K.fixef = "full", G.adj = TRUE)
  )
  seasonal_did_table <- coeftable(seasonal_did_model)
  seasonal_did_estimate <- seasonal_did_table["treated:post", "Estimate"]
  seasonal_did_se <- seasonal_did_table["treated:post", "Std. Error"]
  seasonal_did_rows[[pool_name]] <- data.frame(
    pool = pool_name,
    donors = paste(donors, collapse = "; "),
    donor_count = length(donors),
    observations = nrow(sample),
    log_points = seasonal_did_estimate,
    percent_effect = 100 * (exp(seasonal_did_estimate) - 1),
    std_error = seasonal_did_se,
    ci_low = seasonal_did_estimate - qnorm(0.975) * seasonal_did_se,
    ci_high = seasonal_did_estimate + qnorm(0.975) * seasonal_did_se,
    clustered_p = seasonal_did_table["treated:post", "Pr(>|t|)"]
  )

  seasonal_test_rows[[pool_name]] <- data.frame(
    pool = pool_name,
    donor_count = length(donors),
    observations = nrow(sample),
    did_percent = 100 * (exp(seasonal_did_estimate) - 1),
    did_clustered_p = seasonal_did_table["treated:post", "Pr(>|t|)"],
    average_post_percent = 100 * (exp(as.numeric(seasonal_contrasts$average_post["estimate"])) - 1),
    average_post_p = as.numeric(seasonal_contrasts$average_post["p_value"]),
    joint_post_p = as.numeric(seasonal_contrasts$joint_post["p_value"]),
    pretrend_p = as.numeric(seasonal_contrasts$joint_pre["p_value"])
  )
}

# The three HonestDiD calculations are independent. On non-Windows systems,
# run them concurrently because the conditional confidence-set grid search is
# substantially slower than the underlying fixed-effects regressions.
honest_pool_names <- names(donor_pools)
honest_results <- parallel::mclapply(
  honest_pool_names,
  function(pool_name) run_pool_honest_did(
    event_models[[pool_name]],
    event_test_statistics[[pool_name]],
    pool_name,
    donor_pools[[pool_name]]
  ),
  mc.cores = min(3L, length(honest_pool_names)),
  mc.set.seed = FALSE
)
names(honest_results) <- honest_pool_names
for (pool_name in honest_pool_names) {
  honest_rm_rows[[pool_name]] <- honest_results[[pool_name]]$relative_magnitude
  honest_sd_rows[[pool_name]] <- honest_results[[pool_name]]$smoothness
  honest_summary_rows[[pool_name]] <- honest_results[[pool_name]]$summary
}

event_results <- do.call(rbind, event_rows)
did_results <- do.call(rbind, did_rows)
test_results <- do.call(rbind, test_rows)
seasonal_event_results <- do.call(rbind, seasonal_event_rows)
seasonal_did_results <- do.call(rbind, seasonal_did_rows)
seasonal_test_results <- do.call(rbind, seasonal_test_rows)
honest_rm_results <- do.call(rbind, honest_rm_rows)
honest_sd_results <- do.call(rbind, honest_sd_rows)
honest_summary_results <- do.call(rbind, honest_summary_rows)
rownames(event_results) <- NULL
rownames(did_results) <- NULL
rownames(test_results) <- NULL
rownames(seasonal_event_results) <- NULL
rownames(seasonal_did_results) <- NULL
rownames(seasonal_test_results) <- NULL
rownames(honest_rm_results) <- NULL
rownames(honest_sd_results) <- NULL
rownames(honest_summary_results) <- NULL
write.csv(event_results, event_output, row.names = FALSE, na = "")
write.csv(did_results, did_output, row.names = FALSE, na = "")
write.csv(test_results, test_output, row.names = FALSE, na = "")
write.csv(seasonal_event_results, seasonal_event_output, row.names = FALSE, na = "")
write.csv(seasonal_did_results, seasonal_did_output, row.names = FALSE, na = "")
write.csv(seasonal_test_results, seasonal_test_output, row.names = FALSE, na = "")
write.csv(honest_rm_results, honest_rm_output, row.names = FALSE, na = "")
write.csv(honest_sd_results, honest_sd_output, row.names = FALSE, na = "")
write.csv(honest_summary_results, honest_summary_output, row.names = FALSE, na = "")

cat("Saved:", event_output, "\n")
cat("Saved:", did_output, "\n\n")
cat("Saved:", test_output, "\n\n")
cat("Saved:", seasonal_event_output, "\n")
cat("Saved:", seasonal_did_output, "\n")
cat("Saved:", seasonal_test_output, "\n\n")
cat("Saved:", honest_rm_output, "\n")
cat("Saved:", honest_sd_output, "\n")
cat("Saved:", honest_summary_output, "\n\n")
print(did_results[, c("pool", "donor_count", "observations", "log_points", "percent_effect")], row.names = FALSE)
cat("\nFour-test summary:\n")
print(test_results[, c(
  "pool", "did_percent", "wild_bootstrap_p", "average_post_percent",
  "average_post_p", "joint_post_p", "pretrend_p"
)], row.names = FALSE)
cat("\nPort-by-quarter-of-year robustness summary:\n")
print(seasonal_test_results, row.names = FALSE)
cat("\nLargest donor pool:\n", paste(donor_pools[[length(donor_pools)]], collapse = ", "), "\n")
