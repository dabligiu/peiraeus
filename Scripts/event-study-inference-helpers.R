# Shared inference helpers for the Piraeus event-study scripts.

analysis_outcome_config <- function() {
  key <- Sys.getenv("PIRAEUS_OUTCOME", unset = "total")
  valid <- c("total", "containers", "non-containers")
  if (!key %in% valid) stop("PIRAEUS_OUTCOME must be one of: ", paste(valid, collapse = ", "))
  switch(
    key,
    "total" = list(
      key = key, suffix = "", title = "Total cargo", short_title = "total cargo",
      log_label = "ln(total cargo tonnes)", sample_label = "positive total cargo"
    ),
    "containers" = list(
      key = key, suffix = "-containers", title = "Container cargo", short_title = "container cargo",
      log_label = "ln(container cargo tonnes)", sample_label = "positive container cargo"
    ),
    "non-containers" = list(
      key = key, suffix = "-non-containers", title = "Non-container cargo", short_title = "non-container cargo",
      log_label = "ln(non-container cargo tonnes)", sample_label = "positive non-container cargo"
    )
  )
}

analysis_output_path <- function(project_dir, stem, config, extension = "csv") {
  file.path(project_dir, ".traffic-work", paste0(stem, config$suffix, ".", extension))
}

analysis_treatment_config <- function() {
  # Quarter indices use 2014 Q1 = 0. Piraeus treatment begins in 2016 Q4;
  # the immediately preceding quarter, 2016 Q3, is the event-study base.
  treatment_index <- 11L
  list(
    treatment_quarter = "2016 Q4",
    base_quarter = "2016 Q3",
    treatment_index = treatment_index,
    pre_sample_event_times = -treatment_index:(23L - treatment_index),
    full_sample_event_times = -treatment_index:(47L - treatment_index),
    estimated_pre_event_times = -treatment_index:-2L,
    all_pre_event_times = -treatment_index:-1L,
    pre_sample_post_event_times = 0L:(23L - treatment_index),
    full_sample_post_event_times = 0L:(47L - treatment_index)
  )
}

set_analysis_outcome <- function(panel, config) {
  total_columns <- c(
    "External_Unloaded_Total", "External_Loaded_Total",
    "Internal_Unloaded_Total", "Internal_Loaded_Total"
  )
  container_columns <- c(
    "External_Unloaded_Containers", "External_Loaded_Containers",
    "Internal_Unloaded_Containers", "Internal_Loaded_Containers"
  )
  required <- unique(c(total_columns, container_columns))
  missing <- setdiff(required, names(panel))
  if (length(missing)) stop("Outcome columns missing from input: ", paste(missing, collapse = ", "))
  panel[required] <- lapply(panel[required], function(value) {
    value[is.na(value)] <- 0
    as.numeric(value)
  })
  total <- rowSums(panel[total_columns])
  containers <- rowSums(panel[container_columns])
  panel$total_cargo <- switch(
    config$key,
    "total" = total,
    "containers" = containers,
    "non-containers" = total - containers
  )
  if (any(panel$total_cargo < 0, na.rm = TRUE)) stop("Constructed outcome contains negative cargo")
  panel
}

analysis_donor_pools <- function(config) {
  if (config$key == "containers") {
    return(list(
      "Core container" = c("Volos", "Thessaloniki"),
      "Core + Heraklio" = c("Volos", "Thessaloniki", "Heraklio"),
      "All container-active" = c("Volos", "Thessaloniki", "Heraklio", "Lavrio")
    ))
  }
  list(
    "Narrow" = c("Volos", "Patra", "Igoumenitsa"),
    "Narrow + Thessaloniki" = c("Volos", "Patra", "Igoumenitsa", "Thessaloniki"),
    "Wide commercial" = c(
      "Volos", "Patra", "Igoumenitsa", "Thessaloniki",
      "Alexandroupolis", "Heraklio", "Kavala", "Lavrio",
      "Corfu", "Rodos", "Rafina", "Chios", "Mytiline", "Souda Bay"
    )
  )
}

event_term <- function(k) paste0("event_time::", k, ":treated")

linear_combination_test <- function(model, terms, weights, cluster_count) {
  beta <- coef(model)[terms]
  covariance <- vcov(model)[terms, terms, drop = FALSE]
  estimate <- sum(weights * beta)
  standard_error <- sqrt(as.numeric(t(weights) %*% covariance %*% weights))
  statistic <- estimate / standard_error
  p_value <- 2 * pt(-abs(statistic), df = cluster_count - 1)
  c(
    estimate = estimate,
    std_error = standard_error,
    statistic = statistic,
    p_value = p_value
  )
}

joint_wald_test <- function(model, terms, cluster_count) {
  beta <- coef(model)[terms]
  covariance <- vcov(model)[terms, terms, drop = FALSE]
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  tolerance <- max(decomposition$values) * 1e-8
  keep <- decomposition$values > tolerance
  rank <- sum(keep)
  if (rank == 0) stop("The restriction covariance has zero numerical rank")
  projected <- crossprod(decomposition$vectors[, keep, drop = FALSE], beta)
  chi_square <- sum((projected^2) / decomposition$values[keep])
  f_statistic <- chi_square / rank
  p_value <- pf(f_statistic, df1 = rank, df2 = cluster_count - 1, lower.tail = FALSE)
  c(
    restrictions = length(terms),
    covariance_rank = rank,
    f_statistic = f_statistic,
    df1 = rank,
    df2 = cluster_count - 1,
    p_value = p_value
  )
}

event_study_tests <- function(
  model,
  cluster_count,
  pre_event_times = analysis_treatment_config()$estimated_pre_event_times,
  post_event_times = analysis_treatment_config()$pre_sample_post_event_times
) {
  post_terms <- vapply(post_event_times, event_term, character(1))
  pre_terms <- vapply(pre_event_times, event_term, character(1))
  average_post <- linear_combination_test(
    model,
    post_terms,
    rep(1 / length(post_terms), length(post_terms)),
    cluster_count
  )
  joint_post <- joint_wald_test(model, post_terms, cluster_count)
  joint_pre <- joint_wald_test(model, pre_terms, cluster_count)
  list(average_post = average_post, joint_post = joint_post, joint_pre = joint_pre)
}

# Re-express a saturated seasonal event study relative to the average of all
# untreated observations from the same quarter of the year. The regression
# still needs one computationally omitted observation per season, but these
# arbitrary bases cancel from the reported linear contrasts. Consequently,
# every observed quarter receives an estimate and a covariance-based interval.
seasonal_average_event_study <- function(
  model,
  event_times,
  quarter_of_year,
  pre_event_times,
  cluster_count
) {
  coefficient_names <- names(coef(model))
  event_terms <- vapply(event_times, event_term, character(1))
  selector <- matrix(
    0,
    nrow = length(event_times),
    ncol = length(coefficient_names),
    dimnames = list(as.character(event_times), coefficient_names)
  )
  for (index in seq_along(event_terms)) {
    if (event_terms[index] %in% coefficient_names) {
      selector[index, event_terms[index]] <- 1
    }
  }

  contrast <- selector
  for (index in seq_along(event_times)) {
    seasonal_pre <- which(
      event_times %in% pre_event_times &
        quarter_of_year == quarter_of_year[index]
    )
    if (length(seasonal_pre) == 0) {
      stop("Every quarter of the year must have at least one pre-treatment observation")
    }
    contrast[index,] <- selector[index,] - colMeans(selector[seasonal_pre,, drop = FALSE])
  }

  beta <- coef(model)[coefficient_names]
  covariance <- vcov(model)[coefficient_names, coefficient_names, drop = FALSE]
  estimates <- as.numeric(contrast %*% beta)
  transformed_covariance <- contrast %*% covariance %*% t(contrast)
  standard_errors <- sqrt(pmax(diag(transformed_covariance), 0))
  statistics <- estimates / standard_errors
  p_values <- 2 * pt(-abs(statistics), df = cluster_count - 1)

  linear_test <- function(indexes, weights) {
    estimate <- sum(weights * estimates[indexes])
    variance <- as.numeric(
      t(weights) %*% transformed_covariance[indexes, indexes, drop = FALSE] %*% weights
    )
    standard_error <- sqrt(max(variance, 0))
    statistic <- estimate / standard_error
    c(
      estimate = estimate,
      std_error = standard_error,
      statistic = statistic,
      p_value = 2 * pt(-abs(statistic), df = cluster_count - 1)
    )
  }

  joint_test <- function(indexes) {
    restriction_estimates <- estimates[indexes]
    restriction_covariance <- transformed_covariance[indexes, indexes, drop = FALSE]
    decomposition <- eigen(
      (restriction_covariance + t(restriction_covariance)) / 2,
      symmetric = TRUE
    )
    tolerance <- max(decomposition$values) * 1e-8
    keep <- decomposition$values > tolerance
    rank <- sum(keep)
    if (rank == 0) stop("The seasonal-average restriction covariance has zero numerical rank")
    projected <- crossprod(decomposition$vectors[, keep, drop = FALSE], restriction_estimates)
    chi_square <- sum((projected^2) / decomposition$values[keep])
    f_statistic <- chi_square / rank
    c(
      restrictions = length(indexes),
      covariance_rank = rank,
      f_statistic = f_statistic,
      df1 = rank,
      df2 = cluster_count - 1,
      p_value = pf(f_statistic, df1 = rank, df2 = cluster_count - 1, lower.tail = FALSE)
    )
  }

  pre_indexes <- which(event_times %in% pre_event_times)
  post_indexes <- which(event_times >= 0)
  list(
    estimates = estimates,
    standard_errors = standard_errors,
    ci_low = estimates - qnorm(0.975) * standard_errors,
    ci_high = estimates + qnorm(0.975) * standard_errors,
    p_values = p_values,
    covariance = transformed_covariance,
    average_post = linear_test(post_indexes, rep(1 / length(post_indexes), length(post_indexes))),
    joint_post = joint_test(post_indexes),
    joint_pre = joint_test(pre_indexes)
  )
}

two_way_demean <- function(matrix_value) {
  sweep(sweep(matrix_value, 1, rowMeans(matrix_value)), 2, colMeans(matrix_value)) +
    mean(matrix_value)
}

# Restricted Rademacher wild-cluster bootstrap-t for the single DiD coefficient.
# With 15 or fewer port clusters every sign pattern is enumerated. Larger models
# use a fixed-seed Monte Carlo draw so the results are exactly reproducible.
wild_cluster_did <- function(sample, replications = 9999, seed = 20260820) {
  sample <- sample[order(sample$Port, sample$quarter_index),]
  ports <- sort(unique(sample$Port))
  quarters <- sort(unique(sample$quarter_index))
  cluster_count <- length(ports)
  time_count <- length(quarters)
  if (nrow(sample) != cluster_count * time_count) stop("Wild bootstrap requires a balanced panel")

  outcome <- matrix(sample$log_total, nrow = cluster_count, byrow = TRUE)
  treatment_post <- matrix(sample$treated * sample$post, nrow = cluster_count, byrow = TRUE)
  outcome_tilde <- two_way_demean(outcome)
  treatment_tilde <- two_way_demean(treatment_post)
  denominator <- sum(treatment_tilde^2)
  estimate <- sum(treatment_tilde * outcome_tilde) / denominator
  residual <- outcome_tilde - treatment_tilde * estimate

  observation_count <- nrow(sample)
  parameter_count <- cluster_count + time_count
  cr1 <- (cluster_count / (cluster_count - 1)) *
    ((observation_count - 1) / (observation_count - parameter_count))
  cluster_score <- rowSums(treatment_tilde * residual)
  standard_error <- sqrt(cr1 * sum(cluster_score^2) / denominator^2)
  observed_t <- estimate / standard_error

  restricted_residual <- outcome_tilde
  if (cluster_count <= 15) {
    weight_matrix <- as.matrix(expand.grid(rep(list(c(-1, 1)), cluster_count)))
    method <- paste0("Exact Rademacher wild-cluster bootstrap-t (", nrow(weight_matrix), " assignments)")
    exact <- TRUE
  } else {
    set.seed(seed)
    weight_matrix <- matrix(
      sample(c(-1, 1), replications * cluster_count, replace = TRUE),
      nrow = replications,
      ncol = cluster_count
    )
    method <- paste0("Rademacher wild-cluster bootstrap-t (", replications, " draws; seed ", seed, ")")
    exact <- FALSE
  }

  bootstrap_t <- numeric(nrow(weight_matrix))
  for (draw in seq_len(nrow(weight_matrix))) {
    bootstrap_outcome <- restricted_residual * weight_matrix[draw,]
    bootstrap_tilde <- two_way_demean(bootstrap_outcome)
    bootstrap_estimate <- sum(treatment_tilde * bootstrap_tilde) / denominator
    bootstrap_error <- bootstrap_tilde - treatment_tilde * bootstrap_estimate
    bootstrap_score <- rowSums(treatment_tilde * bootstrap_error)
    bootstrap_se <- sqrt(cr1 * sum(bootstrap_score^2) / denominator^2)
    bootstrap_t[draw] <- if (bootstrap_se > 0) bootstrap_estimate / bootstrap_se else 0
  }

  exceedances <- sum(abs(bootstrap_t) >= abs(observed_t) - sqrt(.Machine$double.eps))
  p_value <- if (exact) {
    exceedances / length(bootstrap_t)
  } else {
    (exceedances + 1) / (length(bootstrap_t) + 1)
  }
  c(
    estimate = estimate,
    std_error = standard_error,
    statistic = observed_t,
    p_value = p_value,
    replications = length(bootstrap_t),
    exact = as.integer(exact),
    method = method
  )
}
