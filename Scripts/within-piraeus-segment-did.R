# Within-Piraeus segment comparison: container versus non-container tonnage.
#
# Treatment begins in 2016 Q4. The potentially partially treated 2016 Q3 is
# excluded. The static estimate is a DiD with segment and quarter fixed effects;
# algebraically, it is the pre/post change in the log container/non-container
# tonnage ratio. Inference is Newey-West HAC on that quarterly log ratio because
# there are only two aggregate segment series and therefore no useful segment
# clustering dimension.

suppressPackageStartupMessages({
  library(readxl)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
project_dir <- dirname(dirname(script_path))
data_path <- file.path(project_dir, "data", "port-cargo-master-2014-2025.xlsx")
template_path <- file.path(project_dir, "Scripts", "within-piraeus-segment-template.html")
output_path <- file.path(project_dir, "piraeus-within-segment-comparison.html")
json_path <- file.path(project_dir, ".traffic-work", "within-piraeus-segment-results.json")
csv_path <- file.path(project_dir, ".traffic-work", "within-piraeus-segment-estimates.csv")

treatment_index <- 11L       # 2016 Q4, with 2014 Q1 = 0
excluded_index <- 10L        # 2016 Q3
q2_index <- 9L               # 2016 Q2
hac_lag <- 4L

raw <- as.data.frame(read_excel(data_path, sheet = "Port Quarterly Data"))
raw <- raw[raw$Port == "Piraeus",]
raw$quarter_index <- (raw$Year - 2014L) * 4L + raw$Quarter_Number - 1L
raw <- raw[order(raw$quarter_index),]

total_columns <- c(
  "External_Unloaded_Total", "External_Loaded_Total",
  "Internal_Unloaded_Total", "Internal_Loaded_Total"
)
container_columns <- c(
  "External_Unloaded_Containers", "External_Loaded_Containers",
  "Internal_Unloaded_Containers", "Internal_Loaded_Containers"
)
required <- c(total_columns, container_columns)
if (length(setdiff(required, names(raw)))) stop("Required cargo columns are missing")
raw[required] <- lapply(raw[required], function(x) { x[is.na(x)] <- 0; as.numeric(x) })
raw$container <- rowSums(raw[container_columns])
raw$total <- rowSums(raw[total_columns])
raw$noncontainer <- raw$total - raw$container
if (any(raw$container <= 0) || any(raw$noncontainer <= 0)) {
  stop("Log specification requires positive container and non-container tonnage")
}
raw$log_ratio <- log(raw$container) - log(raw$noncontainer)
raw$post <- as.integer(raw$quarter_index >= treatment_index)
raw$excluded <- raw$quarter_index == excluded_index

# Newey-West covariance with Bartlett weights. This operates on the quarterly
# log ratio, which is the collapsed form of the two-segment FE DiD.
hac_fit <- function(y, X, lag = 4L) {
  X <- as.matrix(X)
  fit <- lm.fit(X, y)
  u <- as.numeric(fit$residuals)
  n <- nrow(X)
  k <- ncol(X)
  xu <- X * u
  meat <- crossprod(xu)
  if (lag > 0L) {
    for (ell in seq_len(min(lag, n - 1L))) {
      weight <- 1 - ell / (lag + 1)
      gamma <- crossprod(xu[(ell + 1L):n,, drop = FALSE], xu[1L:(n - ell),, drop = FALSE])
      meat <- meat + weight * (gamma + t(gamma))
    }
  }
  bread <- solve(crossprod(X))
  vcov <- (n / (n - k)) * bread %*% meat %*% bread
  list(coefficients = setNames(as.numeric(fit$coefficients), colnames(X)), vcov = vcov,
       residuals = u, df = n - k, n = n, lag = min(lag, n - 1L))
}

coefficient_result <- function(fit, term) {
  position <- match(term, names(fit$coefficients))
  estimate <- fit$coefficients[position]
  se <- sqrt(fit$vcov[position, position])
  critical <- qt(.975, fit$df)
  statistic <- estimate / se
  list(
    estimate = estimate,
    percent = 100 * (exp(estimate) - 1),
    standardError = se,
    ciLow = estimate - critical * se,
    ciHigh = estimate + critical * se,
    percentCiLow = 100 * (exp(estimate - critical * se) - 1),
    percentCiHigh = 100 * (exp(estimate + critical * se) - 1),
    pValue = 2 * pt(-abs(statistic), fit$df)
  )
}

estimate_window <- function(end_year) {
  d <- raw[raw$Year <= end_year & !raw$excluded,]
  X <- cbind(`(Intercept)` = 1, post = d$post)
  static_fit <- hac_fit(d$log_ratio, X, hac_lag)
  static <- coefficient_result(static_fit, "post")

  # Differential seasonality robustness: allow the container/non-container
  # ratio to have its own recurring quarter-of-year pattern.
  season <- model.matrix(~ post + factor(Quarter_Number), data = d)
  seasonal_fit <- hac_fit(d$log_ratio, season, hac_lag)
  seasonal <- coefficient_result(seasonal_fit, "post")

  pre <- d[d$post == 0L,]
  post <- d[d$post == 1L,]
  q2_ratio <- raw$log_ratio[raw$quarter_index == q2_index]
  average_pre_ratio <- mean(pre$log_ratio)
  q2_post_log <- mean(post$log_ratio - q2_ratio)
  average_pre_post_log <- mean(post$log_ratio - average_pre_ratio)

  list(
    label = if (end_year == 2019L) "Short sample: 2014 Q1–2019 Q4" else "Full sample: 2014 Q1–2025 Q4",
    endYear = end_year,
    quarters = nrow(d),
    preQuarters = nrow(pre),
    postQuarters = nrow(post),
    static = static,
    seasonal = seasonal,
    containerGrowth = 100 * (exp(mean(log(post$container)) - mean(log(pre$container))) - 1),
    noncontainerGrowth = 100 * (exp(mean(log(post$noncontainer)) - mean(log(pre$noncontainer))) - 1),
    q2AveragePostLog = q2_post_log,
    q2AveragePostPercent = 100 * (exp(q2_post_log) - 1),
    averagePrePostLog = average_pre_post_log,
    averagePrePostPercent = 100 * (exp(average_pre_post_log) - 1)
  )
}

retained <- raw[!raw$excluded,]
pre <- retained[retained$post == 0L,]
pre_time <- pre$quarter_index - min(pre$quarter_index)
trend_fit <- hac_fit(pre$log_ratio, cbind(`(Intercept)` = 1, time = pre_time), hac_lag)
pretrend <- coefficient_result(trend_fit, "time")

# Verify numerically that the collapsed ratio estimate equals the coefficient
# from the explicit stacked segment-by-quarter FE regression.
verify_stacked <- function(end_year) {
  d <- raw[raw$Year <= end_year & !raw$excluded,]
  stacked <- rbind(
    data.frame(Quarter = d$Quarter, segment = "Container", log_tonnage = log(d$container), post = d$post),
    data.frame(Quarter = d$Quarter, segment = "Non-container", log_tonnage = log(d$noncontainer), post = d$post)
  )
  stacked$container_segment <- as.integer(stacked$segment == "Container")
  model <- lm(log_tonnage ~ factor(segment) + factor(Quarter) + container_segment:post, data = stacked)
  unname(coef(model)["container_segment:post"])
}

q2_ratio <- raw$log_ratio[raw$quarter_index == q2_index]
average_pre_ratio <- mean(pre$log_ratio)
series <- lapply(seq_len(nrow(raw)), function(i) list(
  quarter = raw$Quarter[i],
  year = raw$Year[i],
  quarterNumber = raw$Quarter_Number[i],
  eventTime = raw$quarter_index[i] - treatment_index,
  container = raw$container[i],
  noncontainer = raw$noncontainer[i],
  excluded = raw$excluded[i],
  averagePreEffect = 100 * (exp(raw$log_ratio[i] - average_pre_ratio) - 1),
  q2Effect = 100 * (exp(raw$log_ratio[i] - q2_ratio) - 1)
))

windows <- list(estimate_window(2019L), estimate_window(2025L))
checks <- lapply(c(2019L, 2025L), function(end_year) {
  explicit <- verify_stacked(end_year)
  collapsed <- windows[[if (end_year == 2019L) 1L else 2L]]$static$estimate
  list(endYear = end_year, stackedCoefficient = explicit, collapsedCoefficient = collapsed,
       absoluteDifference = abs(explicit - collapsed))
})
if (any(vapply(checks, `[[`, numeric(1), "absoluteDifference") > 1e-10)) {
  stop("Stacked FE and collapsed ratio estimates do not match")
}

payload <- list(
  meta = list(
    title = "Within-Piraeus segment comparison",
    treatmentQuarter = "2016 Q4",
    excludedQuarter = "2016 Q3",
    primaryNormalization = "Average of fully untreated pre-treatment quarters",
    conventionalNormalization = "2016 Q2",
    outcome = "Log quarterly tonnage",
    inference = paste0("Newey-West HAC, Bartlett kernel, lag ", hac_lag),
    source = "data/port-cargo-master-2014-2025.xlsx"
  ),
  windows = windows,
  pretrend = pretrend,
  series = series,
  equivalenceChecks = checks
)

dir.create(dirname(json_path), showWarnings = FALSE, recursive = TRUE)
write_json(payload, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 15)
estimate_rows <- do.call(rbind, lapply(windows, function(x) data.frame(
  sample = x$label,
  quarters = x$quarters,
  container_growth_percent = x$containerGrowth,
  noncontainer_growth_percent = x$noncontainerGrowth,
  did_log_points = x$static$estimate,
  did_percent = x$static$percent,
  hac_standard_error = x$static$standardError,
  ci_low_percent = x$static$percentCiLow,
  ci_high_percent = x$static$percentCiHigh,
  p_value = x$static$pValue,
  seasonal_did_percent = x$seasonal$percent,
  seasonal_p_value = x$seasonal$pValue,
  q2_normalized_average_post_percent = x$q2AveragePostPercent,
  average_pre_normalized_average_post_percent = x$averagePrePostPercent
)))
write.csv(estimate_rows, csv_path, row.names = FALSE)

template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
json_text <- toJSON(payload, auto_unbox = TRUE, na = "null", digits = 15)
html <- sub("__WITHIN_PIRAEUS_DATA__", json_text, template, fixed = TRUE)
writeLines(html, output_path, useBytes = TRUE)

cat("Created:", output_path, "\n")
cat("Results:", json_path, "\n")
for (x in windows) {
  cat(sprintf("%s: DiD %+.1f%% (95%% CI %+.1f%% to %+.1f%%; p = %.3f); Q2-normalized average post %+.1f%%\n",
              x$label, x$static$percent, x$static$percentCiLow, x$static$percentCiHigh,
              x$static$pValue, x$q2AveragePostPercent))
}
cat(sprintf("Pre-period differential linear trend: %+.2f%% per quarter; p = %.3f\n",
            100 * (exp(pretrend$estimate) - 1), pretrend$pValue))
