# Container dashboard variant 2:
# exclude the partly treated 2016 Q3 observation and normalize event-study
# coefficients to the average of all retained, fully untreated pre-periods.

if (identical(Sys.getenv("CONTAINER_SHOW_WARNINGS", unset = "false"), "true")) options(warn = 1)

variant_config <- list(
  id = "containers-average-pre",
  title = "Piraeus container cargo — fully untreated pre-period normalization",
  output_html = "piraeus-port-traffic-containers-average-pre.html",
  normalization = "average_pre",
  drop_2016_q3 = TRUE
)

wrapper_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
wrapper_dir <- dirname(normalizePath(sub("^--file=", "", wrapper_arg)))
source(file.path(wrapper_dir, "container-variant-core.R"))
