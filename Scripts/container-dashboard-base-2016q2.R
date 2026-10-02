# Container dashboard variant 1:
# retain 2016 Q3 in the sample and use 2016 Q2 as the conventional omitted base.

if (identical(Sys.getenv("CONTAINER_SHOW_WARNINGS", unset = "false"), "true")) options(warn = 1)

variant_config <- list(
  id = "containers-base-2016q2",
  title = "Piraeus container cargo — 2016 Q2 omitted-base specification",
  output_html = "piraeus-port-traffic-containers-base-2016q2.html",
  normalization = "q2",
  drop_2016_q3 = FALSE
)

wrapper_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
wrapper_dir <- dirname(normalizePath(sub("^--file=", "", wrapper_arg)))
source(file.path(wrapper_dir, "container-variant-core.R"))
