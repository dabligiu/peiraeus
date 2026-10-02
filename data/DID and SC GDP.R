rm(list = ls())

pkgs <- c("readxl","readr","dplyr","tidyr","stringr","fixest","broom","ggplot2",
          "knitr","kableExtra","rmarkdown")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(fixest)
library(broom)
library(ggplot2)
library(knitr)
library(kableExtra)
library(rmarkdown)

# --- Paths ---
# Put the Excel file in the same folder as this script, or change this path.
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir  <- dirname(script_path)

data_dir   <- normalizePath(file.path(script_dir))
csv_path <- file.path(data_dir, "gdp pc.csv")
gdp_pc <- readr::read_csv2(csv_path)


# If your first column is not already called "region", rename it
names(gdp_pc)[1] <- "region"

# --- Convert from wide to long ---
gdp_long <- gdp_pc %>%
  pivot_longer(
    cols = -region,
    names_to = "year",
    values_to = "gdp_pc"
  ) %>%
  mutate(
    year = as.integer(year),
    gdp_pc = as.numeric(gdp_pc)
  )

# --- Define treatment ---
treatment_region <- "Peiraias, Nisoi"
treatment_year <- 2014

gdp_long <- gdp_long %>%
  mutate(
    treated = if_else(region == treatment_region, 1, 0),
    post = if_else(year >= treatment_year, 1, 0),
    did = treated * post,
    event_time = year - treatment_year
  )

# Check
glimpse(gdp_long)
table(gdp_long$treated, gdp_long$post)

did_model <- feols(
  gdp_pc ~ did | region + year,
  data = gdp_long,
  cluster = ~region
)

summary(did_model)
event_model <- feols(
  gdp_pc ~ i(event_time, treated, ref = -1) | region + year,
  data = gdp_long,
  cluster = ~region
)

summary(event_model)

iplot(
  event_model,
  ref.line = 0,
  xlab = "Years relative to treatment",
  ylab = "Effect on GDP per capita",
  main = "Event-study: Effect of COSCO investment on Peiraias, Nisoi"
)

gdp_plot <- gdp_long %>%
  mutate(group = if_else(region == treatment_region, "Peiraias, Nisoi", "Other regions")) %>%
  group_by(group, year) %>%
  summarise(
    mean_gdp_pc = mean(gdp_pc, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(gdp_plot, aes(x = year, y = mean_gdp_pc, linetype = group)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = treatment_year, linetype = "dashed") +
  labs(
    title = "GDP per capita before and after COSCO investment",
    x = "Year",
    y = "GDP per capita",
    linetype = ""
  ) +
  theme_minimal()


