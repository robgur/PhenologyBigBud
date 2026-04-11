# ==============================================================================
# 01b_filter_hopkins_outputs_by_latitude.R
# Filter observation- and cell-year outputs to a specified latitude band
# and remove rows with NA in tmeanC_Jan
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

setwd("C:/Users/robgu/Downloads")

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

IN_DIR  <- "hopkins_outputs_v11_gatePP_obsOnly"
OUT_DIR <- "hopkins_outputs_v11_gatePP_obsOnly21to68N"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBS_IN  <- file.path(IN_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")
OBS_OUT  <- file.path(OUT_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")


LAT_MIN <- 21
LAT_MAX <- 68

need <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

need(OBS_IN)


# ------------------------------------------------------------------------------
# Read
# ------------------------------------------------------------------------------

obs  <- read_csv(OBS_IN, show_col_types = FALSE)


# ------------------------------------------------------------------------------
# Filter
# ------------------------------------------------------------------------------

obs_filt <- obs |>
  filter(
    is.finite(latitude),
    latitude >= LAT_MIN,
    latitude <= LAT_MAX
  ) |>
  drop_na(tmeanC_Jan,hfp)


# ------------------------------------------------------------------------------
# Write
# ------------------------------------------------------------------------------

write_csv(obs_filt, OBS_OUT)
cat("\n=== FILTER COMPLETE ===\n")
cat("Latitude band: ", LAT_MIN, " to ", LAT_MAX, " N\n", sep = "")
cat("Obs in / out:  ", nrow(obs), " / ", nrow(obs_filt), "\n", sep = "")
cat("Wrote:\n")
cat("  ", OBS_OUT, "\n", sep = "")
