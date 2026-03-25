# ==============================================================================
# 01b_filter_hopkins_outputs_by_latitude.R
# Filter observation- and cell-year outputs to a specified latitude band
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

setwd("C:/Users/robgu/Downloads")

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

IN_DIR  <- "hopkins_outputs_v11_gatePP"
OUT_DIR <- "hopkins_outputs_v11_gatePP_21to68N"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBS_IN   <- file.path(IN_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")
SYNC_IN  <- file.path(IN_DIR, "sync_cellyear_25km_with_covariates_plus_gatePP.csv")
SPEC_IN  <- file.path(IN_DIR, "hopkins_analysis_summary.csv")

OBS_OUT  <- file.path(OUT_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")
SYNC_OUT <- file.path(OUT_DIR, "sync_cellyear_25km_with_covariates_plus_gatePP.csv")
SPEC_OUT <- file.path(OUT_DIR, "hopkins_analysis_summary.csv")

LAT_MIN <- 21
LAT_MAX <- 68

need <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

need(OBS_IN)
need(SYNC_IN)
need(SPEC_IN)

# ------------------------------------------------------------------------------
# Read
# ------------------------------------------------------------------------------

obs <- read_csv(OBS_IN, show_col_types = FALSE)
sync <- read_csv(SYNC_IN, show_col_types = FALSE)
spec <- read_csv(SPEC_IN, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# Filter
# ------------------------------------------------------------------------------

obs_filt <- obs |>
  filter(
    is.finite(latitude),
    latitude >= LAT_MIN,
    latitude <= LAT_MAX
  )

sync_filt <- sync |>
  filter(
    is.finite(lat_median),
    lat_median >= LAT_MIN,
    lat_median <= LAT_MAX
  )

# keep only species still represented after filtering
species_keep <- union(
  unique(obs_filt$scientificName),
  unique(sync_filt$scientificName)
)

spec_filt <- spec |>
  filter(scientificName %in% species_keep)

# ------------------------------------------------------------------------------
# Write
# ------------------------------------------------------------------------------

write_csv(obs_filt, OBS_OUT)
write_csv(sync_filt, SYNC_OUT)
write_csv(spec_filt, SPEC_OUT)

cat("\n=== FILTER COMPLETE ===\n")
cat("Latitude band: ", LAT_MIN, " to ", LAT_MAX, " N\n", sep = "")
cat("Obs in / out:  ", nrow(obs), " / ", nrow(obs_filt), "\n", sep = "")
cat("Sync in / out: ", nrow(sync), " / ", nrow(sync_filt), "\n", sep = "")
cat("Species in / out: ", nrow(spec), " / ", nrow(spec_filt), "\n", sep = "")
cat("Wrote:\n")
cat("  ", OBS_OUT, "\n", sep = "")
cat("  ", SYNC_OUT, "\n", sep = "")
cat("  ", SPEC_OUT, "\n", sep = "")