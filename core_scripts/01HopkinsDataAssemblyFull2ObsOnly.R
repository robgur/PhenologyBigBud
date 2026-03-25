# ==============================================================================
# 01_build_hopkins_obs_gatePP_inputs.R
# Observation-level Hopkins + gatePP inputs
#
# Produces:
#   1) observation-level output with Hopkins anomaly, ERA covariates,
#      photoperiod covariates, March 15 sensitivity covariates, and cell_id_25km
#   2) species summary table
#
# Does NOT build synchrony output.
#
# FIX INCLUDED:
#   - avoids many-to-many join warning by assigning a stable .obs_id to each
#     record before ERA extraction and joining back by .obs_id only
# ==============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(terra)
  library(readr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 0) USER SETTINGS
# ------------------------------------------------------------------------------

setwd("C:/Users/robgu/Downloads")

# Inputs
INFILE   <- "phenovision_leaves_03_15_2026.csv"
DEM_PATH <- "dtm.bareearth_ensemble_p10_250m_s_2018_go_epsg4326_v20230221.tif"
BIO4_DETRENDED_PATH <- "bio4_lat_elev_detrended_NH.tif"

ERA5_PATH <- "era5land_t2m_2017_2025_combined.nc"
ERA_VAR   <- "t2m"
ORIGIN_DATE <- as.Date("2017-01-01")
ERA_COVERAGE_YEARS <- 2017:2025

# Outputs
OUT_DIR <- "hopkins_outputs_v11_gatePP_obsOnly"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_OBS     <- file.path(OUT_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")
OUT_SPECIES <- file.path(OUT_DIR, "hopkins_analysis_summary.csv")

# Filters
TRAIT_FILTER <- "breaking buds"
MIN_YEAR     <- 2017
NH_ONLY      <- TRUE
DOY_MAX      <- 215

# Species filtering
MIN_N <- 75
MIN_LAT_RANGE_DEG <- 3.0
RANGE_Q_LO <- 0.05
RANGE_Q_HI <- 0.95

# Anchor fitting thresholds
ANCHOR_MIN_N <- 10
ANCHOR_MIN_UNIQUE_LAT <- 3
ANCHOR_MIN_FINITE_ELEV <- 6

# Hopkins constants
HOPKINS_LAT_SLOPE  <- 4
HOPKINS_ELEV_SLOPE <- 4/122

# 25 km cell IDs retained for downstream random effects
SYNC_CELL_KM     <- 25
SYNC_CELL_SIZE_M <- SYNC_CELL_KM * 1000
SYNC_PROJ_EPSG   <- "EPSG:6933"

# Temperature forcing windows
GDD_BASE_C <- 5

# Main windows
EARLY_END_MONTH <- 2
EARLY_END_DAY   <- 28

LATE_START_MONTH <- 3
LATE_START_DAY   <- 1
LATE_END_MONTH   <- 4
LATE_END_DAY     <- 30

# March 15 sensitivity windows
EARLY_END_MONTH_M15 <- 3
EARLY_END_DAY_M15   <- 15

LATE_START_MONTH_M15 <- 3
LATE_START_DAY_M15   <- 16
LATE_END_MONTH_M15   <- 4
LATE_END_DAY_M15     <- 30

PP_GATE_LABEL <- "Mar01"

# ------------------------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------------------------

stopif_missing <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
}

collapse_to_binomial <- function(x) {
  x <- trimws(x)
  parts <- strsplit(x, "\\s+")
  vapply(parts, function(p) {
    if (length(p) >= 2) paste(p[1], p[2]) else p[1]
  }, character(1))
}

remove_outliers_iqr <- function(df, vars, k = 3) {
  keep <- rep(TRUE, nrow(df))
  for (v in vars) {
    x <- df[[v]]
    ok <- is.finite(x)
    if (sum(ok) < 5) next
    q1 <- quantile(x[ok], 0.25, na.rm = TRUE)
    q3 <- quantile(x[ok], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lo <- q1 - k * iqr
    hi <- q3 + k * iqr
    keep <- keep & (!ok | (x >= lo & x <= hi))
  }
  df[keep, , drop = FALSE]
}

robust_range <- function(x, q_lo = 0.05, q_hi = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(list(lo = NA_real_, hi = NA_real_, range = NA_real_))
  lo <- as.numeric(quantile(x, q_lo, na.rm = TRUE))
  hi <- as.numeric(quantile(x, q_hi, na.rm = TRUE))
  list(lo = lo, hi = hi, range = hi - lo)
}

add_elevation <- function(df, dem) {
  df <- df |> mutate(elevation = NA_real_)
  ok <- is.finite(df$longitude) & is.finite(df$latitude)
  if (!any(ok)) return(df)
  coords <- df[ok, c("longitude", "latitude")]
  elev <- terra::extract(dem, coords, method = "bilinear")[, 2]
  df$elevation[ok] <- elev
  df
}

add_raster_values <- function(df, r, out_col, method = "bilinear") {
  df <- df |> mutate(!!out_col := NA_real_)
  ok <- is.finite(df$longitude) & is.finite(df$latitude)
  if (!any(ok)) return(df)
  coords <- df[ok, c("longitude", "latitude")]
  vals <- terra::extract(r, coords, method = method)[, 2]
  df[[out_col]][ok] <- vals
  df
}

estimate_species_anchors <- function(df) {
  df |>
    group_by(scientificName) |>
    group_modify(~{
      x <- .x |>
        filter(is.finite(dayOfYear), is.finite(latitude))
      
      if (nrow(x) == 0) {
        return(tibble(
          start_lat = NA_real_,
          start_elev = NA_real_,
          start_doy = NA_real_,
          anchor_model = NA_character_,
          n_anchor = 0L,
          n_elev_anchor = 0L
        ))
      }
      
      slat <- median(x$latitude, na.rm = TRUE)
      selev <- if (sum(is.finite(x$elevation)) > 0) median(x$elevation, na.rm = TRUE) else NA_real_
      
      enough_lat <- nrow(x) >= ANCHOR_MIN_N &&
        length(unique(x$latitude[is.finite(x$latitude)])) >= ANCHOR_MIN_UNIQUE_LAT
      
      enough_elev <- sum(is.finite(x$elevation)) >= ANCHOR_MIN_FINITE_ELEV
      
      sdoy <- NA_real_
      anchor_model <- NA_character_
      
      if (enough_lat && enough_elev) {
        x_fit <- x |> filter(is.finite(elevation))
        m <- tryCatch(
          lm(dayOfYear ~ latitude + elevation, data = x_fit),
          error = function(e) NULL
        )
        if (!is.null(m) && all(is.finite(coef(m)))) {
          sdoy <- as.numeric(
            predict(m, newdata = data.frame(latitude = slat, elevation = selev))
          )
          anchor_model <- "lat_elev"
        }
      }
      
      if (!is.finite(sdoy) && enough_lat) {
        m <- tryCatch(
          lm(dayOfYear ~ latitude, data = x),
          error = function(e) NULL
        )
        if (!is.null(m) && all(is.finite(coef(m)))) {
          sdoy <- as.numeric(
            predict(m, newdata = data.frame(latitude = slat))
          )
          anchor_model <- "lat_only"
        }
      }
      
      if (!is.finite(sdoy)) {
        sdoy <- median(x$dayOfYear, na.rm = TRUE)
        anchor_model <- "median_doy"
      }
      
      tibble(
        start_lat = slat,
        start_elev = selev,
        start_doy = sdoy,
        anchor_model = anchor_model,
        n_anchor = nrow(x),
        n_elev_anchor = sum(is.finite(x$elevation))
      )
    }) |>
    ungroup()
}

add_hopkins_predictions <- function(df) {
  df |>
    mutate(
      elev_dev = case_when(
        is.finite(elevation) & is.finite(start_elev) ~ elevation - start_elev,
        is.finite(elevation) & !is.finite(start_elev) ~ elevation,
        TRUE ~ 0
      ),
      hopkins_lat = start_doy + (latitude - start_lat) * HOPKINS_LAT_SLOPE,
      hopkins_lat_elev = start_doy +
        (latitude - start_lat) * HOPKINS_LAT_SLOPE +
        elev_dev * HOPKINS_ELEV_SLOPE,
      anom = dayOfYear - hopkins_lat_elev
    )
}

add_sync_grid_25km <- function(df, cell_size_m = 25000, proj_epsg = "EPSG:6933") {
  ok <- is.finite(df$longitude) & is.finite(df$latitude)
  df <- df |> mutate(cell_id_25km = NA_character_, x_m = NA_real_, y_m = NA_real_)
  if (!any(ok)) return(df)
  
  v <- terra::vect(df[ok, ], geom = c("longitude", "latitude"), crs = "EPSG:4326", keepgeom = TRUE)
  v_m <- terra::project(v, proj_epsg)
  xy <- terra::crds(v_m)
  
  cx <- floor(xy[, 1] / cell_size_m)
  cy <- floor(xy[, 2] / cell_size_m)
  cell_id <- paste0(cx, "_", cy)
  
  df$cell_id_25km[ok] <- cell_id
  df$x_m[ok] <- xy[, 1]
  df$y_m[ok] <- xy[, 2]
  df
}

to_0_360 <- function(lon) ifelse(is.na(lon), NA_real_, ifelse(lon < 0, lon + 360, lon))

load_era <- function(path, var = NULL) {
  stopif_missing(file.exists(path), paste0("Missing ERA5 file: ", path))
  tryCatch(
    if (!is.null(var)) terra::rast(path, subds = var) else terra::rast(path),
    error = function(e) terra::rast(path)
  )
}

era_all_dates <- function(era, origin = ORIGIN_DATE) {
  seq.Date(origin, by = "day", length.out = terra::nlyr(era))
}

extract_layers_pts <- function(era, idx, pts) {
  if (length(idx) == 0) {
    return(matrix(NA_real_, nrow = nrow(pts), ncol = 0))
  }
  vals <- terra::extract(era[[idx]], pts)
  as.matrix(vals[, -1, drop = FALSE])
}

idx_between <- function(all_dates, start_date, end_date) {
  which(all_dates >= start_date & all_dates <= end_date)
}

make_year_layer_sets <- function(year, all_dates,
                                 early_end_month, early_end_day,
                                 late_start_month, late_start_day,
                                 late_end_month, late_end_day) {
  
  y_start <- as.Date(sprintf("%04d-01-01", year))
  y_end   <- as.Date(sprintf("%04d-12-31", year))
  
  early_end  <- as.Date(sprintf("%04d-%02d-%02d", year, early_end_month, early_end_day))
  late_start <- as.Date(sprintf("%04d-%02d-%02d", year, late_start_month, late_start_day))
  late_end   <- as.Date(sprintf("%04d-%02d-%02d", year, late_end_month, late_end_day))
  
  early_end  <- max(min(early_end, y_end), y_start)
  late_start <- max(min(late_start, y_end), y_start)
  late_end   <- max(min(late_end, y_end), y_start)
  
  idx_jan <- idx_between(all_dates,
                         as.Date(sprintf("%04d-01-01", year)),
                         as.Date(sprintf("%04d-01-31", year)))
  idx_fma <- idx_between(all_dates,
                         as.Date(sprintf("%04d-02-01", year)),
                         as.Date(sprintf("%04d-04-30", year)))
  idx_ma  <- idx_between(all_dates,
                         as.Date(sprintf("%04d-03-01", year)),
                         as.Date(sprintf("%04d-04-30", year)))
  
  idx_gdd_early <- idx_between(all_dates, y_start, early_end)
  idx_gdd_late  <- idx_between(all_dates, late_start, late_end)
  
  list(
    idx_jan = idx_jan,
    idx_fma = idx_fma,
    idx_ma = idx_ma,
    idx_gdd_early = idx_gdd_early,
    idx_gdd_late = idx_gdd_late
  )
}

extract_era_covariates_for_year <- function(era, all_dates, year, pts,
                                            early_end_month, early_end_day,
                                            late_start_month, late_start_day,
                                            late_end_month, late_end_day,
                                            gdd_base_c = 5) {
  idx <- make_year_layer_sets(year, all_dates,
                              early_end_month, early_end_day,
                              late_start_month, late_start_day,
                              late_end_month, late_end_day)
  
  v_jan <- extract_layers_pts(era, idx$idx_jan, pts)
  v_fma <- extract_layers_pts(era, idx$idx_fma, pts)
  v_ma  <- extract_layers_pts(era, idx$idx_ma,  pts)
  v_e   <- extract_layers_pts(era, idx$idx_gdd_early, pts)
  v_l   <- extract_layers_pts(era, idx$idx_gdd_late,  pts)
  
  v_jan <- v_jan - 273.15
  v_fma <- v_fma - 273.15
  v_ma  <- v_ma  - 273.15
  v_e   <- v_e   - 273.15
  v_l   <- v_l   - 273.15
  
  tibble(
    tmeanC_Jan       = rowMeans(v_jan, na.rm = TRUE),
    tmeanC_FebMarApr = rowMeans(v_fma, na.rm = TRUE),
    tmeanC_MarApr    = rowMeans(v_ma,  na.rm = TRUE),
    gddC_early       = rowSums(pmax(v_e - gdd_base_c, 0), na.rm = TRUE),
    gddC_late        = rowSums(pmax(v_l - gdd_base_c, 0), na.rm = TRUE)
  ) |>
    mutate(gddC_total = gddC_early + gddC_late)
}

daylength_hours <- function(lat_deg, doy) {
  lat_rad <- lat_deg * pi / 180
  decl <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  
  x <- -tan(lat_rad) * tan(decl)
  x <- pmin(pmax(x, -1), 1)
  
  24 * acos(x) / pi
}

add_photoperiod_fixed <- function(df, lat_col = "latitude") {
  lat <- df[[lat_col]]
  
  df |>
    mutate(
      pp_Mar01      = daylength_hours(lat, 60),
      pp_Equinox    = daylength_hours(lat, 79),
      pp_Apr01      = daylength_hours(lat, 91),
      dpp_Mar01     = daylength_hours(lat, 61) - daylength_hours(lat, 59),
      dpp_Equinox   = daylength_hours(lat, 80) - daylength_hours(lat, 78),
      dpp_Apr01     = daylength_hours(lat, 92) - daylength_hours(lat, 90)
    )
}

add_gate_aliases <- function(df, gate_label = "Mar01") {
  gate_label <- match.arg(gate_label, c("Mar01", "Equinox", "Apr01"))
  
  if (gate_label == "Mar01") {
    df |> mutate(pp_gate = pp_Mar01, dpp_gate = dpp_Mar01)
  } else if (gate_label == "Equinox") {
    df |> mutate(pp_gate = pp_Equinox, dpp_gate = dpp_Equinox)
  } else {
    df |> mutate(pp_gate = pp_Apr01, dpp_gate = dpp_Apr01)
  }
}

# ------------------------------------------------------------------------------
# 2) LOAD + FILTER OBSERVATIONS
# ------------------------------------------------------------------------------

ds <- arrow::open_dataset(INFILE, format = "csv")

dat <- ds |>
  select(scientificName, year, dayOfYear, latitude, longitude, verbatimTrait) |>
  filter(verbatimTrait == TRAIT_FILTER, year >= MIN_YEAR) |>
  collect()

dat <- dat |>
  mutate(scientificName = collapse_to_binomial(scientificName))

if (NH_ONLY) dat <- dat |> filter(latitude >= 0)
dat <- dat |> filter(dayOfYear < DOY_MAX)

species_keep <- dat |>
  count(scientificName, name = "n") |>
  filter(n >= MIN_N) |>
  pull(scientificName)

dat <- dat |> filter(scientificName %in% species_keep)

# Stable observation ID for safe joins back after ERA extraction
dat <- dat |> mutate(.obs_id = row_number())

# ------------------------------------------------------------------------------
# 3) ELEVATION + BIO4 DETRENDED
# ------------------------------------------------------------------------------

if (file.exists(DEM_PATH)) {
  dem <- terra::rast(DEM_PATH)
  dat <- dat |> add_elevation(dem)
} else {
  warning("DEM file not found; proceeding without elevation.")
  dat$elevation <- NA_real_
}

if (file.exists(BIO4_DETRENDED_PATH)) {
  bio4_det <- terra::rast(BIO4_DETRENDED_PATH)
  dat <- dat |> add_raster_values(
    bio4_det,
    out_col = "bio4_lat_detrended",
    method = "bilinear"
  )
} else {
  warning("bio4_lat_elev_detrended_NH.tif not found; proceeding without bio4 detrended seasonality.")
  dat$bio4_lat_detrended <- NA_real_
}

# ------------------------------------------------------------------------------
# 4) SPECIES HOPKINS ANCHORS
# ------------------------------------------------------------------------------

species_params <- estimate_species_anchors(dat)
dat <- dat |> inner_join(species_params, by = "scientificName")

# ------------------------------------------------------------------------------
# 5) REMOVE NARROW-LAT-RANGE SPECIES
# ------------------------------------------------------------------------------

lat_ranges_clean <- dat |>
  group_by(scientificName) |>
  group_modify(~{
    x_clean <- remove_outliers_iqr(.x, vars = c("latitude", "dayOfYear"), k = 3)
    lat_q05 <- as.numeric(quantile(x_clean$latitude, RANGE_Q_LO, na.rm = TRUE))
    lat_q95 <- as.numeric(quantile(x_clean$latitude, RANGE_Q_HI, na.rm = TRUE))
    tibble(lat_range_q05_q95_clean = lat_q95 - lat_q05)
  }) |>
  ungroup()

species_keep2 <- lat_ranges_clean |>
  filter(is.finite(lat_range_q05_q95_clean), lat_range_q05_q95_clean >= MIN_LAT_RANGE_DEG) |>
  pull(scientificName)

dat <- dat |> filter(scientificName %in% species_keep2)

# ------------------------------------------------------------------------------
# 6) ADD HOPKINS + ANOM + CELL IDS + PHOTOPERIOD
# ------------------------------------------------------------------------------

dat <- dat |>
  add_hopkins_predictions() |>
  add_sync_grid_25km(
    cell_size_m = SYNC_CELL_SIZE_M,
    proj_epsg = SYNC_PROJ_EPSG
  ) |>
  add_photoperiod_fixed(lat_col = "latitude") |>
  add_gate_aliases(gate_label = PP_GATE_LABEL)

# ------------------------------------------------------------------------------
# 7) ERA5 TEMPERATURE COVARIATES AT OBSERVATION SUPPORT
# ------------------------------------------------------------------------------

era <- load_era(ERA5_PATH, var = ERA_VAR)
all_dates <- era_all_dates(era, ORIGIN_DATE)

obs_for_era <- dat |>
  filter(year %in% ERA_COVERAGE_YEARS, is.finite(latitude), is.finite(longitude)) |>
  mutate(
    lon_era = to_0_360(longitude),
    .row_id = row_number()
  )

pts_obs <- terra::vect(
  obs_for_era |> transmute(lon = lon_era, lat = latitude),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

idx_by_year <- split(seq_len(nrow(obs_for_era)), obs_for_era$year)

era_obs_cov_list <- lapply(names(idx_by_year), function(yy_chr) {
  yy <- as.integer(yy_chr)
  ii <- idx_by_year[[yy_chr]]
  
  cov <- extract_era_covariates_for_year(
    era, all_dates, year = yy, pts = pts_obs[ii],
    early_end_month = EARLY_END_MONTH, early_end_day = EARLY_END_DAY,
    late_start_month = LATE_START_MONTH, late_start_day = LATE_START_DAY,
    late_end_month = LATE_END_MONTH, late_end_day = LATE_END_DAY,
    gdd_base_c = GDD_BASE_C
  )
  
  tibble(.row_id = ii, year = yy) |>
    bind_cols(cov)
})

era_obs_cov <- bind_rows(era_obs_cov_list)

era_obs_cov_list_m15 <- lapply(names(idx_by_year), function(yy_chr) {
  yy <- as.integer(yy_chr)
  ii <- idx_by_year[[yy_chr]]
  
  cov <- extract_era_covariates_for_year(
    era, all_dates, year = yy, pts = pts_obs[ii],
    early_end_month = EARLY_END_MONTH_M15, early_end_day = EARLY_END_DAY_M15,
    late_start_month = LATE_START_MONTH_M15, late_start_day = LATE_START_DAY_M15,
    late_end_month = LATE_END_MONTH_M15, late_end_day = LATE_END_DAY_M15,
    gdd_base_c = GDD_BASE_C
  ) |>
    rename(
      gddC_early_m15 = gddC_early,
      gddC_late_m15  = gddC_late,
      gddC_total_m15 = gddC_total
    )
  
  tibble(.row_id = ii, year = yy) |>
    bind_cols(cov)
})

era_obs_cov_m15 <- bind_rows(era_obs_cov_list_m15)

obs_for_era <- obs_for_era |>
  left_join(era_obs_cov, by = c(".row_id", "year")) |>
  left_join(
    era_obs_cov_m15 |>
      select(.row_id, year, gddC_early_m15, gddC_late_m15, gddC_total_m15),
    by = c(".row_id", "year")
  ) |>
  select(
    .obs_id,
    tmeanC_Jan, tmeanC_FebMarApr, tmeanC_MarApr,
    gddC_early, gddC_late, gddC_total,
    gddC_early_m15, gddC_late_m15, gddC_total_m15
  )

dat_obs <- dat |>
  left_join(obs_for_era, by = ".obs_id") |>
  select(-.obs_id)

# ------------------------------------------------------------------------------
# 8) SPECIES SUMMARY
# ------------------------------------------------------------------------------

dat_for_global <- dat_obs |>
  group_by(scientificName) |>
  group_modify(~ remove_outliers_iqr(.x, vars = c("latitude", "dayOfYear"), k = 3)) |>
  ungroup()

global_lat_model <- lm(dayOfYear ~ latitude, data = dat_for_global)
GLOBAL_ALPHA <- unname(coef(global_lat_model)[["(Intercept)"]])
GLOBAL_BETA  <- unname(coef(global_lat_model)[["latitude"]])

fit_one_species_summary <- function(df) {
  df <- remove_outliers_iqr(df, vars = c("latitude", "dayOfYear"), k = 3)
  if (nrow(df) < 30) return(tibble())
  
  lat_median <- median(df$latitude, na.rm = TRUE)
  elev_median <- median(df$elevation, na.rm = TRUE)
  median_doy_clean <- median(df$dayOfYear, na.rm = TRUE)
  
  lat_rr  <- robust_range(df$latitude, RANGE_Q_LO, RANGE_Q_HI)
  elev_rr <- robust_range(df$elevation, RANGE_Q_LO, RANGE_Q_HI)
  bio4_rr <- robust_range(df$bio4_lat_detrended, RANGE_Q_LO, RANGE_Q_HI)
  
  expected_doy_at_median_lat <- GLOBAL_ALPHA + GLOBAL_BETA * lat_median
  lat_standardized_doy_index <- median_doy_clean - expected_doy_at_median_lat
  
  mod_lat <- lm(dayOfYear ~ latitude, data = df)
  lat_beta <- unname(coef(mod_lat)[["latitude"]])
  lat_r2 <- summary(mod_lat)$r.squared
  
  has_elev <- sum(is.finite(df$elevation)) > 10
  mod_lat_elev <- if (has_elev) lm(dayOfYear ~ latitude + elevation, data = df) else NULL
  elev_beta <- if (!is.null(mod_lat_elev) && "elevation" %in% names(coef(mod_lat_elev))) {
    unname(coef(mod_lat_elev)[["elevation"]])
  } else NA_real_
  lat_elev_r2 <- if (!is.null(mod_lat_elev)) summary(mod_lat_elev)$r.squared else NA_real_
  
  hop_res <- df$dayOfYear - df$hopkins_lat_elev
  hop_rmse <- sqrt(mean(hop_res^2, na.rm = TRUE))
  
  tibble(
    scientificName = df$scientificName[1],
    n_records = nrow(df),
    latitude_median = lat_median,
    elevation_median = elev_median,
    lat_q05 = lat_rr$lo,
    lat_q95 = lat_rr$hi,
    lat_range_q05_q95 = lat_rr$range,
    elev_q05 = elev_rr$lo,
    elev_q95 = elev_rr$hi,
    elev_range_q05_q95 = elev_rr$range,
    bio4_detrended_median = median(df$bio4_lat_detrended, na.rm = TRUE),
    bio4_detrended_q05 = bio4_rr$lo,
    bio4_detrended_q95 = bio4_rr$hi,
    bio4_detrended_range_q05_q95 = bio4_rr$range,
    median_doy_clean = median_doy_clean,
    global_lat_alpha = GLOBAL_ALPHA,
    global_lat_beta = GLOBAL_BETA,
    expected_doy_at_median_lat = expected_doy_at_median_lat,
    lat_standardized_doy_index = lat_standardized_doy_index,
    start_lat = df$start_lat[1],
    start_doy = df$start_doy[1],
    lat_beta = lat_beta,
    elev_beta = elev_beta,
    hopkins_lat_slope = HOPKINS_LAT_SLOPE,
    hopkins_elev_slope = HOPKINS_ELEV_SLOPE,
    lat_r2 = lat_r2,
    lat_elev_r2 = lat_elev_r2,
    hopkins_lat_elev_rmse = hop_rmse
  )
}

species_summary <- dat_obs |>
  group_by(scientificName) |>
  group_modify(~ fit_one_species_summary(.x)) |>
  ungroup()

# ------------------------------------------------------------------------------
# 9) CHECKS + WRITE
# ------------------------------------------------------------------------------

stopif_missing(all(c("pp_gate", "dpp_gate", "cell_id_25km") %in% names(dat_obs)),
               "Required obs-level fields missing from output.")

cat("\n=== GDD ZERO CHECKS ===\n")
cat("Obs zeros, current early GDD: ",
    sum(dat_obs$gddC_early == 0, na.rm = TRUE), " / ",
    sum(is.finite(dat_obs$gddC_early)),
    " (", round(100 * mean(dat_obs$gddC_early == 0, na.rm = TRUE), 1), "%)\n", sep = "")

cat("Obs zeros, M15 early GDD:     ",
    sum(dat_obs$gddC_early_m15 == 0, na.rm = TRUE), " / ",
    sum(is.finite(dat_obs$gddC_early_m15)),
    " (", round(100 * mean(dat_obs$gddC_early_m15 == 0, na.rm = TRUE), 1), "%)\n", sep = "")

readr::write_csv(dat_obs, OUT_OBS)
readr::write_csv(species_summary, OUT_SPECIES)

cat("\n=== COMPLETE: OBS-ONLY BUILD ===\n")
cat("Obs-level file:  ", OUT_OBS, "\n", sep = "")
cat("Species summary: ", OUT_SPECIES, "\n", sep = "")
cat("GDD windows:\n")
cat("  EARLY:     Jan 1 -> ", sprintf("%02d/%02d", EARLY_END_MONTH, EARLY_END_DAY), "\n", sep = "")
cat("  LATE :     ", sprintf("%02d/%02d", LATE_START_MONTH, LATE_START_DAY),
    " -> ", sprintf("%02d/%02d", LATE_END_MONTH, LATE_END_DAY), "\n", sep = "")
cat("  EARLY_M15: Jan 1 -> ", sprintf("%02d/%02d", EARLY_END_MONTH_M15, EARLY_END_DAY_M15), "\n", sep = "")
cat("  LATE_M15 : ", sprintf("%02d/%02d", LATE_START_MONTH_M15, LATE_START_DAY_M15),
    " -> ", sprintf("%02d/%02d", LATE_END_MONTH_M15, LATE_END_DAY_M15), "\n", sep = "")
cat("Canonical photoperiod gate used downstream: ", PP_GATE_LABEL, "\n", sep = "")
cat("Canonical photoperiod gate used downstream: ", PP_GATE_LABEL, "\n", sep = "")