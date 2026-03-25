# ==============================================================================
# 05_map_obs_hfp_gatePP_21to68N_FINAL.R
#
# Map observation-model surfaces from the fitted HFP timing model
#
# Produces:
#   - S_Early_danom_dgddE_sc_<MAP_YEAR>.tif
#   - S_Late_danom_dgddL_sc_<MAP_YEAR>.tif        (optional)
#   - GateStrength_dpp_<MAP_YEAR>.tif
#   - Predicted_anom_<MAP_YEAR>.tif
#
# Core model inputs:
#   gddE_sc, gddL_sc, jan_sc, bio4_sc, hfp_sc, dpp_sc, year_f
#
# Notes:
#   - Domain is cropped to 21–68 N
#   - HFP is aggregated before projection for speed
#   - year_f is restored as a factor during prediction
#   - ERA longitude is normalized to -180..180 before cropping
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(lme4)
})

# ------------------------------------------------------------------------------
# 0) USER SETTINGS
# ------------------------------------------------------------------------------

MAP_YEAR <- 2019
DELTA_SCALED <- 0.05
DO_LATE <- TRUE

LAT_MIN <- 21
LAT_MAX <- 68
LON_MIN <- -180
LON_MAX <- 180
LL_EXT  <- terra::ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX)

MODEL_PATH   <- "models_hfp_gatePP_21to68N/m_anom_hfp_best.rds"
SCALERS_PATH <- "models_hfp_gatePP_21to68N/scalers_obs_best.rds"
Q_PATH       <- "models_hfp_gatePP_21to68N/train_quantiles_obs_best.rds"

ERA5_PATH   <- "era5land_t2m_2017_2025_combined.nc"
ERA_VAR     <- "t2m"
ORIGIN_DATE <- as.Date("2017-01-01")

BIO4_PATH <- "bio4_lat_elev_detrended_NH.tif"
HFP_PATH  <- "hfp_2020_100m_v1-2_cog.tif"

HFP_AGG_FACTOR <- 4
USE_HFP_CACHE  <- TRUE

GDD_BASE_C <- 5
EARLY_END  <- "02-28"
LATE_START <- "03-01"
LATE_END   <- "04-30"

DPP_REF_DOY <- 60

OUT_DIR <- sprintf("map_outputs_hfp_21to68N_%d", MAP_YEAR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------------------------

stopif_missing <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)
need <- function(path) stopif_missing(file.exists(path), paste0("Missing: ", path))

need(MODEL_PATH)
need(SCALERS_PATH)
need(ERA5_PATH)
need(BIO4_PATH)
need(HFP_PATH)

clamp_r <- function(r, lo, hi) {
  r2 <- r
  r2 <- ifel(r2 < lo, lo, r2)
  r2 <- ifel(r2 > hi, hi, r2)
  r2
}

crop_ll_raster <- function(r, ext_ll) {
  tryCatch(
    terra::crop(r, ext_ll),
    error = function(e) r
  )
}

# Fix global 0..360 longitude rasters before cropping to -180..180 extents
normalize_lon_global <- function(r) {
  ex <- terra::ext(r)
  
  # rotate only if it looks like a lon/lat global raster on 0..360
  if (terra::is.lonlat(r) && ex$xmax > 180) {
    message("Rotating global lon/lat raster from 0..360 to -180..180...")
    r <- terra::rotate(r)
  }
  r
}

crop_projected_raster_to_ll_domain <- function(r, ext_ll) {
  poly_ll <- terra::as.polygons(ext_ll, crs = "EPSG:4326")
  poly_r  <- terra::project(poly_ll, terra::crs(r))
  terra::crop(r, poly_r)
}

daylength_hours <- function(lat_deg, doy) {
  lat_rad <- lat_deg * pi / 180
  decl <- 0.409 * sin(2 * pi * (doy - 81) / 365)
  x <- -tan(lat_rad) * tan(decl)
  x <- pmin(pmax(x, -1), 1)
  ha <- acos(x)
  24 * ha / pi
}

d_daylength_hours <- function(lat_deg, doy) {
  daylength_hours(lat_deg, doy) - daylength_hours(lat_deg, pmax(doy - 1, 1))
}

gdd_from_layers <- function(x, base = 5) {
  sum(pmax(x - base, 0), na.rm = TRUE)
}

predict_merMod <- function(model, stk, year_levels, ref_year) {
  f <- function(model, data, ...) {
    d <- as.data.frame(data)
    out <- rep(NA_real_, nrow(d))
    ok <- complete.cases(d)
    if (any(ok)) {
      nd <- d[ok, , drop = FALSE]
      nd$year_f <- factor(ref_year, levels = year_levels)
      
      out[ok] <- as.numeric(
        predict(
          model,
          newdata = nd,
          re.form = NA,
          allow.new.levels = TRUE
        )
      )
    }
    out
  }
  terra::predict(stk, model = model, fun = f)
}

fd_deriv <- function(model, build_stack, x_sc, year_levels, ref_year, delta_sc = 0.05) {
  p_plus  <- predict_merMod(
    model,
    build_stack(x_sc + delta_sc),
    year_levels = year_levels,
    ref_year = ref_year
  )
  p_minus <- predict_merMod(
    model,
    build_stack(x_sc - delta_sc),
    year_levels = year_levels,
    ref_year = ref_year
  )
  (p_plus - p_minus) / (2 * delta_sc)
}

# ------------------------------------------------------------------------------
# 2) LOAD MODEL + SCALERS
# ------------------------------------------------------------------------------

m <- readRDS(MODEL_PATH)
sc <- readRDS(SCALERS_PATH)
q  <- if (file.exists(Q_PATH)) readRDS(Q_PATH) else NULL

stopif_missing(inherits(m, "merMod"),
               "MODEL_PATH must be an lmer merMod object.")

stopif_missing(
  all(c("gddC_early", "gddC_late", "tmeanC_Jan",
        "bio4_lat_detrended", "dpp_gate", "hfp") %in% names(sc)),
  paste(
    "scalers_obs_best.rds must contain:",
    "gddC_early, gddC_late, tmeanC_Jan, bio4_lat_detrended, dpp_gate, hfp"
  )
)

year_levels <- levels(model.frame(m)$year_f)
stopif_missing(length(year_levels) > 0,
               "Could not recover year_f levels from fitted model.")

ref_year <- as.character(MAP_YEAR)
if (!ref_year %in% year_levels) {
  warning("MAP_YEAR is not a year_f level in the model; using first available level instead.")
  ref_year <- year_levels[1]
}

# ------------------------------------------------------------------------------
# 3) BUILD YEAR-SPECIFIC CLIMATE RASTERS FROM ERA5
# ------------------------------------------------------------------------------

message("Loading ERA5 and building annual climate rasters...")
era <- terra::rast(ERA5_PATH, subds = ERA_VAR)

# CRITICAL FIX: normalize longitude BEFORE cropping to the -180..180 domain
era <- normalize_lon_global(era)
era <- crop_ll_raster(era, LL_EXT)

all_dates <- seq.Date(ORIGIN_DATE, by = "day", length.out = nlyr(era))

idx_between <- function(start_date, end_date) {
  which(all_dates >= start_date & all_dates <= end_date)
}

y0 <- as.Date(sprintf("%d-01-01", MAP_YEAR))
y1 <- as.Date(sprintf("%d-12-31", MAP_YEAR))
stopif_missing(any(all_dates >= y0 & all_dates <= y1),
               "MAP_YEAR not covered by ERA5 date range.")

idx_jan   <- idx_between(as.Date(sprintf("%d-01-01", MAP_YEAR)),
                         as.Date(sprintf("%d-01-31", MAP_YEAR)))
idx_early <- idx_between(as.Date(sprintf("%d-01-01", MAP_YEAR)),
                         as.Date(sprintf("%d-%s", MAP_YEAR, EARLY_END)))
idx_late  <- idx_between(as.Date(sprintf("%d-%s", MAP_YEAR, LATE_START)),
                         as.Date(sprintf("%d-%s", MAP_YEAR, LATE_END)))

template <- era[[idx_jan[1]]]

tjan_K <- terra::app(era[[idx_jan]], fun = mean, na.rm = TRUE)
tjan_C <- tjan_K - 273.15

early_C <- era[[idx_early]] - 273.15
gddE_raw <- terra::app(early_C, fun = gdd_from_layers, base = GDD_BASE_C)

late_C <- era[[idx_late]] - 273.15
gddL_raw <- terra::app(late_C, fun = gdd_from_layers, base = GDD_BASE_C)

# ------------------------------------------------------------------------------
# 4) LOAD + ALIGN STATIC RASTERS
# ------------------------------------------------------------------------------

message("Preparing static rasters...")

# BIO4
bio4_raw <- rast(BIO4_PATH)
bio4_raw <- normalize_lon_global(bio4_raw)
bio4_raw <- crop_ll_raster(bio4_raw, LL_EXT)

if (!compareGeom(bio4_raw, template, stopOnError = FALSE)) {
  if (!same.crs(bio4_raw, template)) {
    bio4_raw <- project(bio4_raw, template, method = "bilinear")
  } else {
    bio4_raw <- resample(bio4_raw, template, method = "bilinear")
  }
}

# HFP
hfp_cache_native <- file.path(CACHE_DIR, sprintf("hfp_native_crop_agg%d.tif", HFP_AGG_FACTOR))
hfp_cache_template <- file.path(CACHE_DIR, sprintf("hfp_on_template_agg%d.tif", HFP_AGG_FACTOR))

rebuild_hfp <- TRUE

if (USE_HFP_CACHE && file.exists(hfp_cache_template)) {
  message("Found cached HFP aligned to template; checking geometry...")
  hfp_test <- rast(hfp_cache_template)
  
  if (compareGeom(hfp_test, template, stopOnError = FALSE)) {
    message("Cached HFP matches current template.")
    hfp_raw <- hfp_test
    rebuild_hfp <- FALSE
  } else {
    message("Cached HFP does not match current template; rebuilding.")
  }
}

if (rebuild_hfp) {
  if (USE_HFP_CACHE && file.exists(hfp_cache_native)) {
    message("Found cached native cropped/aggregated HFP; checking/reusing...")
    hfp_raw <- rast(hfp_cache_native)
  } else {
    message("Cropping HFP in native CRS...")
    hfp_raw <- rast(HFP_PATH)
    hfp_raw <- crop_projected_raster_to_ll_domain(hfp_raw, LL_EXT)
    
    message("Aggregating HFP by factor ", HFP_AGG_FACTOR, " before projection...")
    hfp_raw <- aggregate(hfp_raw, fact = HFP_AGG_FACTOR, fun = mean, na.rm = TRUE)
    
    if (USE_HFP_CACHE) {
      writeRaster(hfp_raw, hfp_cache_native, overwrite = TRUE)
    }
  }
  
  if (!compareGeom(hfp_raw, template, stopOnError = FALSE)) {
    message("Projecting aggregated HFP directly to template...")
    hfp_raw <- project(hfp_raw, template, method = "bilinear")
  }
  
  if (USE_HFP_CACHE) {
    writeRaster(hfp_raw, hfp_cache_template, overwrite = TRUE)
  }
}

# DPP
lat_r <- init(template, "y")
dpp_raw <- app(lat_r, fun = function(lat) d_daylength_hours(lat, DPP_REF_DOY))

# Common mask
stopif_missing(compareGeom(gddE_raw, template, stopOnError = FALSE), "gddE_raw does not match template")
stopif_missing(compareGeom(gddL_raw, template, stopOnError = FALSE), "gddL_raw does not match template")
stopif_missing(compareGeom(tjan_C,  template, stopOnError = FALSE), "tjan_C does not match template")
stopif_missing(compareGeom(bio4_raw, template, stopOnError = FALSE), "bio4_raw does not match template")
stopif_missing(compareGeom(hfp_raw,  template, stopOnError = FALSE), "hfp_raw does not match template")
stopif_missing(compareGeom(dpp_raw,  template, stopOnError = FALSE), "dpp_raw does not match template")

ok_mask <- is.finite(gddE_raw) & is.finite(gddL_raw) & is.finite(tjan_C) &
  is.finite(bio4_raw) & is.finite(hfp_raw) & is.finite(dpp_raw)
# ------------------------------------------------------------------------------
# 5) SCALE PREDICTORS
# ------------------------------------------------------------------------------

gddE_sc <- (gddE_raw - sc$gddC_early$mu) / sc$gddC_early$sd
gddL_sc <- (gddL_raw - sc$gddC_late$mu) / sc$gddC_late$sd
jan_sc  <- (tjan_C   - sc$tmeanC_Jan$mu) / sc$tmeanC_Jan$sd
bio4_sc <- (bio4_raw - sc$bio4_lat_detrended$mu) / sc$bio4_lat_detrended$sd
hfp_sc  <- (hfp_raw  - sc$hfp$mu) / sc$hfp$sd
dpp_sc  <- (dpp_raw  - sc$dpp_gate$mu) / sc$dpp_gate$sd

mask01 <- ifel(ok_mask, 1, NA)

gddE_sc <- mask(gddE_sc, mask01)
gddL_sc <- mask(gddL_sc, mask01)
jan_sc  <- mask(jan_sc,  mask01)
bio4_sc <- mask(bio4_sc, mask01)
hfp_sc  <- mask(hfp_sc,  mask01)
dpp_sc  <- mask(dpp_sc,  mask01)

# ------------------------------------------------------------------------------
# 6) OPTIONAL CLAMPING IN SCALED SPACE
# ------------------------------------------------------------------------------

if (!is.null(q)) {
  get_lohi <- function(name, fallback = c(-3, 3)) {
    if (!name %in% names(q)) return(fallback)
    v <- as.numeric(q[[name]])
    if (length(v) >= 2 && all(is.finite(v[1:2]))) c(v[1], v[2]) else fallback
  }
  
  g <- get_lohi("gddE_sc"); gddE_sc <- clamp_r(gddE_sc, g[1], g[2])
  g <- get_lohi("gddL_sc"); gddL_sc <- clamp_r(gddL_sc, g[1], g[2])
  g <- get_lohi("jan_sc");  jan_sc  <- clamp_r(jan_sc,  g[1], g[2])
  g <- get_lohi("bio4_sc"); bio4_sc <- clamp_r(bio4_sc, g[1], g[2])
  g <- get_lohi("hfp_sc");  hfp_sc  <- clamp_r(hfp_sc,  g[1], g[2])
  g <- get_lohi("dpp_sc");  dpp_sc  <- clamp_r(dpp_sc,  g[1], g[2])
}

# ------------------------------------------------------------------------------
# 7) BUILD PREDICTION STACKS
# ------------------------------------------------------------------------------

build_stack <- function(gddE_override = gddE_sc,
                        gddL_override = gddL_sc,
                        dpp_override  = dpp_sc) {
  stk <- c(gddE_override, gddL_override, jan_sc, bio4_sc, hfp_sc, dpp_override)
  names(stk) <- c("gddE_sc", "gddL_sc", "jan_sc", "bio4_sc", "hfp_sc", "dpp_sc")
  
  yr <- template
  values(yr) <- 1
  names(yr) <- "year_f"
  
  stk <- c(stk, yr)
  stk <- mask(stk, mask01)
  stk
}

# ------------------------------------------------------------------------------
# 8) COMPUTE PREDICTION SURFACES
# ------------------------------------------------------------------------------

message("Computing predicted anomaly surface...")
Pred_anom <- predict_merMod(
  m,
  build_stack(),
  year_levels = year_levels,
  ref_year = ref_year
)
Pred_anom <- mask(Pred_anom, mask01)

message("Computing S_Early = d(anom)/d(gddE_sc) ...")
S_E_dpp0 <- fd_deriv(
  m,
  build_stack = function(gE) build_stack(gddE_override = gE, dpp_override = dpp_sc),
  x_sc = gddE_sc,
  year_levels = year_levels,
  ref_year = ref_year,
  delta_sc = DELTA_SCALED
)

if (DO_LATE) {
  message("Computing S_Late = d(anom)/d(gddL_sc) ...")
  S_L <- fd_deriv(
    m,
    build_stack = function(gL) build_stack(gddL_override = gL, dpp_override = dpp_sc),
    x_sc = gddL_sc,
    year_levels = year_levels,
    ref_year = ref_year,
    delta_sc = DELTA_SCALED
  )
} else {
  S_L <- NULL
}

message("Computing GateStrength wrt dpp: S_E(dpp=+1) - S_E(dpp=-1) ...")
S_E_dppP1 <- fd_deriv(
  m,
  build_stack = function(gE) build_stack(gddE_override = gE, dpp_override = dpp_sc * 0 + 1),
  x_sc = gddE_sc,
  year_levels = year_levels,
  ref_year = ref_year,
  delta_sc = DELTA_SCALED
)

S_E_dppM1 <- fd_deriv(
  m,
  build_stack = function(gE) build_stack(gddE_override = gE, dpp_override = dpp_sc * 0 - 1),
  x_sc = gddE_sc,
  year_levels = year_levels,
  ref_year = ref_year,
  delta_sc = DELTA_SCALED
)

GateStrength <- S_E_dppP1 - S_E_dppM1

Pred_anom    <- mask(Pred_anom, mask01)
S_E_dpp0     <- mask(S_E_dpp0, mask01)
if (!is.null(S_L)) S_L <- mask(S_L, mask01)
GateStrength <- mask(GateStrength, mask01)

# ------------------------------------------------------------------------------
# 9) WRITE OUTPUTS
# ------------------------------------------------------------------------------

writeRaster(
  Pred_anom,
  file.path(OUT_DIR, sprintf("Predicted_anom_%d.tif", MAP_YEAR)),
  overwrite = TRUE
)

writeRaster(
  S_E_dpp0,
  file.path(OUT_DIR, sprintf("S_Early_danom_dgddE_sc_%d.tif", MAP_YEAR)),
  overwrite = TRUE
)

if (!is.null(S_L)) {
  writeRaster(
    S_L,
    file.path(OUT_DIR, sprintf("S_Late_danom_dgddL_sc_%d.tif", MAP_YEAR)),
    overwrite = TRUE
  )
}

writeRaster(
  GateStrength,
  file.path(OUT_DIR, sprintf("GateStrength_dpp_%d.tif", MAP_YEAR)),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 10) REPORT
# ------------------------------------------------------------------------------

cat("\nWrote outputs to: ", OUT_DIR, "\n", sep = "")
cat("  - Predicted_anom_", MAP_YEAR, ".tif\n", sep = "")
cat("  - S_Early_danom_dgddE_sc_", MAP_YEAR, ".tif\n", sep = "")
if (DO_LATE) cat("  - S_Late_danom_dgddL_sc_", MAP_YEAR, ".tif\n", sep = "")
cat("  - GateStrength_dpp_", MAP_YEAR, ".tif\n", sep = "")
cat("Domain cropped to: ", LAT_MIN, "–", LAT_MAX, " N\n", sep = "")
cat("Photoperiod progression reference date: DOY ", DPP_REF_DOY, "\n", sep = "")
cat("year_f level used for prediction: ", ref_year, "\n", sep = "")
cat("HFP aggregated before projection by factor: ", HFP_AGG_FACTOR, "\n", sep = "")