# ==============================================================================
# 06_helper_plots_hfp_21to68N.R
#
# Helper plots for mapped HFP observation-model outputs
#
# Expects outputs from:
#   05_map_obs_hfp_gatePP_21to68N.R
#
# Reads:
#   map_outputs_hfp_21to68N_<MAP_YEAR>/
#     - S_Early_danom_dgddE_sc_<MAP_YEAR>.tif
#     - GateStrength_dpp_<MAP_YEAR>.tif
#     - S_Late_danom_dgddL_sc_<MAP_YEAR>.tif   (optional)
#
# Also uses:
#   - original bio4 raster, aligned here on the fly to the map output geometry
#
# Produces:
#   quick map plots
#   gate-vs-seasonality scatter/smoother plots
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------

MAP_YEAR <- 2019

LAT_MIN <- 21
LAT_MAX <- 68
LON_MIN <- -180
LON_MAX <- 180
LL_EXT  <- terra::ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX)

OUT_DIR <- sprintf("map_outputs_hfp_21to68N_%d", MAP_YEAR)

S_EARLY_PATH <- file.path(OUT_DIR, sprintf("S_Early_danom_dgddE_sc_%d.tif", MAP_YEAR))
S_LATE_PATH  <- file.path(OUT_DIR, sprintf("S_Late_danom_dgddL_sc_%d.tif", MAP_YEAR))
GATE_PATH    <- file.path(OUT_DIR, sprintf("GateStrength_dpp_%d.tif", MAP_YEAR))

# Original source raster, not cached/aligned
BIO4_PATH <- "bio4_lat_elev_detrended_NH.tif"

stopif_missing <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)
need <- function(path) stopif_missing(file.exists(path), paste0("Missing: ", path))

need(S_EARLY_PATH)
need(GATE_PATH)
need(BIO4_PATH)

HAS_LATE <- file.exists(S_LATE_PATH)

# ------------------------------------------------------------------------------
# 1) Helpers
# ------------------------------------------------------------------------------

plot_df <- function(r, nm = "value", max_n = 500000) {
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- nm
  
  # wrap longitudes to -180..180 if needed
  if (max(d$x, na.rm = TRUE) > 180) {
    d$x <- ifelse(d$x > 180, d$x - 360, d$x)
  }
  
  d <- d %>%
    filter(is.finite(y), y >= LAT_MIN, y <= LAT_MAX)
  
  if (nrow(d) > max_n) {
    set.seed(1)
    d <- d[sample(seq_len(nrow(d)), max_n), ]
  }
  
  d
}

align_static_to_template <- function(r, template, ext_ll) {
  r <- tryCatch(
    crop(r, ext_ll),
    error = function(e) r
  )
  
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    if (!same.crs(r, template)) {
      r <- project(r, template, method = "bilinear")
    } else {
      r <- resample(r, template, method = "bilinear")
    }
  }
  
  mask(r, template)
}

# ------------------------------------------------------------------------------
# 2) Load rasters
# ------------------------------------------------------------------------------

S_Early <- rast(S_EARLY_PATH)
GateStrength <- rast(GATE_PATH)

if (HAS_LATE) {
  S_Late <- rast(S_LATE_PATH)
}

# Build aligned bio4 directly from original source raster
bio4_raw <- rast(BIO4_PATH)
bio4 <- align_static_to_template(bio4_raw, template = GateStrength, ext_ll = LL_EXT)

# ------------------------------------------------------------------------------
# 3) Quick map data frames
# ------------------------------------------------------------------------------

df_early <- plot_df(S_Early, "S_Early")
df_gate  <- plot_df(GateStrength, "GateStrength")

if (HAS_LATE) {
  df_late <- plot_df(S_Late, "S_Late")
}

# For interpretability: positive = stronger advancement from warming
df_early <- df_early %>%
  mutate(S_Early_plot = -S_Early)

if (HAS_LATE) {
  df_late <- df_late %>%
    mutate(S_Late_plot = -S_Late)
}

# ------------------------------------------------------------------------------
# 4) Quick helper maps
# ------------------------------------------------------------------------------

p_early <- ggplot(df_early, aes(x, y, color = S_Early_plot)) +
  geom_point(size = 0.2) +
  coord_fixed() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Early\nforcing sensitivity"
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    title = "Sensitivity of phenological anomaly to early forcing",
    subtitle = "More positive values indicate stronger advancement under early warming"
  ) +
  theme_bw()

p_gate <- ggplot(df_gate, aes(x, y, color = GateStrength)) +
  geom_point(size = 0.2) +
  coord_fixed() +
  scale_color_gradient2(
    low = "purple", mid = "white", high = "orange", midpoint = 0,
    name = "Gate\nstrength"
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    title = "Strength of the photoperiod gate on early-forcing sensitivity",
    subtitle = "Large absolute values indicate stronger dependence of early-forcing effects on photoperiod environment"
  ) +
  theme_bw()

print(p_early)
print(p_gate)

if (HAS_LATE) {
  p_late <- ggplot(df_late, aes(x, y, color = S_Late_plot)) +
    geom_point(size = 0.2) +
    coord_fixed() +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Late\nforcing sensitivity"
    ) +
    labs(
      x = "Longitude", y = "Latitude",
      title = "Sensitivity of phenological anomaly to late forcing",
      subtitle = "More positive values indicate stronger advancement under late warming"
    ) +
    theme_bw()
  
  print(p_late)
}

# ------------------------------------------------------------------------------
# 5) Early-forcing control classes
# ------------------------------------------------------------------------------

classify_temp_control <- function(x) {
  case_when(
    x <= -20 ~ "strong temperature control",
    x <= -5  ~ "moderate temperature control",
    x > -5   ~ "weak temperature control"
  )
}

df_early_class <- df_early %>%
  mutate(control_class = classify_temp_control(S_Early))

p_class <- ggplot(df_early_class, aes(x, y, color = control_class)) +
  geom_point(size = 0.2) +
  coord_fixed() +
  labs(
    x = "Longitude", y = "Latitude",
    title = "Where is temperature a dominant driver of phenological timing?"
  ) +
  theme_bw()

print(p_class)

# ------------------------------------------------------------------------------
# 6) Gate strength vs relative seasonality
# ------------------------------------------------------------------------------

combo <- c(GateStrength, bio4)
names(combo) <- c("gate_strength", "bio4_det")

df_combo <- as.data.frame(combo, xy = TRUE, na.rm = TRUE) %>%
  filter(is.finite(y), y >= LAT_MIN, y <= LAT_MAX)

if (nrow(df_combo) > 300000) {
  set.seed(1)
  df_combo <- df_combo[sample(seq_len(nrow(df_combo)), 300000), ]
}

# Keep support window broad enough to retain structure, but avoid extreme tails
xlim_use <- quantile(df_combo$bio4_det, probs = c(0.01, 0.99), na.rm = TRUE)
ylim_use <- quantile(df_combo$gate_strength, probs = c(0.005, 0.995), na.rm = TRUE)

df_smooth <- df_combo %>%
  filter(
    bio4_det >= xlim_use[1],
    bio4_det <= xlim_use[2]
  )

# Palette-consistent blue
curve_col <- "#4C78A8"

p_fig2 <- ggplot(df_combo, aes(bio4_det, gate_strength)) +
  geom_point(
    color = "grey35",
    alpha = 0.02,
    size = 0.08
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.4
  ) +
  geom_smooth(
    data = df_smooth,
    method = "gam",
    formula = y ~ s(x, k = 8),
    color = curve_col,
    fill = "grey80",
    linewidth = 1.35,
    se = TRUE
  ) +
  coord_cartesian(
    xlim = xlim_use,
    ylim = c(ylim_use[1], ylim_use[2] + 1.2),
    expand = FALSE
  ) +
  labs(
    x = "Temperature seasonality relative to latitude and elevation",
    y = "Photoperiod gate strength",
    title = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

print(p_fig2)

ggsave(
  filename = file.path(OUT_DIR, "Figure_gate_vs_seasonality_clean.png"),
  plot = p_fig2,
  width = 7.2,
  height = 5.2,
  dpi = 500,
  bg = "white"
)

cat("\n✓ Helper plots written to: ", OUT_DIR, "\n", sep = "")