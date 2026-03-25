# ==============================================================================
# 09-FIGURE: SAMPLING SUPPORT ACROSS GEOGRAPHY, SPECIES, AND CLIMATE SPACE
# Updated for:
#   - filtered 21–68 N obs-only inputs
#   - 2025 ERA coverage
#   - corrected full-longitude climate-space availability
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(terra)
  library(scales)
  library(stringr)
  library(maps)
  library(gridExtra)
  library(grid)
  library(hexbin)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 0) SETTINGS
# ------------------------------------------------------------------------------

LAT_MIN <- 21
LAT_MAX <- 68

# Updated upstream obs-only filtered file
OBS_FILE <- "hopkins_outputs_v11_gatePP_obsOnly21to68N/hopkins_obs_with_covariates_plus_gatePP.csv"

# No change needed here if the anchored seasonality layer filename is unchanged
BIO4_PATH <- "bio4_lat_elev_detrended_NH.tif"

ERA5_PATH <- "era5land_t2m_2017_2025_combined.nc"
ERA_VAR   <- "t2m"
ORIGIN_DATE <- as.Date("2017-01-01")

OUT_DIR <- "figures_sampling"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_FIG_PNG <- file.path(OUT_DIR, "Fig_sampling_support_3panel_relativeSupport_21to68N.png")
OUT_FIG_PDF <- file.path(OUT_DIR, "Fig_sampling_support_3panel_relativeSupport_21to68N.pdf")

GRID_KM <- 100
GRID_M  <- GRID_KM * 1000
EQAREA_CRS <- "EPSG:6933"

N_SHOW_SPECIES <- 40
MAP_YEAR <- 2025

GDD_BASE_C <- 5
EARLY_END_MONTH <- 2
EARLY_END_DAY   <- 28

MAX_N_GRID_CLIMATE <- 250000
MAX_N_OBS_CLIMATE  <- 120000

# ---- climate-space panel controls ----
CLIM_BINS_X <- 55
CLIM_BINS_Y <- 45
EPS_REL_SUPPORT <- 1e-8

# axis trimming for display only
X_PROBS <- c(0.01, 0.99)
Y_PROBS <- c(0.01, 0.995)

# ---- typography / layout controls ----
FIG_WIDTH  <- 11.0
FIG_HEIGHT <- 8.6
FIG_DPI    <- 500

TITLE_SIZE        <- 21
PANEL_TITLE_SIZE  <- 16
PANEL_SUB_SIZE    <- 11.5
AXIS_TITLE_SIZE   <- 14
AXIS_TEXT_SIZE    <- 11.5
LEGEND_TITLE_SIZE <- 12.5
LEGEND_TEXT_SIZE  <- 11
SPECIES_TEXT_SIZE <- 9.2

SHOW_EVERY_OTHER_SPECIES_LABEL <- FALSE

# ------------------------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------------------------

stopif_missing <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)

idx_between <- function(all_dates, start_date, end_date) {
  which(all_dates >= start_date & all_dates <= end_date)
}

gdd_from_layers <- function(x, base = 5) {
  sum(pmax(x - base, 0), na.rm = TRUE)
}

normalize_lon_global <- function(r) {
  ex <- terra::ext(r)
  if (terra::is.lonlat(r) && ex$xmax > 180) {
    r <- terra::rotate(r)
  }
  r
}

crop_ll_raster <- function(r, xmin = -180, xmax = 180, ymin = LAT_MIN, ymax = LAT_MAX) {
  terra::crop(r, terra::ext(xmin, xmax, ymin, ymax))
}

theme_panel <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = PANEL_TITLE_SIZE),
    plot.subtitle = element_text(size = PANEL_SUB_SIZE, margin = margin(b = 6)),
    axis.title = element_text(size = AXIS_TITLE_SIZE),
    axis.text = element_text(size = AXIS_TEXT_SIZE),
    legend.title = element_text(size = LEGEND_TITLE_SIZE),
    legend.text = element_text(size = LEGEND_TEXT_SIZE),
    plot.margin = margin(4, 8, 4, 8)
  )

# ------------------------------------------------------------------------------
# 2) LOAD OBSERVATIONS
# ------------------------------------------------------------------------------

stopif_missing(file.exists(OBS_FILE), paste0("Missing: ", OBS_FILE))
stopif_missing(file.exists(BIO4_PATH), paste0("Missing: ", BIO4_PATH))
stopif_missing(file.exists(ERA5_PATH), paste0("Missing: ", ERA5_PATH))

dat <- read_csv(OBS_FILE, show_col_types = FALSE) %>%
  filter(
    is.finite(longitude),
    is.finite(latitude),
    latitude >= LAT_MIN,
    latitude <= LAT_MAX,
    !is.na(scientificName)
  ) %>%
  mutate(
    scientificName = str_squish(scientificName)
  )

n_total_obs <- nrow(dat)
n_total_species <- n_distinct(dat$scientificName)

# ------------------------------------------------------------------------------
# 3) PANEL A: GEOGRAPHIC SAMPLING DENSITY
# ------------------------------------------------------------------------------

v <- vect(dat[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
v_eq <- project(v, EQAREA_CRS)
xy <- crds(v_eq)

dat_map <- dat %>%
  mutate(
    x_m = xy[, 1],
    y_m = xy[, 2],
    cell_x = floor(x_m / GRID_M),
    cell_y = floor(y_m / GRID_M)
  )

grid_counts <- dat_map %>%
  group_by(cell_x, cell_y) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(scientificName),
    .groups = "drop"
  ) %>%
  mutate(
    x_center = (cell_x + 0.5) * GRID_M,
    y_center = (cell_y + 0.5) * GRID_M
  )

v_cent <- vect(grid_counts[, c("x_center", "y_center")], geom = c("x_center", "y_center"), crs = EQAREA_CRS)
v_cent_ll <- project(v_cent, "EPSG:4326")
xy_ll <- crds(v_cent_ll)

grid_counts <- grid_counts %>%
  mutate(
    lon = xy_ll[, 1],
    lat = xy_ll[, 2],
    log10_records = log10(n_records)
  )

world_df <- ggplot2::map_data("world") %>%
  filter(lat >= LAT_MIN, lat <= LAT_MAX)

tile_width_deg  <- 1.5
tile_height_deg <- 1.1

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(long, lat, group = group),
    fill = "#f7f7f7", color = "#c3c3c3", linewidth = 0.20
  ) +
  geom_tile(
    data = grid_counts,
    aes(lon, lat, fill = log10_records),
    width = tile_width_deg,
    height = tile_height_deg
  ) +
  coord_fixed(xlim = c(-180, 180), ylim = c(LAT_MIN, LAT_MAX), expand = FALSE) +
  scale_fill_gradientn(
    colours = c("#edf8fb", "#bfd3e6", "#8c96c6", "#8856a7", "#810f7c"),
    values = rescale(c(0, 0.15, 0.4, 0.7, 1)),
    breaks = c(0, 1, 2, 3),
    labels = c("1", "10", "100", "1,000"),
    name = "Records per 100-km cell"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(5.0, "cm"),
      barheight = unit(0.65, "cm")
    )
  ) +
  labs(
    title = "A. Geographic sampling density",
    subtitle = paste0(comma(n_total_obs), " observations from ", comma(n_total_species), " species")
  ) +
  theme_panel +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  )

# ------------------------------------------------------------------------------
# 4) PANEL B: LATITUDINAL SUPPORT ACROSS SPECIES
# ------------------------------------------------------------------------------

species_support <- dat %>%
  group_by(scientificName) %>%
  summarise(
    n_total = n(),
    lat_q05 = quantile(latitude, 0.05, na.rm = TRUE),
    lat_q25 = quantile(latitude, 0.25, na.rm = TRUE),
    lat_med = median(latitude, na.rm = TRUE),
    lat_q75 = quantile(latitude, 0.75, na.rm = TRUE),
    lat_q95 = quantile(latitude, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_total)) %>%
  slice_head(n = N_SHOW_SPECIES) %>%
  arrange(desc(lat_med), desc(n_total)) %>%
  mutate(
    species_plot = factor(scientificName, levels = rev(scientificName)),
    alpha_val = rescale(log10(n_total), to = c(0.40, 1))
  )

if (SHOW_EVERY_OTHER_SPECIES_LABEL) {
  sp_levels <- levels(species_support$species_plot)
  sp_labels <- ifelse(seq_along(sp_levels) %% 2 == 1, sp_levels, "")
} else {
  sp_levels <- levels(species_support$species_plot)
  sp_labels <- sp_levels
}

p_lat <- ggplot(species_support, aes(y = species_plot)) +
  geom_segment(
    aes(x = lat_q05, xend = lat_q95, yend = species_plot),
    linewidth = 1.05, color = "grey70", lineend = "round"
  ) +
  geom_segment(
    aes(x = lat_q25, xend = lat_q75, yend = species_plot, alpha = alpha_val),
    linewidth = 3.2, color = "#bf6b17", lineend = "round"
  ) +
  geom_point(
    aes(x = lat_med, alpha = alpha_val),
    size = 2.7, shape = 21, stroke = 0.20,
    fill = "#8c4f12", color = "white"
  ) +
  scale_alpha_continuous(range = c(0.45, 1), guide = "none") +
  scale_x_continuous(
    limits = c(LAT_MIN, LAT_MAX),
    breaks = seq(20, 70, by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_discrete(labels = sp_labels) +
  labs(
    title = "B. Species latitudinal range",
    subtitle = paste0("Top ", N_SHOW_SPECIES, " species by sample size"),
    x = "Latitude (°N)",
    y = NULL
  ) +
  theme_panel +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey89", linewidth = 0.35),
    axis.text.y = element_text(size = SPECIES_TEXT_SIZE, face = "italic"),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# 5) PANEL C: RELATIVE OBSERVATIONAL SUPPORT IN CLIMATE SPACE
# ------------------------------------------------------------------------------

era <- rast(ERA5_PATH, subds = ERA_VAR)
era <- normalize_lon_global(era)
era <- crop_ll_raster(era, ymin = LAT_MIN, ymax = LAT_MAX)

all_dates <- seq.Date(ORIGIN_DATE, by = "day", length.out = nlyr(era))

idx_early <- idx_between(
  all_dates,
  as.Date(sprintf("%d-01-01", MAP_YEAR)),
  as.Date(sprintf("%d-%02d-%02d", MAP_YEAR, EARLY_END_MONTH, EARLY_END_DAY))
)

template <- era[[idx_early[1]]]

early_C <- era[[idx_early]] - 273.15
gddE_r <- app(early_C, fun = gdd_from_layers, base = GDD_BASE_C)

bio4_r <- rast(BIO4_PATH)
bio4_r <- normalize_lon_global(bio4_r)
bio4_r <- crop_ll_raster(bio4_r, ymin = LAT_MIN, ymax = LAT_MAX)

if (!compareGeom(bio4_r, template, stopOnError = FALSE)) {
  if (!same.crs(bio4_r, template)) bio4_r <- project(bio4_r, template, method = "bilinear")
  bio4_r <- resample(bio4_r, template, method = "bilinear")
}

combo_r <- c(bio4_r, gddE_r)
names(combo_r) <- c("bio4", "gddE")

clim_all <- as.data.frame(combo_r, xy = FALSE, na.rm = TRUE)

if (nrow(clim_all) > MAX_N_GRID_CLIMATE) {
  set.seed(1)
  clim_all <- clim_all[sample(seq_len(nrow(clim_all)), MAX_N_GRID_CLIMATE), ]
}

obs_clim <- dat %>%
  filter(
    is.finite(bio4_lat_detrended),
    is.finite(gddC_early)
  ) %>%
  transmute(
    bio4 = bio4_lat_detrended,
    gddE = gddC_early
  )

if (nrow(obs_clim) > MAX_N_OBS_CLIMATE) {
  set.seed(2)
  obs_clim <- obs_clim[sample(seq_len(nrow(obs_clim)), MAX_N_OBS_CLIMATE), ]
}

combo_x <- c(clim_all$bio4, obs_clim$bio4)
combo_y <- c(clim_all$gddE, obs_clim$gddE)

x_lim <- quantile(combo_x, probs = X_PROBS, na.rm = TRUE)
y_lim <- quantile(combo_y, probs = Y_PROBS, na.rm = TRUE)

clim_all_trim <- clim_all %>%
  filter(
    bio4 >= x_lim[1], bio4 <= x_lim[2],
    gddE >= y_lim[1], gddE <= y_lim[2]
  )

obs_clim_trim <- obs_clim %>%
  filter(
    bio4 >= x_lim[1], bio4 <= x_lim[2],
    gddE >= y_lim[1], gddE <= y_lim[2]
  )

x_breaks <- seq(x_lim[1], x_lim[2], length.out = CLIM_BINS_X + 1)
y_breaks <- seq(y_lim[1], y_lim[2], length.out = CLIM_BINS_Y + 1)

x_mids <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2
y_mids <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2

bin2d_df <- function(df, x_breaks, y_breaks) {
  xb <- cut(df$bio4, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  yb <- cut(df$gddE, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  
  tibble(xbin = xb, ybin = yb) %>%
    filter(is.finite(xbin), is.finite(ybin)) %>%
    count(xbin, ybin, name = "n")
}

avail_bins <- bin2d_df(clim_all_trim, x_breaks, y_breaks) %>%
  rename(avail_n = n)

obs_bins <- bin2d_df(obs_clim_trim, x_breaks, y_breaks) %>%
  rename(obs_n = n)

clim_grid <- tidyr::expand_grid(
  xbin = seq_len(CLIM_BINS_X),
  ybin = seq_len(CLIM_BINS_Y)
) %>%
  left_join(avail_bins, by = c("xbin", "ybin")) %>%
  left_join(obs_bins,   by = c("xbin", "ybin")) %>%
  mutate(
    avail_n = ifelse(is.na(avail_n), 0L, avail_n),
    obs_n   = ifelse(is.na(obs_n),   0L, obs_n),
    xmid = x_mids[xbin],
    ymid = y_mids[ybin]
  )

avail_total <- sum(clim_grid$avail_n)
obs_total   <- sum(clim_grid$obs_n)

clim_grid <- clim_grid %>%
  mutate(
    avail_prop = avail_n / avail_total,
    obs_prop   = obs_n   / obs_total,
    rel_support = log10((obs_prop + EPS_REL_SUPPORT) / (avail_prop + EPS_REL_SUPPORT)),
    log10_avail = ifelse(avail_n > 0, log10(avail_n), NA_real_)
  )

overlay_df <- clim_grid %>%
  filter(obs_n > 0)

rel_lim <- quantile(abs(overlay_df$rel_support), probs = 0.98, na.rm = TRUE)
rel_lim <- max(rel_lim, 0.5)

p_clim <- ggplot() +
  geom_tile(
    data = clim_grid %>% filter(avail_n > 0),
    aes(x = xmid, y = ymid, fill = log10_avail),
    alpha = 0.82
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "grey45",
    na.value = "white",
    breaks = c(0, 1, 2, 3, 4),
    labels = c("1", "10", "100", "1000", "10000"),
    name = "Available NH cells"
  ) +
  geom_tile(
    data = overlay_df,
    aes(x = xmid, y = ymid, alpha = pmin(obs_n / max(obs_n), 1), color = rel_support),
    fill = NA,
    linewidth = 0.48
  ) +
  scale_color_gradient2(
    low = "#2166ac",
    mid = "#7a0177",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-rel_lim, rel_lim),
    oob = squish,
    name = "Relative support\n(obs vs available)"
  ) +
  scale_alpha(range = c(0.20, 0.95), guide = "none") +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(
    title = "C. Climate-space support",
    subtitle = "Grey = available climate space; colored cells = observational support relative to availability",
    x = "Temperature seasonality (detrended)",
    y = "Early forcing (GDD)"
  ) +
  theme_panel +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
    legend.position = "bottom"
  ) +
  guides(
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.8, "cm"),
      barheight = unit(0.55, "cm")
    ),
    color = guide_colorbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.8, "cm"),
      barheight = unit(0.55, "cm")
    )
  )

# ------------------------------------------------------------------------------
# 6) COMBINE + SAVE
# ------------------------------------------------------------------------------

title_grob <- textGrob(
  "Sampling support across geography, species, and climate space",
  x = 0, hjust = 0,
  gp = gpar(fontsize = TITLE_SIZE, fontface = "bold")
)

fig_grob <- arrangeGrob(
  title_grob,
  p_map,
  arrangeGrob(p_lat, p_clim, ncol = 2),
  ncol = 1,
  heights = c(0.075, 0.47, 0.50)
)

png(OUT_FIG_PNG, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_DPI)
grid.draw(fig_grob)
dev.off()

pdf(OUT_FIG_PDF, width = FIG_WIDTH, height = FIG_HEIGHT, useDingbats = FALSE)
grid.draw(fig_grob)
dev.off()

cat("Saved figures to:\n", OUT_FIG_PNG, "\n", OUT_FIG_PDF, "\n")