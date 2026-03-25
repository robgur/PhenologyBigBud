# ==============================================================================
# FIGURE 6 (REVISED): 2-PANEL SYNTHESIS MAP
# A = predicted Hopkins deviation
# B = sensitivity to early forcing
# Compact layout with panel-local legends
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(grid)
  library(maps)
})

MAP_YEAR <- 2019

LAT_MIN <- 21
LAT_MAX <- 68

MAP_DIR <- sprintf("map_outputs_hfp_21to68N_%d", MAP_YEAR)

SENS_PATH <- file.path(
  MAP_DIR,
  sprintf("S_Early_danom_dgddE_sc_%d.tif", MAP_YEAR)
)

ANOM_PATH <- file.path(
  MAP_DIR,
  sprintf("Predicted_anom_%d.tif", MAP_YEAR)
)

OUT_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_FIG <- file.path(OUT_DIR, "Figure6_synthesis_maps_2panel_compact.png")

need <- function(path) if (!file.exists(path)) stop("Missing: ", path, call. = FALSE)
need(SENS_PATH)
need(ANOM_PATH)

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

plot_df <- function(r, nm = "value") {
  d <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(d)[3] <- nm
  
  if (max(d$x, na.rm = TRUE) > 180) {
    d$x <- ifelse(d$x > 180, d$x - 360, d$x)
  }
  
  d %>%
    filter(is.finite(y), y >= LAT_MIN, y <= LAT_MAX) %>%
    arrange(y, x)
}

robust_sym_limit <- function(x, prob = 0.98) {
  unname(quantile(abs(x), prob, na.rm = TRUE))
}

# ----------------------------------------------------------------------
# Load rasters
# ----------------------------------------------------------------------

r_sens <- rast(SENS_PATH)
r_anom <- rast(ANOM_PATH)

# convert so positive = stronger advancement / earlier than Hopkins
df_sens <- plot_df(r_sens, "value") %>%
  mutate(value = -value)

df_anom <- plot_df(r_anom, "value") %>%
  mutate(value = -value)

sens_lim <- robust_sym_limit(df_sens$value, prob = 0.98)
anom_lim <- robust_sym_limit(df_anom$value, prob = 0.98)

# shared x extent from actual data
x_rng <- range(c(df_sens$x, df_anom$x), na.rm = TRUE)
x_pad <- 0.01 * diff(x_rng)
xlim_use <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)

world_df <- map_data("world") %>%
  filter(
    long >= xlim_use[1] - 5, long <= xlim_use[2] + 5,
    lat >= LAT_MIN - 5, lat <= LAT_MAX + 5
  )

# ----------------------------------------------------------------------
# Styling
# ----------------------------------------------------------------------

theme_map <- theme_void(base_size = 11) +
  theme(
    panel.background    = element_rect(fill = "#8f8f8f", colour = NA),
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.border        = element_rect(colour = "grey35", fill = NA, linewidth = 0.45),
    plot.title          = element_text(
      size = 12.5, face = "bold", hjust = 0,
      margin = margin(b = 0.5)
    ),
    plot.title.position = "panel",
    plot.margin         = margin(0, 2, 0, 2),
    legend.title        = element_text(size = 8.8, lineheight = 0.95),
    legend.text         = element_text(size = 7.8),
    legend.key.height   = unit(0.48, "cm"),
    legend.key.width    = unit(0.18, "cm"),
    legend.margin       = margin(0, 0, 0, 0),
    legend.box.margin   = margin(0, 0, 0, 0),
    legend.spacing.y    = unit(2, "pt")
  )

scale_div <- function(name, lim) {
  scale_fill_gradient2(
    low = "#D73027",
    mid = "white",
    high = "#3B4CC0",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = squish,
    name = name,
    breaks = pretty(c(-lim, lim), n = 5)
  )
}

make_map <- function(df, title, legend_name, lim) {
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = value), na.rm = FALSE) +
    geom_path(
      data = world_df,
      aes(long, lat, group = group),
      colour = "grey80",
      linewidth = 0.18
    ) +
    coord_quickmap(
      xlim = xlim_use,
      ylim = c(LAT_MIN, LAT_MAX),
      expand = FALSE
    ) +
    scale_div(legend_name, lim) +
    labs(title = title) +
    theme_map +
    theme(
      legend.position = "right",
      legend.justification = "center"
    )
}

# ----------------------------------------------------------------------
# Build panels
# ----------------------------------------------------------------------

pA <- make_map(
  df_anom,
  "A. Predicted Hopkins deviation",
  "Earlier than\nHopkins",
  anom_lim
)

pB <- make_map(
  df_sens,
  "B. Sensitivity to early forcing",
  "Advancement\nsensitivity",
  sens_lim
)

fig <- pA / pB +
  plot_layout(heights = c(1, 1)) &
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(
  OUT_FIG,
  fig,
  width = 9.2,
  height = 4.4,
  dpi = 500,
  bg = "white"
)

print(fig)