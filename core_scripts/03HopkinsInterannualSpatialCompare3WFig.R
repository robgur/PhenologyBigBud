# ==============================================================================
# 03_interannual_direction_test_obsOnly_21to68N.R
#
# Goal:
#   Test whether TEMPORAL (interannual) responses have the same direction as
#   SPATIAL responses implied by the selected observation model.
#
# Uses:
#   - filtered obs CSV from hopkins_outputs_v11_gatePP_obsOnly21to68N
#   - selected spatial anomaly model from Script 02b
#
# Core idea:
#   Build species×cell×year panels, within-panel center predictors and response,
#   fit temporal anomaly model, then compare predicted +/- 1 SD changes against
#   the spatial anomaly model using the same delta-prediction logic.
#
# Outputs:
#   - temporal model summary
#   - temporal partial plot
#   - temporal direction-check csv
#   - spatial vs temporal comparison csv
#   - spatial vs temporal comparison plot
#   - Figure 4 empirical temporal-consistency plot
#
# IMPORTANT:
#   This script uses the canonical gate predictor exported by Script 01:
#     - dpp_gate
#   so temporal tests remain aligned with the main spatial model.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(splines)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------

LAT_MIN <- 21
LAT_MAX <- 68

INPUT_DIR <- "hopkins_outputs_v11_gatePP_obsOnly21to68N"
MODEL_DIR <- "models_gatePP_21to68N_model_selection"
OUT_DIR   <- "interannual_tests_21to68N"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBS_FILE <- file.path(INPUT_DIR, "hopkins_obs_with_covariates_plus_gatePP.csv")
SPATIAL_ANOM_MODEL_FILE <- file.path(MODEL_DIR, "m_anom_best.rds")

OUT_FIG4 <- file.path(OUT_DIR, "Figure4_interannual_test.png")

stopif_missing <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)

stopif_missing(file.exists(OBS_FILE), paste0("Missing obs file: ", OBS_FILE))
stopif_missing(
  file.exists(SPATIAL_ANOM_MODEL_FILE),
  paste0("Missing spatial anomaly model: ", SPATIAL_ANOM_MODEL_FILE)
)

# ------------------------------------------------------------------------------
# 1) Load data + spatial model
# ------------------------------------------------------------------------------

dat_obs <- readr::read_csv(OBS_FILE, show_col_types = FALSE) %>%
  filter(
    is.finite(latitude),
    latitude >= LAT_MIN,
    latitude <= LAT_MAX
  )

m_spatial_anom <- readRDS(SPATIAL_ANOM_MODEL_FILE)

stopif_missing(
  "dpp_gate" %in% names(dat_obs),
  "OBS file must contain canonical gate column 'dpp_gate'."
)

# determine year factor reference level for spatial prediction
year_levels_obs <- sort(unique(dat_obs$year))

# choose a reference year for spatial delta predictions
REF_YEAR <- 2019
if (!(REF_YEAR %in% year_levels_obs)) REF_YEAR <- year_levels_obs[1]

# ------------------------------------------------------------------------------
# 2) Build species×cell×year panel for mean timing anomalies
# ------------------------------------------------------------------------------

panel <- dat_obs %>%
  filter(
    is.finite(anom),
    is.finite(gddC_early), is.finite(gddC_late),
    is.finite(tmeanC_Jan),
    is.finite(dpp_gate),
    is.finite(bio4_lat_detrended),
    !is.na(scientificName), !is.na(cell_id_25km), !is.na(year)
  ) %>%
  group_by(scientificName, cell_id_25km, year) %>%
  summarise(
    n_obs = n(),
    anom_mean = mean(anom, na.rm = TRUE),
    gddE = mean(gddC_early, na.rm = TRUE),
    gddL = mean(gddC_late,  na.rm = TRUE),
    jan  = mean(tmeanC_Jan, na.rm = TRUE),
    dpp  = mean(dpp_gate, na.rm = TRUE),
    bio4 = mean(bio4_lat_detrended, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(scientificName, cell_id_25km) %>%
  mutate(n_years = n()) %>%
  ungroup() %>%
  filter(n_years >= 3, n_obs >= 3)

panel <- panel %>%
  mutate(
    gddE_sc = as.numeric(scale(gddE)),
    gddL_sc = as.numeric(scale(gddL)),
    jan_sc  = as.numeric(scale(jan)),
    dpp_sc  = as.numeric(scale(dpp)),
    bio4_sc = as.numeric(scale(bio4)),
    year_sc = as.numeric(scale(year))
  ) %>%
  group_by(scientificName, cell_id_25km) %>%
  mutate(
    gddE_dev = gddE_sc - mean(gddE_sc, na.rm = TRUE),
    gddL_dev = gddL_sc - mean(gddL_sc, na.rm = TRUE),
    jan_dev  = jan_sc  - mean(jan_sc,  na.rm = TRUE),
    year_dev = year_sc - mean(year_sc, na.rm = TRUE),
    anom_dev = anom_mean - mean(anom_mean, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(w = sqrt(pmin(n_obs, 25)))

# ------------------------------------------------------------------------------
# 3) Fit temporal anomaly model
# ------------------------------------------------------------------------------

m_temporal_anom_bio4 <- lm(
  anom_dev ~
    ns(gddL_dev, 2) * bio4_sc +
    ns(gddE_dev, 2) * dpp_sc * bio4_sc +
    jan_dev +
    year_dev,
  data = panel,
  weights = w
)

sink(file.path(OUT_DIR, "temporal_model_anom_summary.txt"))
cat("=== TEMPORAL MODEL (anom_dev) ===\n")
print(summary(m_temporal_anom_bio4))
sink()

# ------------------------------------------------------------------------------
# 4) Temporal partial plot
# ------------------------------------------------------------------------------

panel <- panel %>%
  mutate(
    dpp_bin = cut(
      dpp_sc,
      breaks = quantile(dpp_sc, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("1", "2", "3")
    )
  )

p_anom_gddE <- ggplot(panel, aes(gddE_dev, anom_dev)) +
  geom_point(alpha = 0.18) +
  geom_smooth(method = "lm", se = FALSE, color = "royalblue", linewidth = 1.3) +
  facet_wrap(~ dpp_bin, nrow = 1) +
  labs(
    title = "Temporal test: within species×cell deviations (anom_dev vs gddE_dev)",
    x = "gddE_dev (within-panel, scaled)",
    y = "anom_dev (within-panel)"
  ) +
  theme_bw(base_size = 14)

ggsave(
  filename = file.path(OUT_DIR, "temporal_partial_anom_gddE.png"),
  plot = p_anom_gddE, width = 14, height = 4.8, dpi = 300
)

# ------------------------------------------------------------------------------
# 5) Delta-prediction helper
# ------------------------------------------------------------------------------

predict_delta <- function(model, base_row, var, delta = 1) {
  nd0 <- base_row
  nd1 <- base_row
  nd1[[var]] <- nd1[[var]] + delta
  
  if (inherits(model, "lmerMod")) {
    p1 <- predict(model, newdata = nd1, re.form = NA, allow.new.levels = TRUE)
    p0 <- predict(model, newdata = nd0, re.form = NA, allow.new.levels = TRUE)
  } else {
    p1 <- predict(model, newdata = nd1)
    p0 <- predict(model, newdata = nd0)
  }
  
  as.numeric(p1 - p0)
}

# ------------------------------------------------------------------------------
# 6) Direction checks for temporal model
# ------------------------------------------------------------------------------

dpp_levels  <- c(-1, 0, 1)
bio4_levels <- c(-1, 0, 1)

context_grid <- expand.grid(
  dpp_sc  = dpp_levels,
  bio4_sc = bio4_levels
)

base_anom <- data.frame(
  gddL_dev = 0,
  gddE_dev = 0,
  jan_dev  = 0,
  year_dev = 0,
  dpp_sc   = 0,
  bio4_sc  = 0
)

dir_anom <- lapply(seq_len(nrow(context_grid)), function(i) {
  br <- base_anom
  br$dpp_sc  <- context_grid$dpp_sc[i]
  br$bio4_sc <- context_grid$bio4_sc[i]
  tibble(
    dpp_sc  = br$dpp_sc,
    bio4_sc = br$bio4_sc,
    delta_anom_gddE = predict_delta(m_temporal_anom_bio4, br, "gddE_dev", delta = 1),
    delta_anom_gddL = predict_delta(m_temporal_anom_bio4, br, "gddL_dev", delta = 1),
    delta_anom_jan  = predict_delta(m_temporal_anom_bio4, br, "jan_dev",  delta = 1)
  )
}) %>% bind_rows()

write_csv(dir_anom, file.path(OUT_DIR, "direction_check_temporal_anom.csv"))

# ------------------------------------------------------------------------------
# 7) Spatial direction checks using same context grid
# ------------------------------------------------------------------------------

base_spatial_anom <- data.frame(
  gddL_sc  = 0,
  gddE_sc  = 0,
  jan_sc   = 0,
  dpp_sc   = 0,
  bio4_sc  = 0,
  year_f = factor(REF_YEAR, levels = sort(unique(dat_obs$year))),
  cell_id_25km = panel$cell_id_25km[1],
  scientificName = panel$scientificName[1]
)

spatial_dir_anom <- lapply(seq_len(nrow(context_grid)), function(i) {
  br <- base_spatial_anom
  br$dpp_sc  <- context_grid$dpp_sc[i]
  br$bio4_sc <- context_grid$bio4_sc[i]
  tibble(
    dpp_sc  = br$dpp_sc,
    bio4_sc = br$bio4_sc,
    delta_anom_gddE = predict_delta(m_spatial_anom, br, "gddE_sc", delta = 1),
    delta_anom_gddL = predict_delta(m_spatial_anom, br, "gddL_sc", delta = 1),
    delta_anom_jan  = predict_delta(m_spatial_anom, br, "jan_sc",  delta = 1)
  )
}) %>% bind_rows()

write_csv(spatial_dir_anom, file.path(OUT_DIR, "direction_check_spatial_anom.csv"))

# ------------------------------------------------------------------------------
# 8) Compare spatial vs temporal direction explicitly
# ------------------------------------------------------------------------------

comparison_anom <- spatial_dir_anom %>%
  rename(
    spatial_gddE = delta_anom_gddE,
    spatial_gddL = delta_anom_gddL,
    spatial_jan  = delta_anom_jan
  ) %>%
  left_join(
    dir_anom %>%
      rename(
        temporal_gddE = delta_anom_gddE,
        temporal_gddL = delta_anom_gddL,
        temporal_jan  = delta_anom_jan
      ),
    by = c("dpp_sc", "bio4_sc")
  ) %>%
  mutate(
    agree_gddE = sign(spatial_gddE) == sign(temporal_gddE),
    agree_gddL = sign(spatial_gddL) == sign(temporal_gddL),
    agree_jan  = sign(spatial_jan)  == sign(temporal_jan)
  )

write_csv(comparison_anom, file.path(OUT_DIR, "spatial_temporal_comparison_anom.csv"))

# ------------------------------------------------------------------------------
# 9) Summary table of sign agreement
# ------------------------------------------------------------------------------

anom_agreement_summary <- tibble(
  variable = c("gddE", "gddL", "jan"),
  n_agree = c(
    sum(comparison_anom$agree_gddE, na.rm = TRUE),
    sum(comparison_anom$agree_gddL, na.rm = TRUE),
    sum(comparison_anom$agree_jan,  na.rm = TRUE)
  ),
  n_total = nrow(comparison_anom)
) %>%
  mutate(prop_agree = n_agree / n_total)

write_csv(anom_agreement_summary, file.path(OUT_DIR, "agreement_summary_anom.csv"))

# ------------------------------------------------------------------------------
# 10) Spatial vs temporal comparison plot
# ------------------------------------------------------------------------------

comparison_anom_long <- comparison_anom %>%
  select(dpp_sc, bio4_sc,
         spatial_gddE, temporal_gddE,
         spatial_gddL, temporal_gddL,
         spatial_jan, temporal_jan) %>%
  pivot_longer(
    cols = -c(dpp_sc, bio4_sc),
    names_to = c("type", "variable"),
    names_sep = "_",
    values_to = "effect"
  ) %>%
  pivot_wider(names_from = type, values_from = effect)

p_compare_anom <- ggplot(comparison_anom_long, aes(spatial, temporal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
  geom_point(size = 2) +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Spatial vs temporal direction check: anomaly model",
    x = "Spatial predicted delta",
    y = "Temporal predicted delta"
  ) +
  theme_bw(base_size = 13)

ggsave(
  filename = file.path(OUT_DIR, "spatial_vs_temporal_anom.png"),
  plot = p_compare_anom, width = 9, height = 4.8, dpi = 300
)

# ------------------------------------------------------------------------------
# 11) Figure 4 empirical temporal consistency plot
# ------------------------------------------------------------------------------

fig4_panel <- dat_obs %>%
  filter(
    is.finite(anom),
    is.finite(gddC_early),
    is.finite(tmeanC_Jan),
    is.finite(dpp_gate),
    is.finite(bio4_lat_detrended),
    !is.na(scientificName),
    !is.na(cell_id_25km),
    !is.na(year)
  ) %>%
  group_by(scientificName, cell_id_25km, year) %>%
  summarise(
    n_obs = n(),
    anom_mean = mean(anom, na.rm = TRUE),
    gddE = mean(gddC_early, na.rm = TRUE),
    jan  = mean(tmeanC_Jan, na.rm = TRUE),
    dpp  = mean(dpp_gate, na.rm = TRUE),
    bio4 = mean(bio4_lat_detrended, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(scientificName, cell_id_25km) %>%
  mutate(n_years = n()) %>%
  ungroup() %>%
  filter(n_years >= 3, n_obs >= 3) %>%
  mutate(
    gddE_sc = as.numeric(scale(gddE)),
    jan_sc  = as.numeric(scale(jan)),
    dpp_sc  = as.numeric(scale(dpp)),
    bio4_sc = as.numeric(scale(bio4)),
    year_sc = as.numeric(scale(year))
  ) %>%
  group_by(scientificName, cell_id_25km) %>%
  mutate(
    gddE_dev = gddE_sc - mean(gddE_sc, na.rm = TRUE),
    jan_dev  = jan_sc  - mean(jan_sc,  na.rm = TRUE),
    year_dev = year_sc - mean(year_sc, na.rm = TRUE),
    anom_dev = anom_mean - mean(anom_mean, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    dpp_bin = cut(
      dpp_sc,
      breaks = quantile(dpp_sc, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c(
        "Slow photoperiod progression",
        "Intermediate photoperiod progression",
        "Fast photoperiod progression"
      )
    )
  )

line_color <- "#3A5FCD"

theme_fig2_style <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    strip.background = element_rect(fill = "grey95", colour = "grey35", linewidth = 0.6),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 15, face = "plain", hjust = 0),
    panel.border = element_rect(colour = "grey35")
  )

p_interannual <- ggplot(fig4_panel, aes(gddE_dev, anom_dev)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    alpha = 0.16,
    size = 1.5,
    colour = "black"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = line_color,
    linewidth = 1.35
  ) +
  facet_wrap(~ dpp_bin, nrow = 1) +
  coord_cartesian(xlim = c(-1.3, 1.05)) +
  labs(
    x = "Early spring forcing deviation",
    y = "Mean Hopkins anomaly deviation"
  ) +
  theme_fig2_style

ggsave(
  filename = OUT_FIG4,
  plot = p_interannual,
  width = 10.6,
  height = 4.2,
  dpi = 400,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 12) Finish
# ------------------------------------------------------------------------------

cat("\n✓ Finished interannual anomaly direction test.\n")
cat("Input obs file: ", OBS_FILE, "\n", sep = "")
cat("Spatial model:  ", SPATIAL_ANOM_MODEL_FILE, "\n", sep = "")
cat("Outputs in:     ", OUT_DIR, "\n", sep = "")