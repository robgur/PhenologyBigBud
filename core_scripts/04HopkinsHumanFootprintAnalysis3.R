# ==============================================================================
# 04_fit_models_hfp_gatePP_21to68N.R
# Observation-only: Human Footprint (HFP) as modifier of cue environment and
# forcing sensitivity
#
# IMPORTANT:
#   This version uses the canonical photoperiod-gate variable exported by Script 01:
#     - dpp_gate
#
# Inputs:
#   hopkins_outputs_v11_gatePP_obsOnly21to68N/hopkins_obs_with_covariates_plus_gatePP.csv
#
# Outputs (folder: models_hfp_gatePP_21to68N/):
#   - obs_model_comparison_hfp.csv
#   - m_anom_hfp_best.rds
#   - dat2_obs_with_hfp.rds
#   - scalers_obs_best.rds
#   - train_quantiles_obs_best.rds
#   - model_summaries_hfp.txt
#   - Figure_HFP_onepanel.png
#
# Model-selection framework:
#   Observation candidates:
#     M0 = baseline spatial model (no HFP)
#     M1 = + additive HFP shift
#     M2 = + general forcing moderation by HFP
#     M3 = HFP modifies the EARLY forcing gate
#
# Final model:
#   Hard-coded REML refit of M3, because that is the chosen interpretive model.
#
# Notes:
#   - HFP is assumed to already be present in the obs CSV as column `hfp`
#   - No raster extraction is performed here
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(splines)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(grid)
})

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------

LAT_MIN <- 21
LAT_MAX <- 68

OBS_FILE <- "hopkins_outputs_v11_gatePP_obsOnly21to68N/hopkins_obs_with_covariates_plus_gatePP.csv"

OUT_DIR <- "models_hfp_gatePP_21to68N"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_FIG <- file.path(OUT_DIR, "Figure_HFP_onepanel.png")

stopif_missing <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)

stopif_missing(file.exists(OBS_FILE), paste0("Missing obs file: ", OBS_FILE))

q2 <- function(x, probs = c(0.01, 0.99)) {
  as.numeric(quantile(x, probs = probs, na.rm = TRUE))
}

ms <- function(x) list(mu = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))

pred_fixef <- function(model, newdata) {
  predict(model, newdata = newdata, re.form = NA, allow.new.levels = TRUE)
}

# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------

dat_obs <- readr::read_csv(OBS_FILE, show_col_types = FALSE) %>%
  filter(
    is.finite(latitude),
    latitude >= LAT_MIN,
    latitude <= LAT_MAX
  )

# ------------------------------------------------------------------------------
# 2) Prepare HFP from precomputed obs column
# ------------------------------------------------------------------------------

needed_obs <- c(
  "anom", "cell_id_25km", "scientificName", "year",
  "gddC_early", "gddC_late", "tmeanC_Jan", "dpp_gate", "bio4_lat_detrended",
  "hfp"
)

stopif_missing(
  all(needed_obs %in% names(dat_obs)),
  paste0(
    "OBS file missing required columns:\n",
    paste(setdiff(needed_obs, names(dat_obs)), collapse = ", ")
  )
)

dat_obs <- dat_obs %>%
  mutate(
    hfp_sc = as.numeric(scale(hfp))
  )

cat("HFP finite values:", sum(is.finite(dat_obs$hfp)), "of", nrow(dat_obs), "\n")
cat("HFP scaled finite values:", sum(is.finite(dat_obs$hfp_sc)), "of", nrow(dat_obs), "\n")

# ------------------------------------------------------------------------------
# 3) Prepare OBS training frame
# ------------------------------------------------------------------------------

dat2 <- dat_obs %>%
  mutate(
    gddE_sc = as.numeric(scale(gddC_early)),
    gddL_sc = as.numeric(scale(gddC_late)),
    jan_sc  = as.numeric(scale(tmeanC_Jan)),
    dpp_sc  = as.numeric(scale(dpp_gate)),
    bio4_sc = as.numeric(scale(bio4_lat_detrended)),
    year_f  = factor(year)
  ) %>%
  filter(
    complete.cases(
      anom, gddE_sc, gddL_sc, jan_sc, dpp_sc, bio4_sc, hfp_sc,
      year_f, cell_id_25km, scientificName
    )
  )

cat("Rows in dat2:", nrow(dat2), "\n")
stopif_missing(nrow(dat2) > 0, "dat2 has zero rows after filtering.")

# ------------------------------------------------------------------------------
# 4) Observation candidate models with HFP
# ------------------------------------------------------------------------------

obs_formula_strings <- list(
  obs_M0_base =
    "anom ~
       ns(gddL_sc, 2) * bio4_sc +
       ns(gddE_sc, 2) * dpp_sc * bio4_sc +
       jan_sc + year_f +
       (1 | cell_id_25km) +
       (1 | scientificName)",
  
  obs_M1_hfp_additive =
    "anom ~
       ns(gddL_sc, 2) * bio4_sc +
       ns(gddE_sc, 2) * dpp_sc * bio4_sc +
       jan_sc + year_f +
       hfp_sc +
       (1 | cell_id_25km) +
       (1 | scientificName)",
  
  obs_M2_hfp_force =
    "anom ~
       ns(gddL_sc, 2) * bio4_sc +
       ns(gddE_sc, 2) * dpp_sc * bio4_sc +
       jan_sc + year_f +
       hfp_sc +
       ns(gddL_sc, 2):hfp_sc +
       ns(gddE_sc, 2):hfp_sc +
       (1 | cell_id_25km) +
       (1 | scientificName)",
  
  obs_M3_hfp_gate =
    "anom ~
       ns(gddL_sc, 2) * bio4_sc +
       ns(gddE_sc, 2) * dpp_sc * bio4_sc * hfp_sc +
       jan_sc + year_f +
       (1 | cell_id_25km) +
       (1 | scientificName)"
)

message("Fitting observation HFP candidate models (ML)...")

obs_candidates <- lapply(names(obs_formula_strings), function(nm) {
  message("  fitting: ", nm)
  lmer(as.formula(obs_formula_strings[[nm]]), data = dat2, REML = FALSE)
})
names(obs_candidates) <- names(obs_formula_strings)

obs_comp <- AIC(
  obs_candidates$obs_M0_base,
  obs_candidates$obs_M1_hfp_additive,
  obs_candidates$obs_M2_hfp_force,
  obs_candidates$obs_M3_hfp_gate
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC) %>%
  mutate(deltaAIC = AIC - min(AIC))

write_csv(obs_comp, file.path(OUT_DIR, "obs_model_comparison_hfp.csv"))

# ------------------------------------------------------------------------------
# 5) Final observation model with HFP (hard-coded chosen model)
# ------------------------------------------------------------------------------

m_anom_hfp_best <- lmer(
  anom ~
    ns(gddL_sc, 2) * bio4_sc +
    ns(gddE_sc, 2) * dpp_sc * bio4_sc * hfp_sc +
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  data = dat2,
  REML = TRUE
)

obs_selected_name <- "obs_M3_hfp_gate"

saveRDS(m_anom_hfp_best, file.path(OUT_DIR, "m_anom_hfp_best.rds"))
saveRDS(dat2, file.path(OUT_DIR, "dat2_obs_with_hfp.rds"))

# ------------------------------------------------------------------------------
# 6) Save best observation model artifacts
# ------------------------------------------------------------------------------

scalers_obs_best <- list(
  gddC_early         = ms(dat2$gddC_early),
  gddC_late          = ms(dat2$gddC_late),
  tmeanC_Jan         = ms(dat2$tmeanC_Jan),
  dpp_gate           = ms(dat2$dpp_gate),
  bio4_lat_detrended = ms(dat2$bio4_lat_detrended),
  hfp                = ms(dat2$hfp)
)

saveRDS(scalers_obs_best, file.path(OUT_DIR, "scalers_obs_best.rds"))

train_quantiles_obs_best <- list(
  gddE_sc = q2(dat2$gddE_sc),
  gddL_sc = q2(dat2$gddL_sc),
  jan_sc  = q2(dat2$jan_sc),
  dpp_sc  = q2(dat2$dpp_sc),
  bio4_sc = q2(dat2$bio4_sc),
  hfp_sc  = q2(dat2$hfp_sc)
)

saveRDS(train_quantiles_obs_best, file.path(OUT_DIR, "train_quantiles_obs_best.rds"))

sink(file.path(OUT_DIR, "model_summaries_hfp.txt"))
cat("=== OBS HFP MODEL COMPARISON ===\n")
print(obs_comp)
cat("\n=== FINAL OBSERVATION HFP MODEL ===\n")
cat("Selected observation HFP model: ", obs_selected_name, "\n\n", sep = "")
print(summary(m_anom_hfp_best))
sink()

# ------------------------------------------------------------------------------
# 7) Single-row three-panel plot: HFP effect on forcing sensitivity
# ------------------------------------------------------------------------------

gddE_low  <- as.numeric(quantile(dat2$gddE_sc, probs = 0.25, na.rm = TRUE))
gddE_high <- as.numeric(quantile(dat2$gddE_sc, probs = 0.75, na.rm = TRUE))

hfp_levels <- as.numeric(quantile(dat2$hfp_sc, probs = c(0.25, 0.75), na.rm = TRUE))
bio4_slices <- c(-1, 0, 1)

year_levels <- levels(dat2$year_f)
ref_year <- if ("2019" %in% year_levels) "2019" else year_levels[1]

gate_grid <- expand.grid(
  dpp_sc  = seq(-2, 2, length.out = 200),
  bio4_sc = bio4_slices,
  hfp_sc  = hfp_levels
) %>%
  mutate(
    gddL_sc = 0,
    jan_sc  = 0,
    year_f = factor(ref_year, levels = year_levels),
    cell_id_25km   = dat2$cell_id_25km[which(!is.na(dat2$cell_id_25km))[1]],
    scientificName = dat2$scientificName[which(!is.na(dat2$scientificName))[1]]
  )

gate_low  <- gate_grid %>% mutate(gddE_sc = gddE_low)
gate_high <- gate_grid %>% mutate(gddE_sc = gddE_high)

plot_df <- gate_grid %>%
  mutate(
    gate_strength = pred_fixef(m_anom_hfp_best, gate_high) - pred_fixef(m_anom_hfp_best, gate_low),
    seasonality = factor(
      bio4_sc,
      levels = bio4_slices,
      labels = c("Low seasonality", "Intermediate", "High seasonality")
    ),
    hfp_label = factor(
      hfp_sc,
      levels = hfp_levels,
      labels = c("Low HFP", "High HFP")
    )
  )

pal_bio4 <- c(
  "Low seasonality" = "#C53E35",
  "Intermediate"    = "#4C78A8",
  "High seasonality" = "#54A24B"
)

theme_fig2 <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey86", linewidth = 0.35),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "plain", size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6)
    )
}

p_hfp <- ggplot(
  plot_df,
  aes(x = dpp_sc, y = gate_strength, color = seasonality, linetype = hfp_label)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ seasonality, nrow = 1) +
  scale_color_manual(
    values = pal_bio4,
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("solid", "22"),
    name = "Human footprint"
  ) +
  labs(
    title = "",
    x = "Photoperiod progression rate",
    y = "Change in Hopkins deviation\n(high minus low early forcing)"
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_fig2() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(2, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

ggsave(
  filename = OUT_FIG,
  plot = p_hfp,
  width = 10.2,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

cat("\n✓ Saved HFP observation model outputs to: ", OUT_DIR, "\n", sep = "")
print(p_hfp)