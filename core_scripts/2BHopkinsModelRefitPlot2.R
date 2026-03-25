# ==============================================================================
# 02b_select_best_models_gatePP_compare_M01_M15_21to68N.R
# Observation models only
#
# Purpose:
#   - remove synchrony analysis entirely
#   - fit 5 observation candidates for March 1 forcing windows
#   - fit 5 observation candidates for March 15 forcing windows
#   - print/save AIC tables for both
#   - hard-code REML refit of M4 for March 1
#   - produce the Figure 2 timing/gate-strength output using March 1 M4
#
# Inputs:
#   models_gatePP_21to68N/dat2_train_obs.rds
#
# Outputs:
#   - obs_model_comparison_AIC_March1.csv
#   - obs_model_comparison_AIC_March15.csv
#   - m_anom_best.rds   (hard-coded March 1 M4 REML refit)
#   - eff_anom_best.rds
#   - figure2_anom_plotdata.csv
#   - figures/Figure2_timing_and_gate_strength_compact.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(splines)
  library(ggeffects)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(grid)
})

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------

IN_DIR  <- "models_gatePP_21to68N"
OUT_DIR <- "models_gatePP_21to68N_model_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBS_TRAIN_FILE <- file.path(IN_DIR, "dat2_train_obs.rds")

need <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}
need(OBS_TRAIN_FILE)

# ------------------------------------------------------------------------------
# 1) Read training data
# ------------------------------------------------------------------------------

dat2_base <- readRDS(OBS_TRAIN_FILE)

if (!"year_f" %in% names(dat2_base)) {
  dat2_base <- dat2_base %>% mutate(year_f = factor(year))
} else {
  dat2_base <- dat2_base %>% mutate(year_f = factor(year_f))
}

# ------------------------------------------------------------------------------
# 2) Build March 1 and March 15 forcing-window datasets
# ------------------------------------------------------------------------------

dat2_m01 <- dat2_base %>%
  mutate(
    gddE_sc = as.numeric(scale(gddC_early)),
    gddL_sc = as.numeric(scale(gddC_late)),
    gddT_sc = as.numeric(scale(gddC_total))
  ) %>%
  filter(
    is.finite(anom),
    is.finite(gddE_sc),
    is.finite(gddL_sc),
    is.finite(gddT_sc),
    is.finite(jan_sc),
    is.finite(bio4_sc),
    is.finite(dpp_sc),
    !is.na(year_f),
    !is.na(cell_id_25km),
    !is.na(scientificName)
  )

dat2_m15 <- dat2_base %>%
  mutate(
    gddE_sc = as.numeric(scale(gddC_early_m15)),
    gddL_sc = as.numeric(scale(gddC_late_m15)),
    gddT_sc = as.numeric(scale(gddC_total_m15))
  ) %>%
  filter(
    is.finite(anom),
    is.finite(gddE_sc),
    is.finite(gddL_sc),
    is.finite(gddT_sc),
    is.finite(jan_sc),
    is.finite(bio4_sc),
    is.finite(dpp_sc),
    !is.na(year_f),
    !is.na(cell_id_25km),
    !is.na(scientificName)
  )

# ------------------------------------------------------------------------------
# 3) Candidate formulas (same model logic, obs only)
# ------------------------------------------------------------------------------

obs_formulas <- list(
  obs_M0_base =
    anom ~
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  
  obs_M1_late =
    anom ~
    ns(gddL_sc, 2) * bio4_sc +
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  
  obs_M2_early_gate =
    anom ~
    ns(gddE_sc, 2) * dpp_sc +
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  
  obs_M3_additive_gate =
    anom ~
    ns(gddL_sc, 2) * bio4_sc +
    ns(gddE_sc, 2) * dpp_sc +
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  
  obs_M4_full =
    anom ~
    ns(gddL_sc, 2) * bio4_sc +
    ns(gddE_sc, 2) * dpp_sc * bio4_sc +
    jan_sc + year_f +
    (1 | cell_id_25km) +
    (1 | scientificName)
)

# ------------------------------------------------------------------------------
# 4) Helper to fit one observation candidate set
# ------------------------------------------------------------------------------

fit_obs_set <- function(dat_use) {
  fits <- lapply(obs_formulas, function(fm) {
    lmer(fm, data = dat_use, REML = FALSE)
  })
  
  comp <- AIC(
    fits$obs_M0_base,
    fits$obs_M1_late,
    fits$obs_M2_early_gate,
    fits$obs_M3_additive_gate,
    fits$obs_M4_full
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("model") %>%
    arrange(AIC) %>%
    mutate(
      deltaAIC = AIC - min(AIC),
      npar = df
    )
  
  list(fits = fits, comp = comp)
}

# ------------------------------------------------------------------------------
# 5) Fit March 1 and March 15 candidate sets (ML)
# ------------------------------------------------------------------------------

message("Fitting March 1 observation candidate set...")
res_m01 <- fit_obs_set(dat2_m01)

message("Fitting March 15 observation candidate set...")
res_m15 <- fit_obs_set(dat2_m15)

obs_comp_m01 <- res_m01$comp
obs_comp_m15 <- res_m15$comp

print(obs_comp_m01)
print(obs_comp_m15)

write_csv(obs_comp_m01, file.path(OUT_DIR, "obs_model_comparison_AIC_March1.csv"))
write_csv(obs_comp_m15, file.path(OUT_DIR, "obs_model_comparison_AIC_March15.csv"))

# ------------------------------------------------------------------------------
# 6) Hard-code March 1 M4 REML refit
# ------------------------------------------------------------------------------

m_anom_best <- lmer(
  anom ~
    ns(gddL_sc, 2) * bio4_sc +
    ns(gddE_sc, 2) * dpp_sc * bio4_sc +
    jan_sc +
    year_f +
    (1 | cell_id_25km) +
    (1 | scientificName),
  data = dat2_m01,
  REML = TRUE
)

saveRDS(m_anom_best, file.path(OUT_DIR, "m_anom_best.rds"))
saveRDS(dat2_m01,    file.path(OUT_DIR, "dat2_train_obs.rds"))

sink(file.path(OUT_DIR, "best_model_summaries.txt"))
cat("=== OBS MODEL COMPARISON: MARCH 1 FORCING WINDOWS ===\n")
print(obs_comp_m01)
cat("\nHard-coded REML refit used for March 1: obs_M4_full\n\n")
print(summary(m_anom_best))

cat("\n\n=== OBS MODEL COMPARISON: MARCH 15 FORCING WINDOWS ===\n")
print(obs_comp_m15)
cat("\nExpected qualitative check: March 15 should also favor obs_M4_full if results are stable.\n")
sink()

# ------------------------------------------------------------------------------
# 7) Effect predictions for downstream figure generation
#    Hard-coded to M4 plotting path
# ------------------------------------------------------------------------------

eff_anom_best <- ggpredict(
  m_anom_best,
  terms = c(
    "gddE_sc [-2:2 by=0.05]",
    "dpp_sc [-1,0,1]",
    "bio4_sc [-1,0,1]"
  )
)

saveRDS(eff_anom_best, file.path(OUT_DIR, "eff_anom_best.rds"))

# ------------------------------------------------------------------------------
# 8) Prepare plotting data in Figure 2 format
# ------------------------------------------------------------------------------

anom_df <- as.data.frame(eff_anom_best) %>%
  mutate(
    dpp_lab = factor(
      group,
      levels = c("-1", "0", "1"),
      labels = c("Slow", "Intermediate", "Fast")
    ),
    bio4_lab = factor(
      facet,
      levels = c("-1", "0", "1"),
      labels = c("Low seasonality", "Intermediate seasonality", "High seasonality")
    )
  )

write_csv(anom_df, file.path(OUT_DIR, "figure2_anom_plotdata.csv"))

# ==============================================================================
# FIGURE 2
# ==============================================================================

MODEL_FILE <- file.path(OUT_DIR, "m_anom_best.rds")
DATA_FILE  <- file.path(OUT_DIR, "dat2_train_obs.rds")

FIG_DIR <- "figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_FIG <- file.path(FIG_DIR, "Figure2_timing_and_gate_strength_compact.png")

stopifnot(file.exists(MODEL_FILE))
stopifnot(file.exists(DATA_FILE))

m_anom_best <- readRDS(MODEL_FILE)
dat2 <- readRDS(DATA_FILE)

# ------------------------------------------------------------------------------
# 9) Panel A prediction data (hard-coded M4 path)
# ------------------------------------------------------------------------------

eff_anom <- ggpredict(
  m_anom_best,
  terms = c(
    "gddE_sc [-2:2 by=0.05]",
    "dpp_sc [-1,0,1]",
    "bio4_sc [-1,0,1]"
  )
)

anom_df <- as.data.frame(eff_anom) %>%
  mutate(
    dpp_lab = factor(
      group,
      levels = c("-1", "0", "1"),
      labels = c("Slow", "Intermediate", "Fast")
    ),
    bio4_lab = factor(
      facet,
      levels = c("-1", "0", "1"),
      labels = c("Low seasonality", "Intermediate seasonality", "High seasonality")
    )
  )

# ------------------------------------------------------------------------------
# 10) Panel B compact gate-strength summary (hard-coded M4 path)
# ------------------------------------------------------------------------------

pred_fixef <- function(model, newdata) {
  predict(model, newdata = newdata, re.form = NA, allow.new.levels = TRUE)
}

gddE_low  <- as.numeric(quantile(dat2$gddE_sc, probs = 0.25, na.rm = TRUE))
gddE_high <- as.numeric(quantile(dat2$gddE_sc, probs = 0.75, na.rm = TRUE))
bio4_slices <- c(-1, 0, 1)

year_levels <- levels(dat2$year_f)
ref_year <- if ("2019" %in% year_levels) "2019" else year_levels[1]

gate_grid <- expand.grid(
  dpp_sc  = seq(-2, 2, length.out = 200),
  bio4_sc = bio4_slices
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

gate_df <- gate_grid %>%
  mutate(
    gate_strength = pred_fixef(m_anom_best, gate_high) - pred_fixef(m_anom_best, gate_low),
    bio4_lab = factor(
      bio4_sc,
      levels = bio4_slices,
      labels = c("Low seasonality", "Intermediate", "High seasonality")
    )
  )

# ------------------------------------------------------------------------------
# 11) Colors
# ------------------------------------------------------------------------------

pal_dpp <- c(
  "Slow" = "#C53E35",
  "Intermediate" = "#4C78A8",
  "Fast" = "#54A24B"
)

pal_bio4 <- c(
  "Low seasonality" = "#C53E35",
  "Intermediate" = "#4C78A8",
  "High seasonality" = "#54A24B"
)

# ------------------------------------------------------------------------------
# 12) Shared theme
# ------------------------------------------------------------------------------

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
      legend.title = element_text(size = 10.5),
      legend.text = element_text(size = 9.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# ------------------------------------------------------------------------------
# 13) Panel A
# ------------------------------------------------------------------------------

p_obs <- ggplot(
  anom_df,
  aes(x = x, y = predicted, color = dpp_lab, fill = dpp_lab)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.12,
    linewidth = 0
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ bio4_lab, nrow = 1) +
  scale_color_manual(
    values = pal_dpp,
    name = "Photoperiod\nprogression",
    guide = guide_legend(order = 1)
  ) +
  scale_fill_manual(
    values = pal_dpp,
    guide = "none"
  ) +
  labs(
    title = "Phenological timing",
    x = "Early spring forcing",
    y = "Hopkins deviation"
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_fig2() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(8, 6, 6, 6),
    panel.spacing = unit(0.6, "lines")
  )

# ------------------------------------------------------------------------------
# 14) Panel B
# ------------------------------------------------------------------------------

p_gate_small <- ggplot(
  gate_df,
  aes(x = dpp_sc, y = gate_strength, color = bio4_lab)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(
    values = pal_bio4,
    name = "Seasonality",
    guide = guide_legend(order = 2)
  ) +
  labs(
    title = "Forcing sensitivity",
    x = "Photoperiod progression rate",
    y = "Change in Hopkins deviation\n(high minus low early forcing)"
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_fig2() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(8, 12, 6, 10)
  )

# ------------------------------------------------------------------------------
# 15) Combine on one row
# ------------------------------------------------------------------------------

fig2 <- (p_obs | p_gate_small) +
  plot_layout(
    widths = c(3.7, 1.8),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = unit(10, "pt"),
    legend.box.spacing = unit(10, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )

# ------------------------------------------------------------------------------
# 16) Save
# ------------------------------------------------------------------------------

ggsave(
  filename = OUT_FIG,
  plot = fig2,
  width = 13,
  height = 4.9,
  units = "in",
  dpi = 600,
  bg = "white"
)

fig2