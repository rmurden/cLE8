# =============================================================================
# point_estimate_maps.R
# 
# cle8_mapping_txgaca_within_state.R
#
# Choropleth and LISA Maps — Within-State cLE8 Scores
# States: California (CA), Texas (TX), Georgia (GA),
#
#  
#  
# Maps produced:
#   01_cle8_point_score_ws.png     -- continuous score choropleth
#   02_cle8_cat_point_ws.png       -- CVH category (Low/Moderate/High)
#   03_lisa_point_score_ws.png     -- LISA clusters
#
# Input:  cle8_three_scores_txgaca_within_state.csv
# Output: figures_maps_within_state/
# =============================================================================

library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(patchwork)
library(scales)
library(RColorBrewer)
library(viridis)

options(tigris_use_cache = TRUE)

dir.create("figures_maps_within_state", showWarnings = FALSE)

# =============================================================================
# SECTION 1: PALETTE SELECTION
# =============================================================================

# score_palette <- RColorBrewer::brewer.pal(11, "RdYlGn")
score_palette <- RColorBrewer::brewer.pal(11, "RdBu")
# score_palette <- viridis::viridis(11)
# score_palette <- RColorBrewer::brewer.pal(9, "Blues")

# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

cat("Loading within-state scored data...\n")

scores <- read_csv(
  "cle8_three_scores_txgaca_within_state.csv",
  show_col_types = FALSE
) |>
  mutate(GEOID = str_pad(as.character(CountyFIPS), width = 5, pad = "0"))

cat(sprintf("  %d counties loaded\n", nrow(scores)))

state_labels <- c(
  CA = "California\n(Outside Stroke Belt)",
  TX = "Texas\n(Partial Stroke Belt)",
  GA = "Georgia\n(Full Stroke Belt)"
)

# =============================================================================
# SECTION 3: LOAD GEOMETRIES
# =============================================================================

cat("Loading county geometries...\n")

counties_ca <- counties(state = "CA", cb = TRUE, resolution = "5m", year = 2022) |>
  st_transform(crs = 4326)
counties_tx <- counties(state = "TX", cb = TRUE, resolution = "5m", year = 2022) |>
  st_transform(crs = 4326)
counties_ga <- counties(state = "GA", cb = TRUE, resolution = "5m", year = 2022) |>
  st_transform(crs = 4326)

join_scores <- function(counties_sf, state_abbr) {
  counties_sf |>
    left_join(scores |> filter(StateAbbr == state_abbr), by = "GEOID")
}

geo_ca <- join_scores(counties_ca, "CA")
geo_tx <- join_scores(counties_tx, "TX")
geo_ga <- join_scores(counties_ga, "GA")

cat("  Geometries joined\n")

# =============================================================================
# SECTION 4: SHARED THEME AND ANNOTATION
# =============================================================================

map_theme <- theme_void(base_size = 10) +
  theme(
    plot.title        = element_text(size = 9, face = "bold", hjust = 0.5,
                                     margin = margin(b = 3)),
    plot.subtitle     = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.width  = unit(1.0, "cm"),
    legend.key.height = unit(0.3, "cm"),
    plot.margin       = margin(4, 4, 4, 4)
  )

fig_annotation_theme <- theme(
  plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 8,  hjust = 0.5, color = "grey40")
)

# =============================================================================
# SECTION 5: CONTINUOUS CHOROPLETH -- POINT ESTIMATE ONLY
# =============================================================================

cat("\nGenerating continuous choropleth map (point estimate)...\n")

make_ws_panel <- function(geo_sf, score_var, state_abbr, panel_title) {

  state_vals <- geo_sf[[score_var]]
  val_min    <- floor(min(state_vals, na.rm = TRUE) / 10) * 10
  val_max    <- ceiling(max(state_vals, na.rm = TRUE) / 10) * 10

  fill_sc <- scale_fill_gradientn(
    colours = score_palette,
    limits  = c(val_min, val_max),
    breaks  = c(val_min, val_max),
    labels  = c(paste0(val_min, "\n(Low CVH)"),
                paste0(val_max, "\n(High CVH)")),
    name    = paste0(state_abbr, " Score"),
    guide   = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(3, "cm"),
      barheight      = unit(0.3, "cm")
    )
  )

  ggplot(geo_sf) +
    geom_sf(aes(fill = .data[[score_var]]), color = "white", linewidth = 0.15) +
    fill_sc +
    labs(title = panel_title) +
    map_theme +
    theme(legend.position = "bottom")
}

p_ca_score <- make_ws_panel(geo_ca, "cle8_point", "CA", state_labels["CA"])
p_tx_score <- make_ws_panel(geo_tx, "cle8_point", "TX", state_labels["TX"])
p_ga_score <- make_ws_panel(geo_ga, "cle8_point", "GA", state_labels["GA"])

p_score <- (p_ca_score | p_tx_score | p_ga_score) +
  plot_layout(widths = c(1, 1.6, 1)) +
  plot_annotation(title = "cLE8 Score \u2014 Point Estimate",
                  theme = fig_annotation_theme)

ggsave(file.path("figures_maps_within_state", "01_cle8_point_score_ws.png"),
       p_score, width = 14, height = 7, dpi = 300)
cat("  Saved: 01_cle8_point_score_ws.png\n")

# =============================================================================
# SECTION 6: CVH CATEGORY MAP -- POINT ESTIMATE ONLY
# =============================================================================

cat("\nGenerating CVH category map (point estimate)...\n")

cvh_cat_colors <- c(
  "Low"      = "#d7191c",
  "Moderate" = "#fdae61",
  "High"     = "#1a9641"
)

cvh_cat_fill <- scale_fill_manual(
  values   = cvh_cat_colors,
  na.value = "grey85",
  name     = "CVH Category",
  limits   = c("Low", "Moderate", "High"),
  guide    = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
)

make_cat_panel <- function(geo_sf, cat_var, state_abbr) {
  geo_sf <- geo_sf |>
    mutate(cvh_cat = factor(.data[[cat_var]], levels = c("Low", "Moderate", "High")))

  ggplot(geo_sf) +
    geom_sf(aes(fill = cvh_cat), color = "white", linewidth = 0.15) +
    cvh_cat_fill +
    labs(title = state_labels[state_abbr]) +
    map_theme +
    theme(legend.position = "none")
}

p_ca_cat <- make_cat_panel(geo_ca, "cvh_cat_point", "CA")
p_tx_cat <- make_cat_panel(geo_tx, "cvh_cat_point", "TX")
p_ga_cat <- make_cat_panel(geo_ga, "cvh_cat_point", "GA")

p_cat <- (p_ca_cat | p_tx_cat | p_ga_cat) +
  plot_layout(widths = c(1, 1.6, 1), guides = "collect") &
  cvh_cat_fill & map_theme & theme(legend.position = "bottom")

p_cat <- p_cat +
  plot_annotation(
    title    = "CVH Category \u2014 Point Estimate Score",
    subtitle = "Low: 0\u201349 | Moderate: 50\u201379 | High: 80\u2013100",
    theme    = fig_annotation_theme
  )

ggsave(file.path("figures_maps_within_state", "02_cle8_cat_point_ws.png"),
       p_cat, width = 14, height = 7, dpi = 300)
cat("  Saved: 02_cle8_cat_point_ws.png\n")

# =============================================================================
# SECTION 7: LISA MAP -- POINT ESTIMATE ONLY
# =============================================================================

cat("\nComputing LISA statistics (point estimate)...\n")

lisa_colors <- c(
  "High-High"       = "#d7191c",
  "Low-Low"         = "#2c7bb6",
  "High-Low"        = "#fdae61",
  "Low-High"        = "#abd9e9",
  "Not Significant" = "grey85"
)

lisa_fill <- scale_fill_manual(
  values   = lisa_colors,
  na.value = "grey85",
  name     = "LISA Cluster",
  limits   = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"),
  guide    = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
)

compute_global_morans <- function(geo_sf, score_var) {
  x   <- geo_sf[[score_var]]
  nb  <- poly2nb(geo_sf, queen = TRUE)
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  mt  <- moran.test(x, lw, randomisation = TRUE, zero.policy = TRUE)
  p_fmt <- if (mt$p.value < 0.001) "< 0.001" else sprintf("= %.3f", mt$p.value)
  sprintf("Moran's I = %.3f (p %s)", mt$estimate["Moran I statistic"], p_fmt)
}

compute_lisa <- function(geo_sf, score_var, n_perm = 999, seed = 20240708) {
  x  <- geo_sf[[score_var]]
  nb <- poly2nb(geo_sf, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

  set.seed(seed)

  lm <- tryCatch(
    localmoran_perm(x, lw, nsim = n_perm, zero.policy = TRUE),
    error = function(e) {
      message("localmoran_perm unavailable; using analytic localmoran")
      localmoran(x, lw, zero.policy = TRUE)
    }
  )

  x_std   <- scale(x)[, 1]
  lag_std <- lag.listw(lw, x_std, zero.policy = TRUE)

  p_col <- if ("Pr(z != E(Ii))" %in% colnames(lm)) {
    "Pr(z != E(Ii))"
  } else if ("Pr(z != E(Ii)) Sim" %in% colnames(lm)) {
    "Pr(z != E(Ii)) Sim"
  } else {
    colnames(lm)[grep("Pr", colnames(lm))[1]]
  }

  p_val <- lm[, p_col]

  cluster <- case_when(
    p_val >= 0.05            ~ "Not Significant",
    x_std > 0 & lag_std > 0 ~ "High-High",
    x_std < 0 & lag_std < 0 ~ "Low-Low",
    x_std > 0 & lag_std < 0 ~ "High-Low",
    x_std < 0 & lag_std > 0 ~ "Low-High",
    TRUE                     ~ "Not Significant"
  )

  geo_sf |>
    mutate(
      lisa_cluster = factor(cluster,
                            levels = c("High-High", "Low-Low",
                                       "High-Low", "Low-High",
                                       "Not Significant")),
      local_i = lm[, "Ii"],
      lisa_p  = p_val
    )
}

make_lisa_panel <- function(geo_sf, state_abbr, morans_label) {
  ggplot(geo_sf) +
    geom_sf(aes(fill = lisa_cluster), color = "white", linewidth = 0.15) +
    lisa_fill +
    labs(title = state_labels[state_abbr], subtitle = morans_label) +
    map_theme +
    theme(legend.position = "none")
}

## Compute
mi_ca <- compute_global_morans(geo_ca, "cle8_point")
mi_tx <- compute_global_morans(geo_tx, "cle8_point")
mi_ga <- compute_global_morans(geo_ga, "cle8_point")

geo_ca_lisa <- compute_lisa(geo_ca, "cle8_point")
geo_tx_lisa <- compute_lisa(geo_tx, "cle8_point")
geo_ga_lisa <- compute_lisa(geo_ga, "cle8_point")

p_ca_lisa <- make_lisa_panel(geo_ca_lisa, "CA", mi_ca)
p_tx_lisa <- make_lisa_panel(geo_tx_lisa, "TX", mi_tx) +
  theme(legend.position = "bottom")
p_ga_lisa <- make_lisa_panel(geo_ga_lisa, "GA", mi_ga)

p_lisa <- (p_ca_lisa | p_tx_lisa | p_ga_lisa) / guide_area() +
  plot_layout(widths = c(1, 1.6, 1), heights = c(6, 1), guides = "collect")

p_lisa <- p_lisa +
  plot_annotation(
    title    = "LISA Clusters \u2014 cLE8 Point Estimate",
    theme    = fig_annotation_theme
  )

ggsave(file.path("figures_maps_within_state", "03_lisa_point_score_ws.png"),
       p_lisa, width = 14, height = 6, dpi = 300)
cat("  Saved: 03_lisa_point_score_ws.png\n")

## Cluster counts
cat("  Cluster counts:\n")
for (state in c("CA", "TX", "GA")) {
  geo_lisa <- get(paste0("geo_", tolower(state), "_lisa"))
  counts   <- table(geo_lisa$lisa_cluster)
  cat(sprintf("    %s: %s\n", state,
              paste(names(counts), counts, sep = "=", collapse = "  ")))
}

cat("\nDone. Three maps saved to figures_maps_within_state/\n")
