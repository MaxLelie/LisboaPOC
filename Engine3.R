library(terra)
library(sf)
library(osmdata)
library(httr2)
library(jsonlite)
library(dplyr)
library(elevatr)
library(gifski)

source("/home/maxosaurus/Desktop/LisboaPOC/functions.R")
source("/home/maxosaurus/Desktop/LisboaPOC/config.R")

# ==== 1. BASE GRID SETUP ====

lisbon_bbox <- c(xmin = -9.30, xmax = -9.00, ymin = 38.60, ymax = 38.85)
grid_250m <- build_spatial_grid(lisbon_bbox, crs_target = 3763, cell_size_m = 250)

# ==== 2. DATA ACQUISITION & RASTERIZATION ====

osm_layers     <- grab_all_osm_and_elevation(lisbon_bbox, grid_250m)
raw_grid_stack <- rasterize_osm_layers(grid_250m, osm_layers)

inspect_raw_osm(osm_layers)
inspect_raster_stack(raw_grid_stack, "Raw Feature")

# ==== 3. SPATIAL DECAY & COMPONENT STACK ====

# Un-collapsed spatial stack containing individual decayed and weighted layers
individual_effects_stack <- build_individual_effect_grids(raw_grid_stack, RISK_PARAMS)
plot(individual_effects_stack$noise__airport_fraction)
plot(individual_effects_stack$noise__road_density)
plot(individual_effects_stack$noise__urban_density)


# Overall normalized static risk layers for display/inspection
static_risk_stack <- build_combined_static_risks(individual_effects_stack, RISK_PARAMS)
plot(static_risk_stack, main = "Combined Static Risk Summaries")

# ==== 4. METEO FETCH & DYNAMIC TIMESTEP DIAGNOSTICS ====

dynamic_meteo <- fetch_meteo_data()

# Diagnostic timestamp matching an active fire / sound festival event
test_time <- as.POSIXct("2026-09-30 20:00:00", tz = "Europe/Lisbon")

test_time <- as.POSIXct("2026-09-29 16:00:00", tz = "Europe/Lisbon")

diag <- compute_dynamic_risk_at_time(
  target_time = test_time,
  spatial_effects_stack = individual_effects_stack,
  meteo_df = dynamic_meteo,
  events_config = events_config,
  params = RISK_PARAMS,
  return_intermediates = TRUE
)

cat("Active Event Overlay Metadata:\n")
print(diag$active_labels)

# par(mfrow = c(2, 2))
par(mfrow = c(1, 1))
plot(diag$additive_events$heat, main = "Heat Event Overlay (Fire)")
plot(diag$mitigation_events$heat, main = "Heat mitigation Overlay (cool zones)")

plot(diag$dynamic_stack$heat, main = "Combined Dynamic Heat Risk (°C)")
plot(diag$additive_events$noise, main = "Noise Event Overlay (Noise)")

plot(diag$dynamic_stack$noise, main = "Combined Dynamic Noise Risk (Dcb)")

plot(diag$additive_events$pollution, main = "Pollution Event Overlay (Fire PM2.5)")
plot(diag$dynamic_stack$pollution, main = "Combined Dynamic PM2.5 (µg/m³)")
par(mfrow = c(1, 1))

# ==== 5. FULL TIME-SERIES RUN & GIF RENDERING ====

time_series_risks <- generate_time_series_stack(
  spatial_effects_stack = individual_effects_stack,
  meteo_df = dynamic_meteo,
  events_config = events_config,
  params = RISK_PARAMS
)

render_risk_gif(time_series_risks, risk_name = "heat", output_gif = "heat_dynamic.gif", fps = 5)
render_risk_gif(time_series_risks, risk_name = "noise", output_gif = "noise_dynamic.gif", fps = 5)

# Optional JSON Export
export_daily_jsons(individual_effects_stack, dynamic_meteo, events_config, RISK_PARAMS)
