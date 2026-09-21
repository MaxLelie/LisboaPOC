library(terra)
library(mgcv)
library(jsonlite)

# ==============================================================================
# 1. SETUP BASE SPATIAL GRID FOR LISBON
# ==============================================================================
cat("Setting up spatial grid for Lisbon...\n")

# Extent roughly covering Central Lisbon to Monsanto & Tagus riverfront
lisbon_extent <- ext(-9.24, -9.09, 38.69, 38.80)
# ~150 x 150 cell resolution for quick rendering
base_grid <- rast(lisbon_extent, res = 0.0025, crs = "EPSG:4326")

grid_coords <- as.data.frame(crds(base_grid))
names(grid_coords) <- c("lon", "lat")

# Key spatial landmarks for realistic GAM patterns
grid_coords$dist_river <- abs(grid_coords$lat - 38.69)
grid_coords$is_monsanto <- ifelse(grid_coords$lon > -9.21 & grid_coords$lon < -9.17 & 
                                    grid_coords$lat > 38.72 & grid_coords$lat < 38.75, 1, 0)
grid_coords$is_baixa    <- ifelse(grid_coords$lon > -9.15 & grid_coords$lon < -9.12 & 
                                    grid_coords$lat > 38.70 & grid_coords$lat < 38.73, 1, 0)

# ==============================================================================
# 2. FIT SPATIAL GAM MODELS FOR HAZARDS
# ==============================================================================
cat("Fitting Spatial GAMs for Heat, PM2.5, and Noise...\n")
set.seed(2026)

# A. Heat / UHI GAM Model
sensors_heat <- data.frame(
  lon = runif(80, -9.24, -9.09),
  lat = runif(80, 38.69, 38.80)
)
sensors_heat$is_monsanto <- ifelse(sensors_heat$lon > -9.21 & sensors_heat$lon < -9.17 & 
                                     sensors_heat$lat > 38.72 & sensors_heat$lat < 38.75, 1, 0)
sensors_heat$is_baixa    <- ifelse(sensors_heat$lon > -9.15 & sensors_heat$lon < -9.12 & 
                                     sensors_heat$lat > 38.70 & sensors_heat$lat < 38.73, 1, 0)
sensors_heat$uhi_val <- 3.5 * sensors_heat$is_baixa - 4.0 * sensors_heat$is_monsanto + rnorm(80, sd=0.4)

gam_heat <- gam(uhi_val ~ s(lon, lat, k = 25) + is_monsanto + is_baixa, data = sensors_heat)

# B. PM2.5 Air Pollution GAM Model
sensors_pm25 <- data.frame(
  lon = runif(60, -9.24, -9.09),
  lat = runif(60, 38.69, 38.80)
)
# Pollution concentrated near major transport corridors (Avenida da Liberdade / Baixa)
sensors_pm25$is_baixa <- ifelse(sensors_pm25$lon > -9.15 & sensors_pm25$lon < -9.12 & 
                                  sensors_pm25$lat > 38.70 & sensors_pm25$lat < 38.73, 1, 0)
sensors_pm25$pm_val <- 18 + 15 * sensors_pm25$is_baixa + rnorm(60, sd=2)

gam_pm25 <- gam(pm_val ~ s(lon, lat, k = 20) + is_baixa, data = sensors_pm25)

# C. Ambient Noise GAM Model
sensors_noise <- data.frame(
  lon = runif(60, -9.24, -9.09),
  lat = runif(60, 38.69, 38.80)
)
sensors_noise$dist_river <- abs(sensors_noise$lat - 38.69)
sensors_noise$noise_val  <- 45 + 20 * (1 - pmin(sensors_noise$dist_river * 10, 1)) + rnorm(60, sd=3)

gam_noise <- gam(noise_val ~ s(lon, lat, k = 20), data = sensors_noise)

# Predict base spatial anomalies across grid
grid_coords$base_heat  <- predict(gam_heat, newdata = grid_coords)
grid_coords$base_pm25  <- predict(gam_pm25, newdata = grid_coords)
grid_coords$base_noise <- predict(gam_noise, newdata = grid_coords)

# ==============================================================================
# 3. HELPER FUNCTION TO EXPORT DAILY MASTER JSON
# ==============================================================================
export_daily_master_json <- function(grid_df, base_raster, date_str, output_path) {
  
  # Normalize vector to 0.0 - 1.0 scale
  norm <- function(v) {
    v_clean <- ifelse(is.na(v), min(v, na.rm=T), v)
    m1 <- min(v_clean, na.rm=T)
    m2 <- max(v_clean, na.rm=T)
    if (m2 == m1) return(rep(0, length(v)))
    return((v_clean - m1) / (m2 - m1))
  }
  
  hours_list <- list()
  
  for (h in 0:23) {
    # Diurnal temperature cycle peaking at 15:00
    heat_factor  <- sin((h - 8) * pi / 12)
    # Traffic surge for PM2.5 and Noise at 08:00 and 18:00
    rush_factor  <- ifelse(h %in% c(7,8,9,17,18,19), 1.4, 0.7)
    night_factor <- ifelse(h >= 22 | h <= 5, 0.3, 1.0)
    
    # Calculate live values per hour
    h_heat  <- pmax(0, grid_df$base_heat * heat_factor)
    h_pm25  <- grid_df$base_pm25 * rush_factor
    h_noise <- grid_df$base_noise * rush_factor * night_factor
    
    hours_list[[as.character(h)]] <- list(
      heat  = round(as.vector(norm(h_heat)), 3),
      pm25  = round(as.vector(norm(h_pm25)), 3),
      noise = round(as.vector(norm(h_noise)), 3)
    )
  }
  
  master_structure <- list(
    date  = date_str,
    bbox  = c(ext(base_raster)[1], ext(base_raster)[3], ext(base_raster)[2], ext(base_raster)[4]),
    rows  = nrow(base_raster),
    cols  = ncol(base_raster),
    hours = hours_list
  )
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_json(master_structure, output_path, auto_unbox = TRUE)
  cat(sprintf("✓ Successfully generated: %s\n", output_path))
}

# ==============================================================================
# 4. EXECUTE GENERATION LOOP ACROSS DATES
# ==============================================================================
dates <- c("2026-10-26", "2026-10-27", "2026-10-28", "2026-10-29", "2026-10-30")

cat("\nGenerating daily files in www/data/...\n")
for (d in dates) {
  out_file <- sprintf("www/data/%s.json", d)
  export_daily_master_json(grid_coords, base_grid, d, out_file)
}

cat("\nAll daily master JSON files generated successfully!\n")