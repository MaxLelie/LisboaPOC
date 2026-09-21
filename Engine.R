library(terra)
library(jsonlite)

# ==============================================================================
# 1. BASE GRID SETUP FOR LISBON
# ==============================================================================
lisbon_ext <- ext(-9.24, -9.09, 38.69, 38.80)
base_grid  <- rast(lisbon_ext, res = 0.0025, crs = "EPSG:4326")
coords     <- as.data.frame(crds(base_grid))
names(coords) <- c("lon", "lat")

# Key locations for spatial effects
monsanto_center <- c(-9.185, 38.735) # Green space / cooling
baixa_center    <- c(-9.138, 38.712) # Urban Heat Island (UHI) high density
airport_center  <- c(-9.135, 38.775) # Humberto Delgado Airport
city_center     <- c(-9.140, 38.725) # City center noise hub
event_center    <- c(-9.120, 38.750) # Evening concert event (e.g., Parque das Nações)

dist_deg <- function(c1, c2) sqrt((coords$lon - c1[1])^2 + (coords$lat - c1[2])^2)

coords$d_monsanto <- dist_deg(monsanto_center, NULL)
coords$d_baixa    <- dist_deg(baixa_center, NULL)
coords$d_airport  <- dist_deg(airport_center, NULL)
coords$d_city     <- dist_deg(city_center, NULL)
coords$d_event    <- dist_deg(event_center, NULL)

# ==============================================================================
# 2. SYNTHETIC LAYER BUILDERS
# ==============================================================================

# A. TEMPERATURE LAYER
# Baseline cosine daily cycle between T_min and T_max + UHI/park anomalies
get_temperature <- function(hour, t_min = 16, t_max = 28) {
  # Diurnal cosine wave (minimum at 06:00, peak at 15:00)
  diurnal_temp <- t_min + (t_max - t_min) * 0.5 * (1 - cos((hour - 6) * pi / 12))
  
  # UHI spatial anomalies (in °C offset)
  uhi_heat <- 3.5 * exp(-((coords$d_baixa / 0.02)^2))      # Red increase (high density)
  uhi_cool <- -2.5 * exp(-((coords$d_monsanto / 0.025)^2))  # Green decrease (park)
  
  # Per-cell temperature
  return(diurnal_temp + uhi_heat + uhi_cool)
}

# B. NOISE LAYER
# Concentric circles around hubs + temporary evening event (18:00–23:00)
get_noise <- function(hour) {
  # Base background ambient noise (dB)
  base_noise <- 45
  
  # Concentric attenuation noise (decaying with distance)
  noise_airport <- 35 * exp(-coords$d_airport / 0.03)
  noise_city    <- 25 * exp(-coords$d_city / 0.02)
  
  # Temporary event (e.g., concert active between 18:00 and 23:00)
  event_active  <- ifelse(hour >= 18 & hour <= 23, 1, 0)
  noise_event   <- 40 * exp(-coords$d_event / 0.015) * event_active
  
  return(pmin(85, base_noise + noise_airport + noise_city + noise_event))
}

# C. POLLUTION (PM2.5) LAYER
# Low spatial frequency (smooth wide wave field) + small hourly variations
set.seed(42)
coords$pollution_macro <- 15 + 10 * sin(coords$lon * 80) * cos(coords$lat * 80)

get_pollution <- function(hour) {
  # Rush hour multipliers at 08:00 and 18:00
  rush_factor <- 1.0 + 0.4 * exp(-((hour - 8)^2) / 4) + 0.5 * exp(-((hour - 18)^2) / 4)
  
  # Small spatial jitter per hour
  micro_variation <- 2 * sin(coords$lon * 300 + hour)
  
  return(pmax(5, coords$pollution_macro * rush_factor + micro_variation))
}

# ==============================================================================
# 3. EXPORT DAILY MASTER JSON
# ==============================================================================
export_day_json <- function(date_str, output_path) {
  hours_list <- list()
  
  for (h in 0:23) {
    hours_list[[as.character(h)]] <- list(
      temp  = round(get_temperature(h), 1), # values in °C
      noise = round(get_noise(h), 1),        # values in dB
      pm25  = round(get_pollution(h), 1)     # values in µg/m³
    )
  }
  
  payload <- list(
    date  = date_str,
    bbox  = c(ext(base_grid)[1], ext(base_grid)[3], ext(base_grid)[2], ext(base_grid)[4]),
    rows  = nrow(base_grid),
    cols  = ncol(base_grid),
    hours = hours_list
  )
  
  dir.create(dirname(output_path), recursive = FALSE, showWarnings = FALSE)
  write_json(payload, output_path, auto_unbox = TRUE)
  cat(sprintf("✓ Exported: %s\n", output_path))
}

# Generate 5 test dates
dates <- c("2026-10-26", "2026-10-27", "2026-10-28", "2026-10-29", "2026-10-30")
for (d in dates) {
  export_day_json(d, sprintf("www/data/%s.json", d))
}

