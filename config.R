library(terra)
library(sf)

# Master Environmental Parameter Matrix
RISK_PARAMS <- list(
  
  # SPATIAL INFLUENCE PARAMETERS
  spatial = list(
    # NOISE (d0 = decay half-distance in meters)
    noise = list(
      road_density        = list(weight = +18.00, d0 = 75,  type = "source"),
      urban_density       = list(weight = +5.00, d0 = 300,  type = "source"),
      rail_density        = list(weight = +18.00, d0 = 150,  type = "source"),
      airport_fraction    = list(weight = +55.00, d0 = 2000, type = "source"),
      industrial_fraction = list(weight = +16.00, d0 = 400,  type = "source"),
      green_fraction      = list(weight = -5.00, d0 = 100,  type = "sink") # Attenuates noise
    ),
    
    # AIR POLLUTION (NO2 / PM2.5 / PM10)
    pollution = list(
      road_density        = list(weight = +20.00, d0 = 250,  type = "source"),
      industrial_fraction = list(weight = +35.00, d0 = 1000,  type = "source"),
      airport_fraction    = list(weight = +15.00, d0 = 500, type = "source"),
      urban_density       = list(weight = +5.00, d0 = 300,  type = "source"), # Street canyon trapping
      green_fraction      = list(weight = -3.40, d0 = 150,  type = "sink")  # Vegetation filtering
    ),
    
    # URBAN HEAT ISLAND (UHI)
    heat = list(
      urban_density       = list(weight = +4.50, d0 = 150,  type = "source"),
      industrial_fraction = list(weight = +3.50, d0 = 200,  type = "source"),
      road_density        = list(weight = +3.50, d0 = 80,  type = "source"),
      airport_fraction    = list(weight = +3.50, d0 = 400, type = "source"),
      green_fraction      = list(weight = -2.2, d0 = 300,  type = "sink"),  # Evapotranspirative cooling
      water_fraction      = list(weight = -1.5, d0 = 500,  type = "sink"),  # Thermal buffering
      elevation_m         = list(weight = -0.01, d0 = 0,    type = "linear") # Lapse rate (~0.65C per 100m)
    ),
    
    # UV EXPOSURE
    uv = list(
      elevation_m         = list(weight = +0.005, d0 = 0,   type = "linear"), # ~5% boost per 1000m
      urban_density       = list(weight = -0.15,  d0 = 50,  type = "sink"),   # Building shade
      green_fraction      = list(weight = -0.20,  d0 = 50,  type = "sink")    # Tree canopy shade
    )
  ),
  
  # TEMPORAL (24-HOUR) DIURNAL PROFILES (0.0 to 1.0)

  temporal = list(
    traffic     = c(0.1,0.1,0.1,0.1,0.1,0.15,0.25,0.75,1,0.85,0.65,0.5,0.5,0.5,0.5,0.65,0.75,1,0.85,0.65,0.5,0.35,0.25,0.15),
    industry    = c(0.1,0.1,0.1,0.1,0.1,0.15,0.2,0.35,0.65,0.85,1,1,1,1,1,1,1,0.85,0.65,0.35,0.2,0.1,0.1,0.1),
    airport     = c(0.25,0.25,0.25,0.25,0.25,0.25,0.45,0.85,1,1,1,1,1,1,1,1,1,1,1,1,1,0.85,0.65,0.45),
    uhi_lag     = c(0.6,0.55,0.5,0.5,0.5,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1,0.95,0.9,0.85,0.8,0.75,0.7,0.65),
    bg_noise_db = c(42, 40, 40, 40, 42, 45, 50, 55, 55, 53, 52, 52, 53, 54, 55, 56, 55, 53, 50, 48, 46, 45, 44, 43)
    
    # traffic    = c(0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 0.9, 1.0, 0.8, 0.6, 0.6, 0.7, 0.7, 0.7, 0.8, 0.9, 1.0, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1),
    # industry   = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5, 0.3, 0.2, 0.2, 0.2, 0.2),
    # solar_heat = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0, 0.95, 0.9, 0.8, 0.65, 0.45, 0.2, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0),
    # solar_uv   = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.6, 0.85, 0.98, 1.0, 0.95, 0.8, 0.55, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  ),
  
  # Physical effect scaling constants
  physical_effects = list(
    baseline_pollution     = 10,
    max_road_rail_noise_db = 1, # Max dBA added by peak traffic
    max_airport_noise_db   = 1, # Max dBA added near flight paths
    max_uhi_temp_c         = 1,  # Max °C added by urban heat island
    max_local_pollution_ug = 1  # Max µg/m³ added by traffic/industry
  )
)

# 1. Event Overlays Configuration
#Events defined in natural units
events_config <- list(
  # Cool zone subtracts 3.5 °C locally
  list(
    id = "municipal_cooling_centers",
    label = "Municipal Cooling Centers Network",
    risks = "heat",
    effect_type = "mitigation",
    coords = list(
      # Building 1: Baixa Library Polygon
      matrix(c(-92200, -105200, -91800, -105200, -91800, -104700, -92200, -104700), ncol = 2, byrow = TRUE),
      # Building 2: Secondary Library Point
      matrix(c(-93100, -104200, -92800, -104200, -92800, -103700, -93100, -103700), ncol = 2, byrow = TRUE),
      #c(-93100, -104200),
      # Building 3: Civic Center Point
      matrix(c(-90800, -106100, -90400, -106100, -90400, -105600, -90800, -105600), ncol = 2, byrow = TRUE)
      
      #c(-90800, -106100)
    ),
    magnitude = -3.5, # -3.5 °C cooling power at each site
    d0_m = 150,
    start_time = as.POSIXct("2026-09-28 00:00:00", tz = "Europe/Lisbon"),
    end_time   = as.POSIXct("2026-10-02 23:59:59", tz = "Europe/Lisbon"),
    active_hours = 10:18
  ),
  
  # Festival adds 30 dBA noise
  list(
    id = "lisbon_sound_fest",
    label = "Lisbon Sound Fest",
    risks = "noise",
    effect_type = "additive",
    coords = c(x = -91500, y = -106000),
    magnitude = 45.0, # +30 dBA
    d0_m = 1100,
    start_time = as.POSIXct("2026-09-28 18:00:00", tz = "UTC"),
    end_time   = as.POSIXct("2026-09-29 23:00:00", tz = "UTC"),
    active_hours = 18:23
  ),
  
  # Airbase fire adds PM2.5 (85 µg/m³) and Heat (4.0 °C)
  list(
    id = "montijo_fire_emergency",
    label = "Montijo Airbase Industrial Fire",
    risks = c("pollution", "heat"),
    effect_type = "additive",
    coords = matrix(c(-79000, -106500, -77000, -106500, -77000, -105500, -79000, -105500), ncol = 2, byrow = TRUE),
    magnitude = list(pollution = 85.0, heat = 14.0),
    d0_m = 1500,
    start_time = as.POSIXct("2026-09-29 06:00:00", tz = "UTC"),
    end_time   = as.POSIXct("2026-09-29 22:00:00", tz = "UTC"),
    active_hours = 0:23
  )
)
