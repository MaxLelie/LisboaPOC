

# ==== 1. GRID & DATA ACQUISITION ====

build_spatial_grid <- function(bbox, crs_target = 3763, cell_size_m = 250) {
  cat("Creating target grid at EPSG:", crs_target, "\n")
  
  box_polygon <- st_polygon(list(matrix(c(
    bbox["xmin"], bbox["ymin"],
    bbox["xmax"], bbox["ymin"],
    bbox["xmax"], bbox["ymax"],
    bbox["xmin"], bbox["ymax"],
    bbox["xmin"], bbox["ymin"]
  ), ncol = 2, byrow = TRUE)))
  
  area_sf <- st_sf(geometry = st_sfc(box_polygon, crs = 4326))
  area_projected <- st_transform(area_sf, crs_target)
  
  grid_geom <- st_make_grid(area_projected, cellsize = cell_size_m, square = TRUE)
  grid_sf <- st_sf(cell_id = seq_along(grid_geom), geometry = grid_geom)
  
  intersects_mask <- as.logical(st_intersects(grid_sf, area_projected, sparse = FALSE))
  grid_sf <- grid_sf[intersects_mask, ]
  
  cat("  Grid built with", nrow(grid_sf), "cells (", cell_size_m, "m resolution)\n")
  return(grid_sf)
}

grab_osm_with_retry <- function(query, max_attempts = 5, initial_wait = 3) {
  for (attempt in 1:max_attempts) {
    res <- tryCatch({
      osmdata_sf(query)
    }, error = function(e) {
      cat(sprintf("  [Attempt %d/%d] OSM download failed: %s. Retrying in %d seconds...\n", 
                  attempt, max_attempts, e$message, initial_wait * attempt))
      Sys.sleep(initial_wait * attempt)
      return(NULL)
    })
    if (!is.null(res)) return(res)
  }
  return(NULL)
}

grab_all_osm_and_elevation <- function(bbox_vector, grid_sf, cache_dir = "data/raw") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  fetch_or_load <- function(filename, query_fn) {
    filepath <- file.path(cache_dir, filename)
    if (file.exists(filepath)) {
      cat("Loading cached", filename, "...\n")
      return(readRDS(filepath))
    } else {
      cat("Fetching", filename, "...\n")
      data <- query_fn()
      saveRDS(data, filepath)
      return(data)
    }
  }
  
  roads <- fetch_or_load("osm_roads.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary", "tertiary"))
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    res$osm_lines
  })
  
  green <- fetch_or_load("osm_green.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_features(features = c('"leisure"="park"', '"landuse"="forest"', '"landuse"="grass"', '"natural"="wood"'))
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    c(st_geometry(res$osm_polygons), st_geometry(res$osm_multipolygons))
  })
  
  water <- fetch_or_load("osm_water.rds", function() {
    cat("  Loading vector coastline water via Natural Earth...\n")
    suppressMessages(sf_use_s2(FALSE))
    ocean <- rnaturalearth::ne_download(scale = 10, type = 'ocean', category = 'physical', returnclass = "sf")
    
    box_poly <- st_as_sfc(st_bbox(c(
      xmin = unname(bbox_vector["xmin"]), 
      ymin = unname(bbox_vector["ymin"]), 
      xmax = unname(bbox_vector["xmax"]), 
      ymax = unname(bbox_vector["ymax"])
    ), crs = 4326))
    
    water_clipped <- st_intersection(st_make_valid(ocean), box_poly)
    return(st_sf(geometry = st_geometry(water_clipped)))
  })
  
  urban_landuse <- fetch_or_load("osm_urban_landuse.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_feature(key = "landuse", value = c("commercial", "industrial", "residential", "retail"))
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    c(st_geometry(res$osm_polygons), st_geometry(res$osm_multipolygons))
  })
  
  industry <- fetch_or_load("osm_industry.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_feature(key = "landuse", value = "industrial")
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    c(st_geometry(res$osm_polygons), st_geometry(res$osm_multipolygons))
  })
  
  airports <- fetch_or_load("osm_airports.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_feature(key = "aeroway", value = c("aerodrome", "airport"))
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    geoms <- list()
    if (!is.null(res$osm_polygons) && nrow(res$osm_polygons) > 0) geoms[[length(geoms) + 1]] <- st_geometry(res$osm_polygons)
    if (!is.null(res$osm_multipolygons) && nrow(res$osm_multipolygons) > 0) geoms[[length(geoms) + 1]] <- st_geometry(res$osm_multipolygons)
    if (length(geoms) == 0) return(NULL)
    st_sf(geometry = do.call(c, geoms), crs = 4326)
  })
  
  railways <- fetch_or_load("osm_railways.rds", function() {
    q <- opq(bbox = bbox_vector) |>
      add_osm_feature(key = "railway", value = c("rail", "light_rail"))
    res <- grab_osm_with_retry(q)
    if (is.null(res)) return(NULL)
    res$osm_lines
  })
  
  elevation <- fetch_or_load("dem_elevation.rds", function() {
    cat("Fetching DEM Elevation Data via AWS...\n")
    elev_sp <- get_elev_raster(locations = grid_sf, z = 10, clip = "bbox")
    rast(elev_sp)
  })
  
  return(list(
    roads = roads, green = green, water = water, 
    urban_landuse = urban_landuse, industry = industry,
    airports = airports, railways = railways, elevation = elevation
  ))
}

inspect_raw_osm <- function(osm_list) {
  cat("Rendering fast vector preview of raw OSM features...\n")
  par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
  
  safe_plot <- function(sf_obj, title_str, col_str, add_mode = FALSE) {
    if (!is.null(sf_obj) && length(sf_obj) > 0) {
      v <- vect(sf_obj)
      if (length(v) > 0) {
        plot(v, main = if(!add_mode) title_str else "", col = col_str, border = NA, add = add_mode)
        return(TRUE)
      }
    }
    if (!add_mode) plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = paste(title_str, "(Empty)"))
    return(FALSE)
  }
  
  safe_plot(osm_list$roads, "Roads (Lines)", "grey30")
  safe_plot(osm_list$green, "Green Spaces (Polygons)", "forestgreen")
  safe_plot(osm_list$water, "Water Bodies (Polygons)", "dodgerblue")
  safe_plot(osm_list$airports, "Airports (Polygons)", "darkorange")
  safe_plot(osm_list$railways, "Railways (Lines)", "firebrick")
  safe_plot(osm_list$industry, "Industrial (Polygons)", "gold")
  safe_plot(osm_list$urban_landuse, "urban_landuse", "royalblue")
  
  plot(1, type = "n", xlim = c(-9.30, -9.00), ylim = c(38.60, 38.85), 
       main = "ALL LAYERS COMBINED", xlab = "", ylab = "", axes = FALSE)
  
  safe_plot(osm_list$water, "", "dodgerblue", add_mode = TRUE)
  safe_plot(osm_list$green, "", "forestgreen", add_mode = TRUE)
  safe_plot(osm_list$industry, "", "gold", add_mode = TRUE)
  safe_plot(osm_list$airports, "", "darkorange", add_mode = TRUE)
  safe_plot(osm_list$roads, "", "grey30", add_mode = TRUE)
  safe_plot(osm_list$railways, "", "firebrick", add_mode = TRUE)
  
  par(mfrow = c(1, 1))
}

# ==== 2. RASTERIZATION & SPATIAL DECAY ====

rasterize_osm_layers <- function(grid_sf, osm_list, target_crs = 3763, cell_res_m = 250) {
  cat("\n[STEP 1] Generating raw gridded features...\n")
  
  v_grid <- vect(grid_sf)
  r_template <- rast(v_grid, res = cell_res_m)
  cell_area_ha <- (cell_res_m * cell_res_m) / 10000
  
  r_id <- r_template
  values(r_id) <- 1:ncell(r_id)
  grid_polys <- st_as_sf(as.polygons(r_id, aggregate = FALSE))
  names(grid_polys)[1] <- "cell_id"
  values(r_template) <- 0
  
  safe_transform <- function(obj, target_type = "POLYGON") {
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "osmdata") || inherits(obj, "osmdata_sf")) {
      geoms <- list()
      if (target_type == "POLYGON") {
        if (!is.null(obj$osm_polygons) && nrow(obj$osm_polygons) > 0) geoms[[length(geoms) + 1]] <- st_geometry(obj$osm_polygons)
        if (!is.null(obj$osm_multipolygons) && nrow(obj$osm_multipolygons) > 0) geoms[[length(geoms) + 1]] <- st_geometry(obj$osm_multipolygons)
      } else if (target_type == "LINESTRING") {
        if (!is.null(obj$osm_lines) && nrow(obj$osm_lines) > 0) geoms[[length(geoms) + 1]] <- st_geometry(obj$osm_lines)
        if (!is.null(obj$osm_multilines) && nrow(obj$osm_multilines) > 0) geoms[[length(geoms) + 1]] <- st_geometry(obj$osm_multilines)
      }
      if (length(geoms) == 0) return(NULL)
      geom <- do.call(c, geoms)
    } else if (inherits(obj, "sf") || inherits(obj, "sfc")) {
      geom <- st_geometry(obj)
    } else {
      return(NULL)
    }
    
    if (length(geom) == 0) return(NULL)
    sf_out <- st_sf(geometry = geom) |> st_make_valid()
    geom_types <- unique(as.character(st_geometry_type(sf_out)))
    if (any(geom_types %in% c("GEOMETRYCOLLECTION", "GEOMETRY"))) {
      sf_out <- suppressWarnings(st_collection_extract(sf_out, target_type))
    }
    if (nrow(sf_out) == 0) return(NULL)
    return(st_transform(sf_out, target_crs))
  }
  
  # ==== UPDATED: Line Rasterization with Gamma Skewing ====
  rasterize_lines_clipped <- function(lines_sf, template_r, grid_polys, cell_area_ha, 
                                      max_cap_m_ha = 150, gamma = 0.35) {
    if (is.null(lines_sf) || nrow(lines_sf) == 0) {
      r_empty <- template_r
      values(r_empty) <- 0
      return(r_empty)
    }
    
    lines_union <- st_union(lines_sf) |> st_make_valid()
    lines_clean <- st_sf(geometry = st_cast(lines_union, "LINESTRING"))
    
    suppressWarnings({
      lines_clipped <- st_intersection(lines_clean, grid_polys)
    })
    
    if (nrow(lines_clipped) == 0) {
      r_empty <- template_r
      values(r_empty) <- 0
      return(r_empty)
    }
    
    lines_clipped$clip_len_m <- as.numeric(st_length(lines_clipped))
    cell_lengths <- aggregate(clip_len_m ~ cell_id, data = lines_clipped, FUN = sum)
    
    r_out <- template_r
    values(r_out) <- 0
    r_out[cell_lengths$cell_id] <- cell_lengths$clip_len_m
    
    # 1. Calculate linear density (m/ha)
    r_density <- r_out / cell_area_ha
    
    # 2. Clamp at saturation cap and normalize to [0, 1]
    r_norm <- clamp(r_density / max_cap_m_ha, upper = 1.0)
    
    # 3. Non-linear gamma curve to skew presence towards 1.0
    r_skewed <- r_norm^gamma
    
    return(r_skewed)
  }  
  green_p <- safe_transform(osm_list$green, "POLYGON")
  r_green <- if (!is.null(green_p) && nrow(green_p) > 0) rasterize(vect(green_p), r_template, cover = TRUE, background = 0) else r_template
  names(r_green) <- "green_fraction"
  
  water_p <- safe_transform(osm_list$water, "POLYGON")
  r_water <- if (!is.null(water_p) && nrow(water_p) > 0) rasterize(vect(water_p), r_template, cover = TRUE, background = 0) else r_template
  names(r_water) <- "water_fraction"
  
  airports_p <- safe_transform(osm_list$airports, "POLYGON")
  if (!is.null(airports_p) && nrow(airports_p) > 0) {
    airports_union <- st_union(airports_p) |> st_make_valid()
    r_airports <- rasterize(vect(airports_union), r_template, cover = TRUE, background = 0)
    r_airports <- clamp(r_airports, lower = 0, upper = 1)
  } else { 
    r_airports <- r_template 
  }
  names(r_airports) <- "airport_fraction"
  
  ind_p <- safe_transform(osm_list$industry, "POLYGON")
  if (!is.null(ind_p) && !is.null(airports_p) && nrow(ind_p) > 0 && nrow(airports_p) > 0) {
    suppressWarnings({ ind_p <- st_difference(ind_p, st_union(airports_p)) })
  }
  r_ind <- if (!is.null(ind_p) && nrow(ind_p) > 0) rasterize(vect(ind_p), r_template, cover = TRUE, background = 0) else r_template
  names(r_ind) <- "industrial_fraction"
  
  urban_p <- safe_transform(osm_list$urban_landuse, "POLYGON")
  r_urban <- if (!is.null(urban_p) && nrow(urban_p) > 0) rasterize(vect(urban_p), r_template, cover = TRUE, background = 0) else r_template
  names(r_urban) <- "urban_density"
  
  roads_l <- safe_transform(osm_list$roads, "LINESTRING")
  r_roads <- rasterize_lines_clipped(roads_l, r_template, grid_polys, cell_area_ha, max_cap_m_ha = 1200)
  names(r_roads) <- "road_density"
  
  railways_l <- safe_transform(osm_list$railways, "LINESTRING")
  r_railways <- rasterize_lines_clipped(railways_l, r_template, grid_polys, cell_area_ha, max_cap_m_ha = 800)
  names(r_railways) <- "rail_density"
  
  r_elev <- project(osm_list$elevation, r_template)
  names(r_elev) <- "elevation_m"
  
  return(c(r_green, r_water, r_urban, r_ind, r_roads, r_railways, r_airports, r_elev))
}

inspect_raster_stack <- function(r_stack, title_prefix = "Layer") {
  n <- nlyr(r_stack)
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  par(mfrow = c(rows, cols), mar = c(3, 3, 2, 4))
  for (i in 1:n) {
    plot(r_stack[[i]], main = paste(title_prefix, ":", names(r_stack)[i]), col = rev(terrain.colors(50)))
  }
  par(mfrow = c(1, 1))
}

# apply_spatial_decay <- function(rast_layer, d0_meters) {
#   if (is.null(d0_meters) || d0_meters <= 0) return(rast_layer)
#   
#   res_m <- res(rast_layer)[1]
#   radius_cells <- max(3, min(50, ceiling((d0_meters * 6.6) / res_m)))
#   window_dim <- radius_cells * 2 + 1
#   max_dist_m <- radius_cells * res_m
#   
#   x <- seq(-radius_cells, radius_cells) * res_m
#   grid_coords <- expand.grid(x = x, y = x)
#   dist_matrix <- matrix(sqrt(grid_coords$x^2 + grid_coords$y^2), 
#                         nrow = window_dim, ncol = window_dim)
#   
#   weights <- exp(-log(2) * dist_matrix / d0_meters)
#   weights[dist_matrix > max_dist_m] <- 0
#   
#   return(focal(rast_layer, w = weights, fun = "max", na.rm = TRUE))
# }

apply_spatial_decay <- function(rast_layer, d0_meters) {
  if (is.null(d0_meters) || d0_meters <= 0) return(rast_layer)
  
  res_m <- res(rast_layer)[1]
  
  # 1. Expand window radius so edge weights naturally drop to ~0.01 (approx 6.6 * d0)
  #    Capped at 50 cells to prevent extreme memory usage on large grids
  radius_cells <- max(3, min(50, ceiling((d0_meters * 6.6) / res_m)))
  window_dim <- radius_cells * 2 + 1
  max_dist_m <- radius_cells * res_m
  
  # 2. Build coordinate distance grid
  x <- seq(-radius_cells, radius_cells) * res_m
  grid_coords <- expand.grid(x = x, y = x)
  dist_matrix <- matrix(sqrt(grid_coords$x^2 + grid_coords$y^2), 
                        nrow = window_dim, ncol = window_dim)
  
  # 3. Calculate exponential decay
  weights <- exp(-log(2) * dist_matrix / d0_meters)
  
  # 4. CIRCULAR MASK: Zero out anything beyond the circular radius cutoff
  weights[dist_matrix > max_dist_m] <- 0
  
  # 5. Apply max focal window
  smoothed <- focal(rast_layer, w = weights, fun = "max", na.rm = TRUE)
  return(smoothed)
}

build_individual_effect_grids <- function(raw_stack, params) {
  effect_layers <- list()
  
  for (risk_name in names(params$spatial)) {
    risk_config <- params$spatial[[risk_name]]
    
    for (layer_name in names(risk_config)) {
      if (!layer_name %in% names(raw_stack)) next
      
      cfg <- risk_config[[layer_name]]
      raw_layer <- raw_stack[[layer_name]]
      
      if (cfg$type == "linear") {
        processed <- raw_layer * cfg$weight
      } else {
        smoothed <- apply_spatial_decay(raw_layer, cfg$d0)
        processed <- smoothed * cfg$weight
      }
      
      layer_key <- paste0(risk_name, "__", layer_name)
      effect_layers[[layer_key]] <- processed
    }
  }
  
  return(rast(effect_layers))
}

build_combined_static_risks <- function(indiv_stack, params) {
  combined_list <- list()
  
  for (risk_name in names(params$spatial)) {
    pattern <- paste0("^", risk_name, "__")
    matching_layers <- indiv_stack[[grep(pattern, names(indiv_stack))]]
    print(pattern)
    print(matching_layers)
    
    raw_sum <- sum(matching_layers, na.rm = TRUE)
    min_val <- global(raw_sum, "min", na.rm = TRUE)$min
    max_val <- global(raw_sum, "max", na.rm = TRUE)$max
    
    normalized_risk <- (raw_sum - min_val) / (max_val - min_val)
    names(normalized_risk) <- risk_name
    
    combined_list[[risk_name]] <- normalized_risk
  }
  
  return(rast(combined_list))
}

# ==== 3. METEO & EVENT PROCESSING ====

fetch_meteo_data <- function(lat = 38.7223, lon = -9.1393) {
  cat("\nFetching Open-Meteo Weather & Air Quality Data...\n")
  
  url_meteo <- sprintf(
    "https://api.open-meteo.com/v1/forecast?latitude=%.4f&longitude=%.4f&hourly=temperature_2m,relative_humidity_2m,shortwave_radiation,uv_index&forecast_days=7&timezone=Europe%%2FLisbon",
    lat, lon
  )
  resp_meteo <- req_perform(request(url_meteo))
  j_meteo <- resp_body_json(resp_meteo)$hourly
  
  url_aq <- sprintf(
    "https://air-quality-api.open-meteo.com/v1/air-quality?latitude=%.4f&longitude=%.4f&hourly=pm2_5,nitrogen_dioxide&forecast_days=7&timezone=Europe%%2FLisbon",
    lat, lon
  )
  resp_aq <- req_perform(request(url_aq))
  j_aq <- resp_body_json(resp_aq)$hourly
  
  clean_vec <- function(lst) {
    if (is.null(lst)) return(rep(NA_real_, length(j_meteo$time)))
    sapply(lst, function(x) if (is.null(x)) NA_real_ else as.numeric(x))
  }
  
  data.frame(
    datetime  = as.POSIXct(unlist(j_meteo$time), format = "%Y-%m-%dT%H:%M", tz = "Europe/Lisbon"),
    temp      = clean_vec(j_meteo$temperature_2m),
    solar_rad = clean_vec(j_meteo$shortwave_radiation),
    uv        = clean_vec(j_meteo$uv_index),
    pm25      = clean_vec(j_aq$pm2_5),
    no2       = clean_vec(j_aq$nitrogen_dioxide)
  )
}

# Helper function to extract spatial effect layers safely
get_layer_or_zero <- function(stk, layer_name, template_r) {
  if (!is.null(stk) && layer_name %in% names(stk)) {
    return(stk[[layer_name]])
  } else {
    return(template_r)
  }
}

generate_event_raster <- function(event_cfg, target_time, r_template, target_crs = 3763) {
  t_posix <- as.POSIXct(target_time, tz = "Europe/Lisbon")
  t_hour  <- as.POSIXlt(t_posix, tz = "Europe/Lisbon")$hour
  
  start_p <- as.POSIXct(event_cfg$start_time, tz = "Europe/Lisbon")
  end_p   <- as.POSIXct(event_cfg$end_time, tz = "Europe/Lisbon")
  
  is_within_date <- t_posix >= start_p && t_posix <= end_p
  is_within_hour <- t_hour %in% event_cfg$active_hours
  
  r_event <- r_template
  values(r_event) <- 0
  
  if (!is_within_date || !is_within_hour) return(r_event)
  
  coords <- event_cfg$coords
  
  # --- Flexible Geometry Builder (Points, Polygons, & Multi-location Lists) ---
  geoms <- list()
  
  if (is.list(coords) && !is.matrix(coords)) {
    # List of multiple locations (points or polygon matrices)
    for (item in coords) {
      if (is.matrix(item)) {
        if (!identical(item[1, ], item[nrow(item), ])) item <- rbind(item, item[1, ])
        geoms[[length(geoms) + 1]] <- st_polygon(list(item))
      } else {
        geoms[[length(geoms) + 1]] <- st_point(as.numeric(item))
      }
    }
  } else if (is.matrix(coords)) {
    # Single Polygon
    if (!identical(coords[1, ], coords[nrow(coords), ])) coords <- rbind(coords, coords[1, ])
    geoms[[1]] <- st_polygon(list(coords))
  } else {
    # Single Point
    geoms[[1]] <- st_point(as.numeric(coords))
  }
  
  geom_sf <- st_sf(geometry = st_sfc(geoms, crs = target_crs))
  
  # Rasterize footprint & apply focal distance decay
  r_src <- rasterize(vect(geom_sf), r_template, field = 1, background = 0)
  decay_r <- apply_spatial_decay(r_src, d0_meters = event_cfg$d0_m)
  
  return(decay_r)
}

# ==== 4. DYNAMIC ENVIRONMENTAL RISK ENGINE ====

compute_dynamic_risk_at_time <- function(target_time, spatial_effects_stack, meteo_df, 
                                         events_config, params, 
                                         return_intermediates = FALSE) {
  
  t_posix <- as.POSIXct(target_time, tz = "Europe/Lisbon")
  t_str   <- format(t_posix, "%Y-%m-%d %H:%M:%S")
  
  meteo_df$dt_str <- format(as.POSIXct(meteo_df$datetime, tz = "Europe/Lisbon"), "%Y-%m-%d %H:%M:%S")
  m_row <- meteo_df[meteo_df$dt_str == t_str, ]
  
  if (nrow(m_row) == 0 || any(is.na(m_row[c("temp", "uv", "pm25", "no2")]))) {
    return(NULL)
  }
  m_row <- m_row[1, ]
  
  t_hour <- as.POSIXlt(t_posix, tz = "Europe/Lisbon")$hour
  h_idx  <- t_hour + 1
  
  mult_traffic  <- params$temporal$traffic[h_idx]
  mult_industry <- params$temporal$industry[h_idx]
  mult_airport  <- params$temporal$airport[h_idx]
  mult_uhi      <- params$temporal$uhi_lag[h_idx]
  bg_noise_db   <- params$temporal$bg_noise_db[h_idx]
  
  r_template <- spatial_effects_stack[[1]]
  values(r_template) <- 0
  
  add_effects <- list(noise = r_template, pollution = r_template, heat = r_template)
  mit_effects <- list(noise = r_template, pollution = r_template, heat = r_template)
  active_metadata <- list()
  
  # Process Active Events
  for (ev in events_config) {
    decay_surface <- generate_event_raster(ev, t_posix, r_template)
    
    if (global(decay_surface, "max", na.rm = TRUE)[1, 1] > 0) {
      active_metadata[[ev$id]] <- list(
        label = ev$label, 
        type = ev$effect_type,
        coords = ev$coords, 
        risks = ev$risks
      )
      
      for (rk in ev$risks) {
        mag <- 0
        if (is.list(ev$magnitude) && !is.null(ev$magnitude[[rk]])) {
          mag <- as.numeric(ev$magnitude[[rk]])
        } else if (is.numeric(ev$magnitude) && !is.null(names(ev$magnitude)) && rk %in% names(ev$magnitude)) {
          mag <- as.numeric(ev$magnitude[rk])
        } else if (is.numeric(ev$magnitude) && length(ev$magnitude) == 1) {
          mag <- as.numeric(ev$magnitude)
        }
        
        effect_val <- decay_surface * mag
        
        if (ev$effect_type == "mitigation") {
          mit_effects[[rk]] <- mit_effects[[rk]] + effect_val
        } else {
          add_effects[[rk]] <- add_effects[[rk]] + effect_val
        }
      }
    }
  }
  
  # -------------------------------------------------------------
  # A. NOISE (dBA)
  # -------------------------------------------------------------
  n_road  <- get_layer_or_zero(spatial_effects_stack, "noise__road_density", r_template)
  n_rail  <- get_layer_or_zero(spatial_effects_stack, "noise__rail_density", r_template)
  n_air   <- get_layer_or_zero(spatial_effects_stack, "noise__airport_fraction", r_template)
  n_ind   <- get_layer_or_zero(spatial_effects_stack, "noise__industrial_fraction", r_template)
  n_urb   <- get_layer_or_zero(spatial_effects_stack, "noise__urban_density", r_template)
  n_green <- get_layer_or_zero(spatial_effects_stack, "noise__green_fraction", r_template)
  
  spatial_noise_delta <- (n_road + n_rail) * mult_traffic + 
    (n_air * mult_airport) + 
    (n_ind * mult_industry) + 
    n_urb + n_green
  
  noise_dBA <- bg_noise_db + spatial_noise_delta + add_effects$noise - mit_effects$noise
  noise_dBA <- clamp(noise_dBA, lower = 30)
  
  # -------------------------------------------------------------
  # B. AIR POLLUTION (PM2.5 in µg/m³)
  # -------------------------------------------------------------
  p_road  <- get_layer_or_zero(spatial_effects_stack, "pollution__road_density", r_template)
  p_ind   <- get_layer_or_zero(spatial_effects_stack, "pollution__industrial_fraction", r_template)
  p_air   <- get_layer_or_zero(spatial_effects_stack, "pollution__airport_fraction", r_template)
  p_urb   <- get_layer_or_zero(spatial_effects_stack, "pollution__urban_density", r_template)
  p_green <- get_layer_or_zero(spatial_effects_stack, "pollution__green_fraction", r_template)
  
  spatial_pol_delta <- (p_road * mult_traffic) + 
    (p_ind * mult_industry) + 
    (p_air * mult_airport) + 
    p_urb + p_green
  
  pollution_pm25 <- m_row$pm25 + spatial_pol_delta + add_effects$pollution - mit_effects$pollution
  pollution_pm25 <- clamp(pollution_pm25, lower = 0)
  
  # -------------------------------------------------------------
  # C. URBAN HEAT ISLAND (°C)
  # -------------------------------------------------------------
  h_urb   <- get_layer_or_zero(spatial_effects_stack, "heat__urban_density", r_template)
  h_ind   <- get_layer_or_zero(spatial_effects_stack, "heat__industrial_fraction", r_template)
  h_road  <- get_layer_or_zero(spatial_effects_stack, "heat__road_density", r_template)
  h_air   <- get_layer_or_zero(spatial_effects_stack, "heat__airport_fraction", r_template)
  h_green <- get_layer_or_zero(spatial_effects_stack, "heat__green_fraction", r_template)
  h_water <- get_layer_or_zero(spatial_effects_stack, "heat__water_fraction", r_template)
  h_elev  <- get_layer_or_zero(spatial_effects_stack, "heat__elevation_m", r_template)
  
  spatial_heat_delta <- (h_urb + h_ind + h_road + h_air) * mult_uhi + h_green + h_water + h_elev
  
  heat_c <- m_row$temp + spatial_heat_delta + add_effects$heat - mit_effects$heat
  
  # -------------------------------------------------------------
  # D. UV EXPOSURE INDEX
  # -------------------------------------------------------------
  u_elev  <- get_layer_or_zero(spatial_effects_stack, "uv__elevation_m", r_template)
  u_urb   <- get_layer_or_zero(spatial_effects_stack, "uv__urban_density", r_template)
  u_green <- get_layer_or_zero(spatial_effects_stack, "uv__green_fraction", r_template)
  
  uv_factor <- clamp(1.0 + u_elev + u_urb + u_green, lower = 0.1, upper = 2.0)
  uv_index  <- m_row$uv * uv_factor
  
  dynamic_stack <- c(noise_dBA, pollution_pm25, heat_c, uv_index)
  names(dynamic_stack) <- c("noise", "pollution", "heat", "uv")
  
  if (return_intermediates) {
    return(list(
      datetime = t_posix,
      dynamic_stack = dynamic_stack,
      meteo_raw = m_row,
      multipliers = c(traffic = mult_traffic, industry = mult_industry, airport = mult_airport, uhi = mult_uhi),
      additive_events = add_effects,
      mitigation_events = mit_effects,
      active_labels = active_metadata
    ))
  } else {
    return(dynamic_stack)
  }
}

generate_time_series_stack <- function(spatial_effects_stack, meteo_df, events_config, params) {
  valid_meteo <- meteo_df[!is.na(meteo_df$pm25) & !is.na(meteo_df$no2), ]
  time_steps  <- valid_meteo$datetime
  cat(sprintf("Processing %d hourly dynamic risk stacks...\n", length(time_steps)))
  
  results <- list()
  for (i in seq_along(time_steps)) {
    t_val <- time_steps[i]
    stk <- compute_dynamic_risk_at_time(t_val, spatial_effects_stack, valid_meteo, events_config, params)
    results[[format(as.POSIXct(t_val, tz="Europe/Lisbon"), "%Y-%m-%d %H:%M:%S")]] <- stk
  }
  
  return(results)
}

# ==== 5. EXPORTS & GIF VISUALIZATION ====

render_risk_gif <- function(time_series, risk_name = "heat", output_gif = "heat_risk.gif", fps = 4) {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  
  img_files <- c()
  timestamps <- names(time_series)
  
  # Auto-detect global range across time series to fix blank frame issue
  global_min <- Inf
  global_max <- -Inf
  
  for (t_str in timestamps) {
    stk <- time_series[[t_str]]
    if (!is.null(stk) && risk_name %in% names(stk)) {
      val_min <- global(stk[[risk_name]], "min", na.rm = TRUE)[1, 1]
      val_max <- global(stk[[risk_name]], "max", na.rm = TRUE)[1, 1]
      if (!is.na(val_min) && val_min < global_min) global_min <- val_min
      if (!is.na(val_max) && val_max > global_max) global_max <- val_max
    }
  }
  
  if (is.infinite(global_min) || is.infinite(global_max) || global_min == global_max) {
    global_min <- 0
    global_max <- 1
  }
  
  unit_str <- switch(risk_name, "heat" = "°C", "noise" = "dBA", "pollution" = "µg/m³", "uv" = "Index", "")
  palette_name <- switch(risk_name, "heat" = "YlOrRd", "noise" = "Inferno", "pollution" = "YlGnBu", "uv" = "Viridis", "YlOrRd")
  
  cat(sprintf("Rendering %d frames for '%s' [%.1f to %.1f %s]...\n", 
              length(timestamps), risk_name, global_min, global_max, unit_str))
  
  for (i in seq_along(timestamps)) {
    t_str <- timestamps[i]
    stk <- time_series[[t_str]]
    if (is.null(stk) || !risk_name %in% names(stk)) next
    
    r_layer <- stk[[risk_name]]
    file_path <- file.path(tmp_dir, sprintf("frame_%04d.png", i))
    
    png(file_path, width = 800, height = 800, res = 100)
    plot(r_layer, 
         main = sprintf("%s (%s) — %s", toupper(risk_name), unit_str, t_str),
         range = c(global_min, global_max), 
         col = hcl.colors(100, palette_name, rev = FALSE))
    dev.off()
    
    img_files <- c(img_files, file_path)
  }
  
  gifski(img_files, gif_file = output_gif, width = 800, height = 800, delay = 1 / fps)
  unlink(tmp_dir, recursive = TRUE)
  cat("Successfully saved GIF to:", output_gif, "\n")
}

export_daily_jsons <- function(spatial_effects_stack, meteo_df, events_config, params, output_dir = "json_outputs") {
  dir.create(output_dir, showWarnings = FALSE)
  
  valid_meteo <- meteo_df[!is.na(meteo_df$pm25) & !is.na(meteo_df$no2), ]
  valid_meteo$date_str <- as.character(as.Date(valid_meteo$datetime))
  
  for (d in unique(valid_meteo$date_str)) {
    day_rows <- valid_meteo[valid_meteo$date_str == d, ]
    hourly_payload <- list()
    
    cat(sprintf("Exporting daily JSON for date: %s (%d hours)...\n", d, nrow(day_rows)))
    
    for (i in 1:nrow(day_rows)) {
      t_val <- day_rows$datetime[i]
      time_key <- format(t_val, "%H:%M")
      
      diag <- compute_dynamic_risk_at_time(
        target_time = t_val,
        spatial_effects_stack = spatial_effects_stack,
        meteo_df = valid_meteo,
        events_config = events_config,
        params = params,
        return_intermediates = TRUE
      )
      
      if (is.null(diag)) next
      
      hourly_payload[[time_key]] <- list(
        timestamp = format(as.POSIXct(t_val, tz = "Europe/Lisbon"), "%Y-%m-%dT%H:%M:%SZ"),
        meteo = list(
          temp = diag$meteo_raw$temp,
          solar_rad = diag$meteo_raw$solar_rad,
          uv = diag$meteo_raw$uv,
          pm25 = diag$meteo_raw$pm25,
          no2 = diag$meteo_raw$no2
        ),
        multipliers = as.list(diag$multipliers),
        active_events = diag$active_labels,
        spatial_summary = list(
          mean_noise     = round(global(diag$dynamic_stack$noise, "mean", na.rm = TRUE)[1, 1], 4),
          mean_pollution = round(global(diag$dynamic_stack$pollution, "mean", na.rm = TRUE)[1, 1], 4),
          mean_heat      = round(global(diag$dynamic_stack$heat, "mean", na.rm = TRUE)[1, 1], 4),
          mean_uv        = round(global(diag$dynamic_stack$uv, "mean", na.rm = TRUE)[1, 1], 4)
        )
      )
    }
    
    daily_structure <- list(date = d, total_hours = length(hourly_payload), hours = hourly_payload)
    out_file <- file.path(output_dir, paste0("risk_data_", d, ".json"))
    write_json(daily_structure, path = out_file, auto_unbox = TRUE, pretty = TRUE)
  }
  
  cat("All daily JSON files exported successfully to folder:", output_dir, "\n")
}