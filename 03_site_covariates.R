# 03_site_covariates.R

# Purpose
# Build the site-level covariates used in the final spatial
# and temporal analyses.

# Covariates include:
# habitat
# elevation
# sheep capture rate per 100 trap nights
# marked trail length within camera buffers

# The final recreational predictor is total marked-trail length
# within 500 m of each camera. Nearest trail distance is kept
# as an intermediate variable from the exploratory workflow but is
# not used in the final candidate model set.

# Outputs
# trail_500
# site_covs
# site_covs_umf

library(sf)
library(tidyverse)
library(elevatr)

# Trail layer and nearest trail covariate

trails_disturb <- st_read("data_raw/Trail_disturbance_correct/Trail_disturbance_correct.shp")

names(trails_disturb)
table(trails_disturb$Disturb)
st_crs(trails_disturb)

camera_sites <- master_presence %>%
  select(site_id, longitude, latitude) %>%
  distinct(site_id, .keep_all = TRUE) %>%
  filter(!is.na(longitude), !is.na(latitude))

camera_sites_sf <- camera_sites %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE)

camera_sites_sf <- st_transform(
  camera_sites_sf,
  st_crs(trails_disturb))

nearest_trail_index <- st_nearest_feature(
  camera_sites_sf,
  trails_disturb)

trail_covariates <- camera_sites_sf %>%
  mutate(
    nearest_trail_name = trails_disturb$NAME[nearest_trail_index],
    trail_disturbance = trails_disturb$Disturb[nearest_trail_index],
    distance_to_trail_m = as.numeric(
      st_distance(
        camera_sites_sf,
        trails_disturb[nearest_trail_index, ],
        by_element = TRUE))) %>%
  st_drop_geometry() %>%
  select(site_id, nearest_trail_name, trail_disturbance, distance_to_trail_m)

# Standardize habitat and feature type text

# Convert values such as "Coniferous_forest" to "Coniferous forest".
master_presence <- master_presence %>%
  mutate(
    habitat = gsub("_", " ", habitat) %>% trimws(),
    feature_type = gsub("_", " ", feature_type) %>% trimws())

# Group habitat into the final analysis categories

site_habitat_chamois <- master_presence %>%
  select(site_id, habitat) %>%
  distinct(site_id, .keep_all = TRUE) %>%
  filter(site_id %in% rownames(chamois_matrix_cov)) %>%
  arrange(match(site_id, rownames(chamois_matrix_cov))) %>%
  mutate(
    habitat_simple = case_when(
      habitat %in% c("Coniferous forest", "Mixed forest","Deciduous forest") ~ "Forest",
      habitat %in% c("Alpine pasture", "Mountain pasture","Meadow") ~ "Open_pasture",
      habitat == "Dwarf pine scrubs" ~ "Dwarf_pine",
      habitat == "Rocks" ~ "Rocks", TRUE ~ "Other"),
    habitat_simple = factor(habitat_simple, levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest", "Other")))

# Sheep capture rate per 100 trap nights

sheep_rate <- all_trap_rates_3 %>%
  filter(species == "sheep") %>%
  group_by(site_id) %>%
  summarise(
    events = sum(events, na.rm = TRUE),
    trap_nights = sum(trap_nights, na.rm = TRUE),
    #drop grouping after the site level summary.
    .groups = "drop") %>%
  mutate(
    sheep_rate_100tn = events / trap_nights * 100) %>%
  select(site_id, sheep_rate_100tn)

# Extract elevation for each camera site

coords_elev <- master_presence %>%
  distinct(site_id, .keep_all = TRUE) %>%
  filter(site_id %in% rownames(chamois_matrix_cov)) %>%
  arrange(match(site_id, rownames(chamois_matrix_cov))) %>%
  select(x = longitude, y = latitude) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

elev_data <- get_elev_point(coords_elev, src = "aws")

# Build the initial site covariate table

site_covs <- site_habitat_chamois %>%
  mutate(elevation_m = elev_data$elevation) %>%
  left_join(sheep_rate, by = "site_id") %>%
  left_join(trail_covariates %>%
      select(site_id, distance_to_trail_m), by = "site_id") %>%
  mutate(sheep_rate_100tn = replace_na(sheep_rate_100tn, 0))

site_covs <- site_covs %>%
  mutate(sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    log_distance_trail = log1p(distance_to_trail_m),
    log_distance_trail_sc = as.numeric(scale(log1p(distance_to_trail_m))),
    elevation_sc = as.numeric(scale(elevation_m)))

# Reclassify the two sites assigned to "Other" after field data review

other_sites <- site_habitat_chamois %>%
  filter(habitat_simple == "Other") %>%
  pull(site_id)

site_habitat_chamois <- site_habitat_chamois %>%
  mutate(
    habitat_simple = case_when(
      site_id == "1RET022" ~ "Open_pasture",
      site_id == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>%
      factor(
        #set the intended habitat reference level and category order.
        levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest", "Other")))

site_habitat_chamois <- site_habitat_chamois %>%
  mutate(habitat_simple = droplevels(habitat_simple))

# Rebuild site_covs after the two habitat reclassifications.
site_covs <- site_habitat_chamois %>%
  mutate(elevation_m = elev_data$elevation) %>%
  left_join(sheep_rate, by = "site_id") %>%
  left_join(trail_covariates %>%
      select(site_id, distance_to_trail_m),
    by = "site_id") %>%
  mutate(
    sheep_rate_100tn = replace_na(sheep_rate_100tn, 0))

site_covs <- site_covs %>%
  mutate(
    sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    log_distance_trail = log1p(distance_to_trail_m),
    log_distance_trail_sc = as.numeric(scale(log1p(distance_to_trail_m))),
    elevation_sc = as.numeric(scale(elevation_m)))

# Calculate marked trail length within multiple buffer sizes

camera_points <- master_presence %>%
  distinct(site_id, .keep_all = TRUE) %>%
  select(site_id, longitude, latitude) %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326) %>%
  st_transform(st_crs(trails_disturb))

buffer_sizes <- c(100, 250, 500, 1000)

trail_buffer_results <- list()

for (buf in buffer_sizes) {

  camera_buffers <- camera_points %>%
    st_buffer(dist = buf) %>%
    select(site_id)

  trails_clipped <- st_intersection(trails_disturb, camera_buffers)

  trail_lengths <- trails_clipped %>%
    mutate(trail_length_m = as.numeric(
        st_length(.))) %>%
    group_by(site_id, Disturb) %>%
    summarise(trail_length_m = sum(trail_length_m, na.rm = TRUE),
      .groups = "drop") %>%
    st_drop_geometry() %>%
    pivot_wider(names_from = Disturb, values_from = trail_length_m, values_fill = 0) %>%
    rename(high_use_trail_m = High, low_use_trail_m = Low) %>%
    mutate(total_trail_m = high_use_trail_m + low_use_trail_m, buffer_m = buf)

  all_cameras <- camera_points %>%
    st_drop_geometry() %>%
    select(site_id)

  trail_lengths_complete <- all_cameras %>%
    left_join(trail_lengths, by = "site_id") %>%
    mutate(buffer_m = buf, high_use_trail_m = replace_na(high_use_trail_m, 0),
      low_use_trail_m = replace_na(low_use_trail_m, 0),
      total_trail_m = replace_na(total_trail_m, 0))

  trail_buffer_results[[as.character(buf)]] <-trail_lengths_complete
  }

trail_buffer_table <- bind_rows(
  trail_buffer_results) %>%
  select(site_id, buffer_m, high_use_trail_m, low_use_trail_m, total_trail_m) %>%
  arrange(site_id, buffer_m)

# Select and scale the final 500 m trail length variable

trail_500 <- trail_buffer_table %>%
  filter(buffer_m == 500)

trail_500 <- trail_500 %>%
  mutate(trail_500_sc = as.numeric(scale(total_trail_m)))

site_covs <- site_covs %>%
  left_join(trail_500 %>%
      select(site_id, total_trail_m, trail_500_sc), by = "site_id")

# Site-covariate table supplied to the final spatial models.
site_covs_umf <- site_covs %>%
  select(habitat_simple, sheep_rate_sc, log_distance_trail_sc, trail_500_sc, elevation_sc)
