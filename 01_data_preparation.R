# 01_data_preparation.R

# Purpose
# Prepare deployment metadata, dog classifications, species
# presence indicators, and capture rates per 100 trap-nights.

# Inputs
# 2024 and 2025 deployment metadata
# Agouti deployment export3
# observation export containing dog classifications

# Outputs
# master_presence
# all_trap_rates_3

# Run this script first.

library(sf)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)

# Read deployment metadata

deployment_2025 <- st_read(
  "data_raw/antropic_disturbance_study/2025_deployment.shp")

deployment_2024 <- read_excel(
  "data_raw/antropic_disturbance_study/2024_Deployment.xlsx")

agouti_deployment <- read.csv(
  "data_raw/retezat_2024_2025/deployments.csv")


# Clean Agouti deployment dates

agouti_deployment <- agouti_deployment %>%
  mutate(
    deploymentStart_clean = na_if(deploymentStart, ""),
    deploymentStart_date = ymd_hms(deploymentStart_clean),
    year = year(deploymentStart_date))

agouti_deployment <- agouti_deployment %>%
  mutate(
    deploymentEnd_clean = na_if(deploymentEnd, ""),
    deploymentEnd_date = ymd_hms(deploymentEnd_clean),
    start_year = year(deploymentStart_date),
    end_year = year(deploymentEnd_date))

agouti_dep <- agouti_deployment %>%
  clean_names()

dep_2024 <- deployment_2024 %>%
  clean_names()

dep_2025_sf <- deployment_2025 %>%
  clean_names()

# Convert 2025 deployment coordinates to WGS84
dep_2025_sf_wgs84 <- st_transform(dep_2025_sf, 4326)

dep_2025_sf <- dep_2025_sf_wgs84 %>%
  mutate(
    longitude_from_shp = st_coordinates(.)[, 1],
    latitude_from_shp = st_coordinates(.)[, 2]
  ) %>% st_drop_geometry()



# Create master deployment table

agouti_dep <- agouti_dep %>%
  mutate(site_id = location_name %>% str_squish() %>% str_to_upper())

dep_2024 <- dep_2024 %>%
  mutate(site_id = location_id %>% str_squish() %>% str_to_upper())

dep_2025_sf <- dep_2025_sf %>%
  mutate(site_id = location_id %>% str_squish() %>% str_to_upper())

dep_2024_meta <- dep_2024 %>%
  transmute(
    site_id,
    longitude,
    latitude,
    feature_type = feature_typ,
    habitat = habitat_st)

dep_2025_meta <- dep_2025_sf %>%
  transmute(
    site_id,
    habitat,
    longitude = longitude_from_shp,
    latitude = latitude_from_shp,
    feature_type = feature_typ)

field_metadata <- bind_rows(dep_2024_meta, dep_2025_meta)

master_dep <- agouti_dep %>%
  left_join(field_metadata, by = "site_id")

master_dep_condense_2 <- master_dep %>%
  transmute(
    deployment_id,
    site_id,
    location_name,
    deployment_start_date,
    deployment_end_date,
    year,
    trap_nights = round(as.numeric(difftime(
          deployment_end_date,
          deployment_start_date,
          units = "days")), 1),
    latitude = latitude.y,
    longitude = longitude.y,
    habitat = habitat.y,
    feature_type = feature_type.y,
    agouti_feature_type = deployment_comments)

# Flag deployments with missing dates or non-positive camera effort
master_dep_condense_2 <- master_dep_condense_2 %>%
  mutate(
    effort_status = case_when(
      is.na(deployment_start_date) | is.na(deployment_end_date) ~ "missing_dates",
      trap_nights <= 0 ~ "zero_effort",
      TRUE ~ "usable"))

master_dep_condense_2 <- master_dep_condense_2 %>%
  mutate(
    start_month = month(deployment_start_date, label = TRUE, abbr = TRUE),
    end_month = month(deployment_end_date, label = TRUE, abbr = TRUE))

master_usable <- master_dep_condense_2 %>%
  filter(
    !is.na(trap_nights),
    trap_nights > 0)


# Read observations containing dog classifications (if using dog classifications in analysis)

obs_dogs <- read_csv(
  "data_raw/retezat_2024_2025_dogclasses/observations.csv")

# Parse dog types from the comments column
dogclasses <- obs_dogs %>%
  filter(scientificName == "Canis lupus familiaris") %>%
  mutate(
    dog_type_text = observationComments %>%
      str_remove("\\s*\\(.*\\)") %>%
      str_remove(regex("Dog\\s*type\\s*:\\s*", ignore_case = TRUE)) %>%
      str_trim(),
    has_L = replace_na(str_detect(dog_type_text, "\\bL\\b"), FALSE),
    has_H = replace_na(str_detect(dog_type_text, "\\bH\\b"), FALSE),
    has_P = replace_na(str_detect(dog_type_text, "\\bP\\b"), FALSE),
    has_U = replace_na(str_detect(dog_type_text, "\\bU\\b"), FALSE),
    pastoral_dog = has_L | has_H,
    pet_dog = has_P,
    unknown_dog = has_U)


# Species presence / absence by deployment
# Create species presence indicators by deployment

chamois_presence <- obs_dogs %>%
  filter(scientificName == "Rupicapra rupicapra") %>%
  distinct(deploymentID) %>%
  mutate(chamois_present = 1) %>%
  rename(deployment_id = deploymentID)

reddeer_presence <- obs_dogs %>%
  filter(scientificName == "Cervus elaphus") %>%
  distinct(deploymentID) %>%
  mutate(reddeer_present = 1) %>%
  rename(deployment_id = deploymentID)

roedeer_presence <- obs_dogs %>%
  filter(scientificName == "Capreolus capreolus") %>%
  distinct(deploymentID) %>%
  mutate(roedeer_present = 1) %>%
  rename(deployment_id = deploymentID)

redfox_presence <- obs_dogs %>%
  filter(scientificName == "Vulpes vulpes") %>%
  distinct(deploymentID) %>%
  mutate(redfox_present = 1) %>%
  rename(deployment_id = deploymentID)

brownbear_presence <- obs_dogs %>%
  filter(scientificName == "Ursus arctos") %>%
  distinct(deploymentID) %>%
  mutate(brownbear_present = 1) %>%
  rename(deployment_id = deploymentID)

sheep_presence <- obs_dogs %>%
  filter(scientificName == "Ovis aries") %>%
  distinct(deploymentID) %>%
  mutate(sheep_present = 1) %>%
  rename(deployment_id = deploymentID)

human_presence <- obs_dogs %>%
  filter(scientificName == "Homo sapiens") %>%
  distinct(deploymentID) %>%
  mutate(human_present = 1) %>%
  rename(deployment_id = deploymentID)

dog_any_presence <- obs_dogs %>%
  filter(scientificName == "Canis lupus familiaris") %>%
  distinct(deploymentID) %>%
  mutate(dog_any_present = 1) %>%
  rename(deployment_id = deploymentID)

pastoral_dog_presence <- dogclasses %>%
  filter(pastoral_dog == TRUE) %>%
  distinct(deploymentID) %>%
  mutate(pastoral_dog_present = 1) %>%
  rename(deployment_id = deploymentID)

pet_dog_presence <- dogclasses %>%
  filter(pet_dog == TRUE) %>%
  distinct(deploymentID) %>%
  mutate(pet_dog_present = 1) %>%
  rename(deployment_id = deploymentID)

unknown_dog_presence <- dogclasses %>%
  filter(unknown_dog == TRUE) %>%
  distinct(deploymentID) %>%
  mutate(unknown_dog_present = 1) %>%
  rename(deployment_id = deploymentID)


# Master presence / absence table

master_presence <- master_usable %>%
  left_join(brownbear_presence, by = "deployment_id") %>%
  left_join(chamois_presence, by = "deployment_id") %>%
  left_join(dog_any_presence, by = "deployment_id") %>%
  left_join(human_presence, by = "deployment_id") %>%
  left_join(pastoral_dog_presence, by = "deployment_id") %>%
  left_join(pet_dog_presence, by = "deployment_id") %>%
  left_join(unknown_dog_presence, by = "deployment_id") %>%
  left_join(reddeer_presence, by = "deployment_id") %>%
  left_join(roedeer_presence, by = "deployment_id") %>%
  left_join(redfox_presence, by = "deployment_id") %>%
  left_join(sheep_presence, by = "deployment_id") %>%
  mutate(
    across(ends_with("_present"),
           ~ replace_na(.x, 0)))

# Summarise the number of deployments with each presence indicator
master_presence %>%
  select(ends_with("_present")) %>%
  summarise(across(everything(), sum))


# Trap-rate functions
# Functions used to calculate capture rates per 100 trap-nights

make_trap_rate <- function(species_name, species_label) {

  species_events <- obs_dogs %>%
    filter(scientificName == species_name) %>%
    mutate(
      event_key = if_else(
        is.na(eventID) | eventID == "",
        as.character(observationID),
        as.character(eventID))) %>%
    group_by(deploymentID) %>%
    summarise(
      events = n_distinct(event_key),
      .groups = "drop") %>%
    rename(deployment_id = deploymentID)

  master_presence %>%
    select(deployment_id, site_id, trap_nights) %>%
    left_join(species_events, by = "deployment_id") %>%
    mutate(
      events = replace_na(events, 0),
      rate_100tn = events / trap_nights * 100,
      species = species_label) %>%
    select(species, deployment_id, site_id, trap_nights, events, rate_100tn)}

all_trap_rates <- bind_rows(
  make_trap_rate("Rupicapra rupicapra", "chamois"),
  make_trap_rate("Cervus elaphus", "red_deer"),
  make_trap_rate("Capreolus capreolus", "roe_deer"),
  make_trap_rate("Ursus arctos", "brown_bear"),
  make_trap_rate("Vulpes vulpes", "red_fox"),
  make_trap_rate("Homo sapiens", "human"),
  make_trap_rate("Ovis aries", "sheep"),
  make_trap_rate("Canis lupus familiaris", "any_dog"))

make_dogtype_trap_rate <- function(dog_column, dog_label) {

  dog_events <- dogclasses %>%
    filter(.data[[dog_column]] == TRUE) %>%
    mutate(
      event_key = if_else(
        is.na(eventID) | eventID == "",
        as.character(observationID),
        as.character(eventID))) %>%
    group_by(deploymentID) %>%
    summarise(
      events = n_distinct(event_key),
      .groups = "drop") %>%
    rename(deployment_id = deploymentID)

  master_presence %>%
    select(deployment_id, site_id, trap_nights) %>%
    left_join(dog_events, by = "deployment_id") %>%
    mutate(
      events = replace_na(events, 0),
      rate_100tn = events / trap_nights * 100,
      species = dog_label) %>%
    select(species, deployment_id, site_id, trap_nights, events, rate_100tn)}

dogtype_trap_rates <- bind_rows(
  make_dogtype_trap_rate("pastoral_dog", "pastoral_dog"),
  make_dogtype_trap_rate("pet_dog", "pet_dog"),
  make_dogtype_trap_rate("unknown_dog", "unknown_dog"))

all_trap_rates_2 <- bind_rows(
  all_trap_rates, dogtype_trap_rates)

# Add an all-dog flag and calculate the final dog-type capture rates
dogclasses <- dogclasses %>%
  mutate(any_dog = TRUE)

dogtype_trap_rates <- bind_rows(
  make_dogtype_trap_rate("pastoral_dog", "pastoral_dog"),
  make_dogtype_trap_rate("pet_dog", "pet_dog"),
  make_dogtype_trap_rate("unknown_dog", "unknown_dog"),
  make_dogtype_trap_rate("any_dog", "any_dog"))

all_trap_rates_3 <- bind_rows(
  all_trap_rates %>%
    filter(species != "any_dog"),
  dogtype_trap_rates)

trap_rate_summary_final <- all_trap_rates_3 %>%
  group_by(species) %>%
  summarise(
    n_sites = n(),
    n_sites_detected = sum(events > 0),
    total_events = sum(events),
    mean_rate_100tn = mean(rate_100tn, na.rm = TRUE),
    median_rate_100tn = median(rate_100tn, na.rm = TRUE),
    max_rate_100tn = max(rate_100tn, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(total_events))

trap_rate_summary_final
