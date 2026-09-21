# 05_temporal_models.R

# Purpose
# Fit circular hierarchical trigonometric models describing diel activity
# patterns and their relationships with pastoral disturbance
# and recreational trail infrastructure.

# Inputs
# cleaned Agouti deployments and observations
# sheep capture rate
# total marked trail length within 500 m
# habitat

# Brown bear, red deer, chamois, and red fox use 104 sites.
# Roe deer temporal models models use 90 sites after excluding rocky
# sites with no roe deer detections.


library(dplyr)
library(tidyr)
library(lubridate)
library(GLMMadaptive)

# Marcus Rowcliffe CHM helper functions
source("https://raw.githubusercontent.com/MarcusRowcliffe/chm_helpers/main/chm_helpers.R")

# Prepare deployment and observation tables for temporal models

# Sheep detection rate per deployment
sheep_rate_per_deployment <- all_trap_rates_3 %>%
  filter(species == "sheep") %>%
  group_by(deployment_id) %>%
  summarise(
    sheep_rate_100tn = sum(events) / sum(trap_nights) * 100,
    .groups = "drop") %>%
  rename(deploymentID = deployment_id)

# Deployment table used by make_chm_data
deployment_covariates_chm <- ag_clean$data$deployments %>%
  select(
    deploymentID,
    locationName,
    deploymentStart,
    deploymentEnd) %>%
  left_join(
    sheep_rate_per_deployment,
    by = "deploymentID") %>%
  mutate(
    sheep_rate_100tn = replace_na(sheep_rate_100tn, 0))

# Observation table used by make_chm_data
observations_chm <- ag_clean$data$observations %>%
  filter(!is.na(eventStart)) %>%
  select(
    deploymentID,
    timestamp = eventStart,
    scientificName)

# Add the final habitat grouping to the deployment table
deployment_covariates_chm <- deployment_covariates_chm %>%
  left_join(
    master_presence %>%
      distinct(site_id, .keep_all = TRUE) %>%
      select(site_id, habitat) %>%
      mutate(
        habitat = gsub("_", " ", habitat) %>%
          trimws(),
        habitat_simple = case_when(
          habitat %in% c(
            "Coniferous forest", "Mixed forest", "Deciduous forest"
          ) ~ "Forest",
          habitat %in% c(
            "Alpine pasture", "Mountain pasture", "Meadow"
          ) ~ "Open_pasture",
          habitat == "Dwarf pine scrubs" ~ "Dwarf_pine",
          habitat == "Rocks" ~ "Rocks",
          TRUE ~ "Other") %>%
          factor(levels = c(
              "Open_pasture", "Rocks", "Dwarf_pine", "Forest", "Other"))), 
    by = c("locationName" = "site_id"))

# Add the final 500 m trail variable and remove 2RET071

deployment_covariates_chm <- deployment_covariates_chm %>%
  left_join(trail_500 %>%
      select(site_id, total_trail_m), by = c("locationName" = "site_id"))

deployment_covariates_chm <- deployment_covariates_chm %>%
  filter(locationName != "2RET071")


# Brown bear temporal model

chm_data_bear <- make_chm_data(
  deployments = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Ursus arctos"),
  nBins = 24, covs = c(
    "sheep_rate_100tn",
    "total_trail_m",
    "habitat_simple"), collapse = TRUE) %>%
  mutate(
    sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    trail_500_sc = as.numeric(scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>%
      factor(levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest")))

chm_null_bear <- mixed_model(
  fixed = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_unimodal_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_bimodal_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_bimodal_hab_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_sheep_hab_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_trail_hab_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

chm_trail_sheep_hab_bear <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_bear)

AIC(
  chm_null_bear,
  chm_unimodal_bear,
  chm_bimodal_bear,
  chm_bimodal_hab_bear,
  chm_sheep_hab_bear,
  chm_trail_hab_bear,
  chm_trail_sheep_hab_bear)

summary(chm_sheep_hab_bear)


# Red deer temporal model

chm_data_red_deer <- make_chm_data(
  deployments = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Cervus elaphus"),
  nBins = 24, covs = c(
    "sheep_rate_100tn",
    "total_trail_m",
    "habitat_simple"), collapse = TRUE) %>%
  mutate(
    sheep_rate_sc = as.numeric(
      scale(sheep_rate_100tn)),
    trail_500_sc = as.numeric(
      scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>%
      factor(levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest")))

chm_null_rd <- mixed_model(
  fixed = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_unimodal_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_bimodal_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_bimodal_hab_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_sheep_hab_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_trail_hab_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

chm_trail_sheep_hab_rd <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_red_deer)

AIC(
  chm_null_rd,
  chm_unimodal_rd,
  chm_bimodal_rd,
  chm_bimodal_hab_rd,
  chm_sheep_hab_rd,
  chm_trail_hab_rd,
  chm_trail_sheep_hab_rd)

summary(chm_sheep_hab_rd)


# Chamois temporal model

chm_data_chamois <- make_chm_data(
  deployments = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Rupicapra rupicapra"),
  nBins = 24,
  covs = c(
    "sheep_rate_100tn",
    "total_trail_m",
    "habitat_simple"), collapse = TRUE) %>%
  mutate(
    sheep_rate_sc = as.numeric(
      scale(sheep_rate_100tn)),
    trail_500_sc = as.numeric(scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>%
      factor(levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest")))

chm_null_chamois <- mixed_model(
  fixed = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_unimodal_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_bimodal_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_bimodal_hab_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_sheep_hab_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_trail_hab_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

chm_trail_sheep_hab_chamois <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_chamois)

AIC(
  chm_null_chamois,
  chm_unimodal_chamois,
  chm_bimodal_chamois,
  chm_bimodal_hab_chamois,
  chm_sheep_hab_chamois,
  chm_trail_hab_chamois,
  chm_trail_sheep_hab_chamois)

summary(chm_trail_hab_chamois)


# Red fox temporal model

chm_data_fox <- make_chm_data(
  deployments = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Vulpes vulpes"),
  nBins = 24,
  covs = c(
    "sheep_rate_100tn",
    "total_trail_m",
    "habitat_simple"), collapse = TRUE) %>%
  mutate(
    sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    trail_500_sc = as.numeric(scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple) ) %>%
      factor(levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest")))

chm_null_fox <- mixed_model(
  fixed = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_unimodal_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_bimodal_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_bimodal_hab_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_sheep_hab_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_trail_hab_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

chm_trail_sheep_hab_fox <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2 * timeRadian) * sheep_rate_sc +
    sin(2 * timeRadian) * sheep_rate_sc +
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_fox)

AIC(
  chm_null_fox,
  chm_unimodal_fox,
  chm_bimodal_fox,
  chm_bimodal_hab_fox,
  chm_sheep_hab_fox,
  chm_trail_hab_fox,
  chm_trail_sheep_hab_fox)

summary(chm_sheep_hab_fox)


# Roe deer temporal model

# Roe deer had no detections at the 14 rocky sites; exclude these
# sites from the temporal analysis.
rocks_sites <- site_habitat_chamois %>%
  filter(habitat_simple == "Rocks") %>%
  pull(site_id)

# Sheep is not modelled temporally for roe deer because roe deer
# and sheep co-occurred at only three sites.
chm_data_roe <- make_chm_data(
  deployments = deployment_covariates_chm %>%
    filter(!locationName %in% rocks_sites),
  observations = observations_chm %>%
    filter(scientificName == "Capreolus capreolus"),
  nBins = 24,
  covs = c(
    "total_trail_m",
    "habitat_simple"), collapse = TRUE) %>%
  mutate(
    trail_500_sc = as.numeric(
      scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)),
    habitat_roe = case_when(
      habitat_simple == "Forest" ~ "Forest",
      TRUE ~ "Non-forest") %>%
      factor(levels = c("Non-forest", "Forest")))

chm_null_roe <- mixed_model(
  fixed = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_roe)

chm_unimodal_roe <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_roe)

chm_bimodal_roe <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_roe)

chm_bimodal_hab_roe <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2 * timeRadian) + sin(2 * timeRadian) +
    habitat_roe,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_roe)

chm_trail_hab_roe <- mixed_model(
  fixed = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2 * timeRadian) * trail_500_sc +
    sin(2 * timeRadian) * trail_500_sc +
    habitat_roe,
  random = ~ 1 | locationName,
  family = binomial(),
  data = chm_data_roe)

AIC(
  chm_null_roe,
  chm_unimodal_roe,
  chm_bimodal_roe,
  chm_bimodal_hab_roe,
  chm_trail_hab_roe)

summary(chm_bimodal_hab_roe)
