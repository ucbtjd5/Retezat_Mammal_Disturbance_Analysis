# 02_detection_matrices.R

# Purpose
# Create 7-day count matrices and species-specific effort
# matrices for chamois, red deer, roe deer, red fox, and
# brown bear.

# Inputs
# all_trap_rates_3 from 01_data_preparation.R
# Agouti camtrapDP datapackage

# Outputs
# five 104-site x 34-occasion count matrices
# species-specific log-effort matrices

# Final spatial models are fitted in 04_spatial_models.R.

library(camtrapdp)
library(dplyr)
library(tidyr)
library(lubridate)

# Marcus Rowcliffe's helper for constructing detection/count matrices
source("https://raw.githubusercontent.com/MarcusRowcliffe/make_detection_matrix/refs/heads/main/make_detection_matrix.R")


# Read Agouti camtrapDP data

ag <- read_camtrapdp(
  "data_raw/retezat_2024_2025/datapackage.json")


# Clean deployment dates inside the camtrap object

ag_clean <- ag

ag_clean$data$deployments <- ag_clean$data$deployments %>%
  mutate(
    deploymentStart = na_if(as.character(deploymentStart), ""),
    deploymentEnd = na_if(as.character(deploymentEnd), ""),
    deploymentStart = ymd_hms(deploymentStart, quiet = TRUE),
  #convert deployment dates to date-time format
    deploymentEnd = ymd_hms(deploymentEnd, quiet = TRUE)) %>%
  filter(
    !is.na(deploymentStart),
    !is.na(deploymentEnd),
    deploymentEnd > deploymentStart)

# Keep observations only from valid deployments
ag_clean$data$observations <- ag_clean$data$observations %>%
  filter(
    deploymentID %in% ag_clean$data$deployments$deploymentID)


# CHAMOIS: binary matrix used only to define final sites/occasions

chamois <- make_detection_matrix(
  ag_clean,
  species = "Rupicapra rupicapra",
  interval = 7)

chamois_matrix <- chamois$matrix$`Rupicapra rupicapra`

# Remove fully empty occasions
empty_occasions <- colSums(!is.na(chamois_matrix)) == 0
chamois_matrix_model <- chamois_matrix[, !empty_occasions]

# Remove 2RET071 (approximately 8-minute deployment)
chamois_matrix_cov <- chamois_matrix_model[rownames(chamois_matrix_model) != "2RET071", ]

# Site-level disturbance rates used to retain the valid final site set
site_rate_covs <- all_trap_rates_3 %>%
  filter(
    species %in% c("human", "sheep", "pastoral_dog")) %>%
  group_by(site_id, species) %>%
  summarise(
    events = sum(events, na.rm = TRUE),
    trap_nights = sum(trap_nights, na.rm = TRUE),
    rate_100tn = events / trap_nights * 100,
    .groups = "drop") %>%
  select(site_id, species, rate_100tn) %>%
  pivot_wider(
    names_from = species,
    values_from = rate_100tn,
    names_glue = "{species}_rate_100tn")

keep_sites <- rownames(chamois_matrix_model)[
  rownames(chamois_matrix_model) %in% site_rate_covs$site_id]

chamois_matrix_cov <- chamois_matrix_model[keep_sites, , drop = FALSE]

# After removing 2RET071, occasion 51 becomes empty
chamois_matrix_cov <- chamois_matrix_cov[, colSums(!is.na(chamois_matrix_cov)) > 0, drop = FALSE]


# CHAMOIS effort matrix

effort_matrix <- chamois$effort

effort_clean <- effort_matrix[rownames(effort_matrix) != "2RET071", ]

effort_clean <- effort_clean[, !empty_occasions]

# Remove occasion 51 / align exactly to final chamois matrix
effort_clean <- effort_clean[, colnames(chamois_matrix_cov)]

log_effort <- log(effort_clean + 0.01)


# CHAMOIS count matrix

chamois_count <- make_detection_matrix(
  ag_clean,
  species = "Rupicapra rupicapra",
  interval = 7,
  type = "count")

chamois_count_matrix <- chamois_count$matrix$`Rupicapra rupicapra`

empty_occ <- colSums(!is.na(chamois_count_matrix)) == 0
chamois_count_matrix <- chamois_count_matrix[, !empty_occ]

chamois_count_matrix <- chamois_count_matrix[rownames(chamois_count_matrix) != "2RET071", ]

# Match the final 34 occasions
chamois_count_matrix <- chamois_count_matrix[, colnames(chamois_count_matrix) %in% colnames(chamois_matrix_cov)]

# Keep exact occasion order
chamois_count_matrix <- chamois_count_matrix[, colnames(chamois_matrix_cov)]

# RED DEER

red_deer <- make_detection_matrix(
  ag_clean,
  species = "Cervus elaphus",
  interval = 7,
  type = "count")

red_deer_matrix <- red_deer$matrix$`Cervus elaphus`

empty_occ_rd <- colSums(!is.na(red_deer_matrix)) == 0
red_deer_matrix <- red_deer_matrix[, !empty_occ_rd]

red_deer_matrix <- red_deer_matrix[rownames(red_deer_matrix) != "2RET071", ]

# Match final chamois occasions
red_deer_matrix <- red_deer_matrix[, colnames(chamois_count_matrix)]

# Species-specific camera-effort matrix
effort_rd <- red_deer$effort

effort_rd <- effort_rd[rownames(effort_rd) != "2RET071", ]

effort_rd <- effort_rd[, colnames(red_deer_matrix)]

log_effort_rd <- log(effort_rd + 0.01)

# ROE DEER

roe_deer <- make_detection_matrix(
  ag_clean,
  species = "Capreolus capreolus",
  interval = 7,
  type = "count")

roe_deer_matrix <- roe_deer$matrix$`Capreolus capreolus`

empty_occ_roe <- colSums(!is.na(roe_deer_matrix)) == 0
roe_deer_matrix <- roe_deer_matrix[, !empty_occ_roe]

roe_deer_matrix <- roe_deer_matrix[
  rownames(roe_deer_matrix) != "2RET071", ]

roe_deer_matrix <- roe_deer_matrix[, colnames(chamois_count_matrix)]

effort_roe <- roe_deer$effort

effort_roe <- effort_roe[rownames(effort_roe) != "2RET071", ]

effort_roe <- effort_roe[, colnames(roe_deer_matrix)]

log_effort_roe <- log(effort_roe + 0.01)

# RED FOX

red_fox <- make_detection_matrix(
  ag_clean,
  species = "Vulpes vulpes",
  interval = 7,
  type = "count")

red_fox_matrix <- red_fox$matrix$`Vulpes vulpes`

empty_occ_red_fox <- colSums(!is.na(red_fox_matrix)) == 0
red_fox_matrix <- red_fox_matrix[, !empty_occ_red_fox]

red_fox_matrix <- red_fox_matrix[rownames(red_fox_matrix) != "2RET071", ]

red_fox_matrix <- red_fox_matrix[, colnames(chamois_count_matrix)]

effort_red_fox <- red_fox$effort

effort_red_fox <- effort_red_fox[rownames(effort_red_fox) != "2RET071", ]

effort_red_fox <- effort_red_fox[, colnames(red_fox_matrix)]

log_effort_red_fox <- log(effort_red_fox + 0.01)

# BROWN BEAR

bear <- make_detection_matrix(
  ag_clean,
  species = "Ursus arctos",
  interval = 7,
  type = "count")

bear_matrix <- bear$matrix$`Ursus arctos`

empty_occ_bear <- colSums(!is.na(bear_matrix)) == 0
bear_matrix <- bear_matrix[, !empty_occ_bear]

bear_matrix <- bear_matrix[rownames(bear_matrix) != "2RET071", ]

bear_matrix <- bear_matrix[, colnames(chamois_count_matrix)]

effort_bear <- bear$effort

effort_bear <- effort_bear[rownames(effort_bear) != "2RET071", ]

effort_bear <- effort_bear[, colnames(bear_matrix)]

log_effort_bear <- log(effort_bear + 0.01)
