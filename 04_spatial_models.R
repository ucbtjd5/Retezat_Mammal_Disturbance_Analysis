# 04_spatial_models.R

# Purpose
# Fit the final N mixture candidate model set separately for
# brown bear, red deer, red fox, chamois, and roe deer.

# Predictors:
# habitat and elevation (environment)
# sheep capture rate (pastoral disturbance)
# marked trail length within 500 m (recreation)

# Detection is modelled using log camera effort.

# Interpretation
# Although unmarked labels the latent state as "Abundance",
# it is interpreted as relative site-use intensity rather than
# absolute abundance.


library(unmarked)

# Brown bear spatial model

umf_bear <- unmarkedFramePCount(
  y = as.matrix(bear_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_bear))

nm0_bear <- pcount(
  ~ log_effort ~ 1,
  data = umf_bear, K = 100)

nm_envr_bear <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc,
  data = umf_bear, K = 100)

nm_pastoral_bear <- pcount(
  ~ log_effort ~ sheep_rate_sc,
  data = umf_bear, K = 100)

nm_recreation_bear <- pcount(
  ~ log_effort ~ trail_500_sc,
  data = umf_bear, K = 100)

nm_envr_pastoral_bear <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + sheep_rate_sc,
  data = umf_bear, K = 100)

nm_envr_recreation_bear <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + trail_500_sc,
  data = umf_bear, K = 100)

nm_recreation_pastoral_bear <- pcount(
  ~ log_effort ~ trail_500_sc + sheep_rate_sc,
  data = umf_bear, K = 100)

nm_full_bear <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc +
    trail_500_sc + sheep_rate_sc,
  data = umf_bear, K = 100)

nm_model_list_bear <- fitList(
  "null"                = nm0_bear,
  "envr"                = nm_envr_bear,
  "pastoral"            = nm_pastoral_bear,
  "recreation"          = nm_recreation_bear,
  "envr+pastoral"       = nm_envr_pastoral_bear,
  "envr+recreation"     = nm_envr_recreation_bear,
  "recreation+pastoral" = nm_recreation_pastoral_bear,
  "full"                = nm_full_bear)

nm_model_table_bear <- modSel(nm_model_list_bear)
print(nm_model_table_bear)

# Final interpretation used the simpler envr+recreation model because
# the full model was only 1.34 AIC lower and added sheep rate.
summary(nm_envr_recreation_bear)


# Red deer spatial model

umf_rd <- unmarkedFramePCount(
  y = as.matrix(red_deer_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_rd))

nm0_rd <- pcount(
  ~ log_effort ~ 1,
  data = umf_rd, K = 100)

nm_envr_rd <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc,
  data = umf_rd, K = 100)

nm_pastoral_rd <- pcount(
  ~ log_effort ~ sheep_rate_sc,
  data = umf_rd, K = 100)

nm_recreation_rd <- pcount(
  ~ log_effort ~ trail_500_sc,
  data = umf_rd, K = 100)

nm_envr_pastoral_rd <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + sheep_rate_sc,
  data = umf_rd, K = 100)

nm_envr_recreation_rd <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + trail_500_sc,
  data = umf_rd, K = 100)

nm_recreation_pastoral_rd <- pcount(
  ~ log_effort ~ trail_500_sc + sheep_rate_sc,
  data = umf_rd, K = 100)

nm_full_rd <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc +
    trail_500_sc + sheep_rate_sc,
  data = umf_rd, K = 100)

nm_model_list_rd <- fitList(
  "null"                = nm0_rd,
  "envr"                = nm_envr_rd,
  "pastoral"            = nm_pastoral_rd,
  "recreation"          = nm_recreation_rd,
  "envr+pastoral"       = nm_envr_pastoral_rd,
  "envr+recreation"     = nm_envr_recreation_rd,
  "recreation+pastoral" = nm_recreation_pastoral_rd,
  "full"                = nm_full_rd)

nm_model_table_rd <- modSel(nm_model_list_rd)
print(nm_model_table_rd)

summary(nm_envr_recreation_rd)


# Red fox spatial model

umf_red_fox <- unmarkedFramePCount(
  y = as.matrix(red_fox_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_red_fox))

nm0_red_fox <- pcount(
  ~ log_effort ~ 1,
  data = umf_red_fox, K = 100)

nm_envr_red_fox <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc,
  data = umf_red_fox, K = 100)

nm_pastoral_red_fox <- pcount(
  ~ log_effort ~ sheep_rate_sc,
  data = umf_red_fox, K = 100)

nm_recreation_red_fox <- pcount(
  ~ log_effort ~ trail_500_sc,
  data = umf_red_fox, K = 100)

nm_envr_pastoral_red_fox <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + sheep_rate_sc,
  data = umf_red_fox, K = 100)

nm_envr_recreation_red_fox <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + trail_500_sc,
  data = umf_red_fox, K = 100)

nm_recreation_pastoral_red_fox <- pcount(
  ~ log_effort ~ trail_500_sc + sheep_rate_sc,
  data = umf_red_fox, K = 100)

nm_full_red_fox <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc +
    trail_500_sc + sheep_rate_sc,
  data = umf_red_fox, K = 100)

nm_model_list_red_fox <- fitList(
  "null"                = nm0_red_fox,
  "envr"                = nm_envr_red_fox,
  "pastoral"            = nm_pastoral_red_fox,
  "recreation"          = nm_recreation_red_fox,
  "envr+pastoral"       = nm_envr_pastoral_red_fox,
  "envr+recreation"     = nm_envr_recreation_red_fox,
  "recreation+pastoral" = nm_recreation_pastoral_red_fox,
  "full"                = nm_full_red_fox)

nm_model_table_red_fox <- modSel(nm_model_list_red_fox)
print(nm_model_table_red_fox)

# envr+pastoral has the lowest AIC, but the environment-only model is
# retained for interpretation under the stated parsimony rule (Delta AIC = 0.52).
# The summary below is retained to inspect the sheep coefficient in the more
# complex model.
summary(nm_envr_pastoral_red_fox)


# Chamois spatial model

umf_chamois_count <- unmarkedFramePCount(
  y = as.matrix(chamois_count_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort))

nm0 <- pcount(
  ~ log_effort ~ 1,
  data = umf_chamois_count, K = 100)

nm_envr <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc,
  data = umf_chamois_count, K = 100)

nm_pastoral <- pcount(
  ~ log_effort ~ sheep_rate_sc,
  data = umf_chamois_count, K = 100)

nm_recreation <- pcount(
  ~ log_effort ~ trail_500_sc,
  data = umf_chamois_count, K = 100)

nm_envr_pastoral <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + sheep_rate_sc,
  data = umf_chamois_count, K = 100)

nm_envr_recreation <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + trail_500_sc,
  data = umf_chamois_count, K = 100)

nm_recreation_pastoral <- pcount(
  ~ log_effort ~ trail_500_sc + sheep_rate_sc,
  data = umf_chamois_count, K = 100)

nm_full <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc +
    trail_500_sc + sheep_rate_sc,
  data = umf_chamois_count, K = 100)

nm_model_list_chamois <- fitList(
  "null"                = nm0,
  "envr"                = nm_envr,
  "pastoral"            = nm_pastoral,
  "recreation"          = nm_recreation,
  "envr+pastoral"       = nm_envr_pastoral,
  "envr+recreation"     = nm_envr_recreation,
  "recreation+pastoral" = nm_recreation_pastoral,
  "full"                = nm_full)

nm_model_table_chamois <- modSel(nm_model_list_chamois)
print(nm_model_table_chamois)

summary(nm_envr_recreation)


# Roe deer spatial model

# Final spatial analysis retains all 104 sites and four habitat categories.
# The 90 site Forest/Non-forest habitat formula was used only as a sensitivity check.

umf_roe <- unmarkedFramePCount(
  y = as.matrix(roe_deer_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_roe))

nm0_roe <- pcount(
  ~ log_effort ~ 1,
  data = umf_roe, K = 100)

nm_envr_roe <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc,
  data = umf_roe, K = 100)

nm_pastoral_roe <- pcount(
  ~ log_effort ~ sheep_rate_sc,
  data = umf_roe, K = 100)

nm_recreation_roe <- pcount(
  ~ log_effort ~ trail_500_sc,
  data = umf_roe, K = 100)

nm_envr_pastoral_roe <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + sheep_rate_sc,
  data = umf_roe, K = 100)

nm_envr_recreation_roe <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc + trail_500_sc,
  data = umf_roe, K = 100)

nm_recreation_pastoral_roe <- pcount(
  ~ log_effort ~ trail_500_sc + sheep_rate_sc,
  data = umf_roe, K = 100)

nm_full_roe <- pcount(
  ~ log_effort ~ habitat_simple + elevation_sc +
    trail_500_sc + sheep_rate_sc,
  data = umf_roe, K = 100)

nm_model_list_roe <- fitList(
  "null"                = nm0_roe,
  "envr"                = nm_envr_roe,
  "pastoral"            = nm_pastoral_roe,
  "recreation"          = nm_recreation_roe,
  "envr+pastoral"       = nm_envr_pastoral_roe,
  "envr+recreation"     = nm_envr_recreation_roe,
  "recreation+pastoral" = nm_recreation_pastoral_roe,
  "full"                = nm_full_roe)

nm_table_roe <- modSel(nm_model_list_roe)
print(nm_table_roe)

summary(nm_full_roe)
