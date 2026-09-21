# 00_run_all.R

# Convenience script for reproducing the core analysis.
# Run from the project root containing r_clean/ and data_raw/.
# No additional analysis is performed here.

source("r_clean/01_data_preparation.R")
source("r_clean/02_detection_matrices.R")
source("r_clean/03_site_covariates.R")
source("r_clean/04_spatial_models.R")
source("r_clean/05_temporal_models.R")