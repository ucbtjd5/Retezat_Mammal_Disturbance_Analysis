# Retezat mammal disturbance analysis

R code for the final analyses on large mammal responses to pastoral and recreational disturbance in Retezat National Park, Romania.

The cleaned scripts were rerun from a fresh R session and reproduce the final spatial and temporal results used in the thesis.

## Scripts

- `01_data_preparation.R` combines deployment metadata, calculates trap nights, and prepares site level detection/capture rate information used later in the analysis.
- `02_detection_matrices.R` creates 7 day count matrices and species specific camera effort matrices for chamois, red deer, roe deer, red fox, and brown bear.
- `03_site_covariates.R` prepares habitat, elevation, sheep capture rate, and marked trail length around camera sites. The final recreational variable is total marked trail length within 500 m.
- `04_spatial_models.R` fits the final N mixture candidate models for the five focal species using `unmarked::pcount`.
- `05_temporal_models.R` fits the final diel activity models using circular trigonometric terms.
- `00_run_all.R` runs scripts 01–05 in order.

## Final analysis

The final dataset contains 104 camera sites and 34 seven day sampling occasions.

The five focal species are:
- chamois (*Rupicapra rupicapra*)
- red deer (*Cervus elaphus*)
- roe deer (*Capreolus capreolus*)
- red fox (*Vulpes vulpes*)
- brown bear (*Ursus arctos*)

The main site covariates are habitat, elevation, sheep capture rate per 100 trap nights, and total marked trail length within 500 m of each camera.

Sheep capture rate is used as the pastoral disturbance proxy. Pastoral dog detections were not included as separate predictors in the final models.

For the spatial analysis, camera effort is included in the detection model. Although `unmarked` labels the parameter as abundance, it is interpreted here as relative site use intensity rather than absolute abundance.

The temporal analysis uses 24 hourly bins. Roe deer temporal models use 90 sites because there were no roe deer detections at the 14 rocky sites.

## Data

Raw camera trap and spatial data are not included in this repository.

To rerun the analysis, the following inputs are expected in `data_raw/`:

- `antropic_disturbance_study/2024_Deployment.xlsx`
- `antropic_disturbance_study/2025_deployment.shp` and its associated shapefile files
- `retezat_2024_2025/datapackage.json` and the complete camtrapDP package it references
- `retezat_2024_2025/deployments.csv`
- `retezat_2024_2025_dogclasses/observations.csv`
- `Trail_disturbance_correct/Trail_disturbance_correct.shp` and its associated shapefile files

## R packages

The workflow uses:

- `sf`
- `tidyverse`
- `readxl`
- `lubridate`
- `janitor`
- `camtrapdp`
- `elevatr`
- `unmarked`
- `GLMMadaptive`

The scripts also source Marcus Rowcliffe's `make_detection_matrix.R` and `chm_helpers.R` helper files from GitHub.

## Running the analysis

Run at once:

```r
source("r_clean/00_run_all.R")
```

Alternatively, scripts 01–05 can be used individually in numerical order.

## Assisted cleanup

Claude, Sonnet 5 was used to help reorganise and document the original R scripts and to troubleshoot reproducibility during cleanup. The analytical choices, model specifications, covariates, model selection, and interpretation were based on the original dissertation analysis.

## Scope

This repository contains the core code needed to reproduce the final analyses. Earlier exploratory code, superseded model versions, scratch diagnostics, and mapping experiments are not included.
