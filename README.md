# GBD-CKD-T2DM-analysis
R code for burden and risk-attributable fraction analysis of diabetic kidney disease (DKD) among people with type 2 diabetes mellitus, based on the Global Burden of Disease 2023 dataset.

## Overview
This repository contains reproducible analytical scripts to estimate DALY, YLL, YLD, population attributable fraction (PAF), sensitivity analysis, health inequality metrics, and global spatial mapping for DKD attributable to modifiable risk factors.

## Scripts
- `gbd_burden_analysis.R`: Main disease burden estimation (DALY/YLL/YLD)
- `case_control_paf.R`: Calculation of population attributable fractions
- `joint_paf_sensitivity.R`: Sensitivity analysis for joint risk factor PAF
- `inequality_metrics.R`: Health inequality assessment
- `mapping_global.R`: Global geographical mapping and visualization
- `bpc_forecast.R`: Burden forecasting

## Environment & Dependencies
R (≥ 4.1). Required packages include `tidyverse`, `ggplot2`, `sf`, `gbdtools`, and related analytical libraries.

## Data availability
GBD raw data are publicly available from the Institute for Health Metrics and Evaluation (IHME). Individual-level clinical data cannot be shared publicly owing to ethical and privacy restrictions.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use these codes in your work, please cite this repository and our associated article upon publication.
