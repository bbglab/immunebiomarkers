# Exhaustive study

## Goal
* Run an exhaustive association study relating survival and response outcomes to all features.

## Methodology
* Best overall response associations modelled with logistic regressions.
* Overall survival and progression free surival modeled with cox proportional hazards models. 
* Testing is done pan-cancer and with 4 cohorts (melanoma, lung, bladder, other).
* Testing is done across a variety of covariate values (e.g. tumor location, age, purity). 
* Multiple testing is adjusted using the Benjamini–Yekutieli FDR procedure.
* Compute correlations of all features with five factors (TMB, T-cell, Pretreated, TGFB, Proliferation).

## Organization
* 0_filters filters features based on a variety of criteria (i.e. number non-zero, variability, mean value). 
* 1_prep prepares data structures for running exhaustive analysis.
* 2_run-far takes as inputs the outcome, cohort, and covariates and run's exhaustive testing.
* 3_rbind combine results from all combinations of outcomes, cohort, covariates.
* 4_get-cors computes correlations of all features to 5 latent factors.
* 5_combine combines results of testing and computed correlations.
* 6_*, 7_*, 8_* curate the results for display in figures (better names, groupings, etc.). 

## Requirements to run

### Pre-processed Hartwig Medical Foundation Data
* Study relies on pre-processed data from: https://github.com/bbglab/hartwig_biomarkers
* Requires data access to Hartwig data: https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/ 
* Hartwig data used for analyses was from data requst shared on August 9, 2022.
