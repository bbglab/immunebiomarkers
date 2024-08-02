# Prediction Study

## Goal
* Evaluate XGBoost models of OS and best overall response.

## Methodology 
* Models fit across grid of tuning parameters (sub-sampling, eta, depth)
* Trained pan-cancer models used as initial inputs to cohort based models (Skin, Lung, Bladder, Other).
* Models highlighted in paper have depth = 1 (i.e. no interaction), see similar performance.
* Evaluate with cross-validation on 1000 repeated samples of 80/20 training/test data splits.  
* Run models with different feature combinations
* For more details see paper supplement Section 5.

## Organization
* 0_prep prepares data for evaluation study.
* 1_study runs model evaluations for specific outcome (e.g. response) and features.
* 2_combine combines results from 1_study across all outcome and feature settings.
* Paper settings in ~/mission_control/paper_settings.R
* Logic for fitting and evaluation ~/mission_control/eval_help.R  

## Requirements to run

### Pre-processed Hartwig Medical Foundation Data
* Main analyses relie on output: https://github.com/bbglab/hartwig_biomarkers
* Requires data access approval: https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/ 
* Hartwig data used for analyses was from data requst shared on August 9, 2022.
