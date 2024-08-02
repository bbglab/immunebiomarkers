# Validation

## Goal
* Share code used to run validation analyses.

## Methodology
* Download public data from non-Hartwig CPI PD/PDL1 studies.
* For each study, extract response and survival; TMB and pretreatment; compute T-cell, TGFB, Proliferation gene sets 
* Fit multivariate/univariate logistic regression and cox-ph regressions within/across studies (Figure 4, supplement)
* Apply XGBoost models build on Hartwig data to external studies (Supplement)

## Organization
* 0_* files each prepare the underyling data downloaded for each study (see requirements).   
* 1_combied_studies.ipynb combines the prepared datasets into one dataframe. 
* 2_* build Hartwig XGBoost models based on all CPI patients.
* 3_* get leave-one-out estimates for Hartwig XGBoost models and extract dependence plots.  
* 4_apply_hmf_model_external apply the Hartwig trained XGBoost models to other external studies.
* 5_hmf_shaps extracts shapley values from Hartwig XGBoost models.
* 6_hmf_non_cpi_apply applies XGBoost models to non-CPI treated patients in Hartwig. 
* 7_measure evaluates the predictive performance of XGBoost models on external studies. 

## Requirements to run

### Pre-processed Hartwig Medical Foundation Data
* Main analyses relie on output: https://github.com/bbglab/hartwig_biomarkers
* Requires data access approval: https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/ 
* Hartwig data used for analyses was from data requst shared on August 9, 2022.

### Downloaded External Data
* See methods section -> Validation Cohorts for description of external studies
* Need to download these external data to reproduce analyses. 
* INSPIRE: https://github.com/pughlab/inspire-genomics/.
* Lyon: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159067,
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161537,
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162519,
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162520.
* MARIATHASAN: http://research-pub.gene.com/IMvigor210CoreBiologies.
* PARKER ICI: https://github.com/ParkerICI/MORRISON-1-public.
* RAVI: https://zenodo.org/records/7625517.
* VHIO: We share pre-processed data in the ref directory of repository. Repo has code to show data was pre-processed. 
