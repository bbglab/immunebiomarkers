# Immune Biomarkers

## Goal
* Share code used to produce results for "Five factors underlie response to immunotherapy" paper. 

## Organization
* 0_analyses contains all analyses.  
* 1_figures produces paper figures (0_analyses run first).
* mission_control specifies directory locations and contains helper files.
* launch_pad contains scripts run all the analyses.
* env contains the conda environments .yaml files used to analyses and figures. 

## Requirements to run

### Pre-processed Hartwig Medical Foundation Data
* Main analyses relie on output: https://github.com/bbglab/hartwig_biomarkers
* Requires data access approval: https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/ 
* Hartwig data used for analyses was from data requst shared on August 9, 2022.

### Downloaded External Data
* For validation we downloaded external public data from other studies.
* Need to download these external data and specify directory locations. 

### QMap and conda environments
* Jobs are submitted with the QMap tool from the bbglab: https://github.com/bbglab/qmap
* QMap tool calls conda environments stored in env folder.
 
## Run it
* From this root directory:
```
$ python run.py
```
