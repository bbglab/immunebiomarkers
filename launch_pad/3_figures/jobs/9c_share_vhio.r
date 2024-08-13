wd <- dirname(dirname(getwd()))
setwd(wd)

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

source(paste0(wd,"/mission_control/treasure_map.R"))

VHIO_DATA <- "/workspace/projects/immune_biomarkers/repo/immune_biomarkers/external_data/val_hebron/"

clinical_ready <- readRDS(paste0(VHIO_DATA, "clinical/clean/clinical_ready.Rds")) 

extra_clinical <- clinical_ready %>% transmute(patient_id, clinical_age, clinical_sex = clinical_gender, clinical_biopsy_location, clinical_tumor_location, clinical_mechanism, clinical_recist)

mini_share <- readRDS(paste0(VHIO_DATA, "all_ready/vhio_ready.Rds")) %>% select(-pdl1, -LAB.ID)

rna_full <- 
readRDS(paste0(VHIO_DATA, "rna/clean/rna_full.Rds")) %>% 
  select(-LAB.ID, -Cohort, -BIOPSY, -Biopsy_date, -Comments)

full_share <- 
mini_share %>% 
  select(-tcell, -tgfb, -prolif) %>% 
  left_join(extra_clinical, by = "patient_id") %>% 
  left_join(rna_full, by = "patient_id")

fwrite( full_share, paste0(VHIO_DATA, "share/vhio_full_share.csv" ))
fwrite( mini_share, paste0(VHIO_DATA, "share/vhio_mini_share.csv" ))
