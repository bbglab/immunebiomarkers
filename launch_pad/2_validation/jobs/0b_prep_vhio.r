wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))

library(tidyverse)

CLIN_DIR <- paste0(E_DIR, "/val_hebron/clinical/clean/")
RNA_DIR <- paste0(E_DIR, "/val_hebron/rna/clean/")
TMB_DIR <- paste0(E_DIR, "/val_hebron/somatic/clean/")
ANNOT_DIR <- paste0(E_DIR, "/val_hebron/rna/raw/")
O_DIR <- paste0(E_DIR, "/val_hebron/all_ready/")

clinical <- readRDS( paste0(CLIN_DIR, "clinical_go.Rds"))
somatic <-  readRDS( paste0(TMB_DIR,  "tmb_go.Rds"))
rna <-      readRDS( paste0(RNA_DIR,  "rna_go.Rds"))
annotations <- read.csv( paste0(ANNOT_DIR, "annotations_share.csv"), sep = ";", stringsAsFactors = FALSE)

vhio_data <- 
(
clinical 
    %>% inner_join( rna, by = "patient_id")
    %>% left_join( annotations, by = "LAB.ID")
    %>% left_join( somatic, by = "LAB.ID")
)

vhio_select <- (
    vhio_data 
        %>% transmute( 
            patient_id, 
            bor,
            os = ifelse( os_event == 0, -os_days, os_days),
            os_event,
            os_days,
            age = clinical_age, 
            gender = clinical_gender, 
            tissue = Cohort.x,
            tissue_full = Cohort.x, 
            tmb, 
            tcell,
            prolif, 
            tgfb,
            pdl1, 
            pretreat = clinical_pretreat,
            pretreat_comp = NA,
            purity = NA,
            Study = "VHIO"
        ) %>% drop_na(os_days)
          %>% group_by(tissue)
          %>% mutate(tmb=ifelse(is.na(tmb),median(tmb,na.rm=TRUE),tmb))  ### add median imputation
          %>% ungroup()
)

saveRDS( vhio_select, paste0(TMP_DIR, "validation-vhio-go.Rds"))
