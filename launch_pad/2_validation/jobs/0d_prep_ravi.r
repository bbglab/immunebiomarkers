options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))

library(tidyverse)
library(stringr)

hmf_tmb_ref <- readRDS(paste0(I_DIR, "cpi_go.Rds")) %>% filter( clinical_tumor_location_group == "lung") %>% transmute(bor = Y_best_response_binary, tmb = somatic_summary_tmbPerMb)

I_DIR <- paste0(E_DIR, "ravi_lung/")
clinical <- read.csv( paste0(I_DIR, "Table_S1_Clinical_Annotations.csv"), sep = ";")
somatic <- read.csv( paste0(I_DIR, "Table_S5_Mutation_Burden.csv"), sep = ";")
rna <- read.csv( paste0(I_DIR, "Table_S13_RNA_TPM.csv"), sep = ";")

binary_response <- function( recist ){
    recist = as.character(recist)
    if( recist %in% c("CR","PR")){
        1
    } else if ( recist %in% c("SD", "PD")){
        0
    } else {
        NA
    }
}

clinical$bor <- unlist(lapply( as.character(clinical$"Harmonized_Confirmed_BOR"), binary_response))

clinical_go <- (
    clinical
        %>% filter(Pre.treatment_RNA_Sample_QC != "Flag", Agent_PD1 != "Atezolizumab")
        %>% transmute(
          patient_id = Harmonized_SU2C_Participant_ID_v2, 
          bor,
          os = ifelse(Harmonized_OS_Event == 1, Harmonized_OS_Days, -Harmonized_OS_Days), 
          os_event = Harmonized_OS_Event, 
          os_days = Harmonized_OS_Days, 
          age = Patient_Age_at_Diagnosis, 
          gender = Patient_Sex,
          tissue = "lung", 
          tissue_full = Institution,
          pretreat = ifelse(Line_of_Therapy > 1, 1, 0),
          pretreat_comp = NA,
          purity = NA,
          Study = "RAVI",
          extra = Initial_Stage, 
          extra2 = Pre.treatment_RNA_Sample_QC
        )
)

hmf_quants <- data.frame( vals = quantile(hmf_tmb_ref$tmb, probs = seq(0,1,.01))) %>% rownames_to_column("quantile")

somatic$patient_id <- unlist(lapply( somatic$Harmonized_SU2C_WES_Tumor_Sample_ID_v2, function(i) str_split(i, "-T1")[[1]][1] ))

somatic <- somatic %>% arrange(TMB)

hmf_quants <- (
    data.frame( vals = quantile(hmf_tmb_ref$tmb, probs = seq(0,1,.01))) 
        %>% rownames_to_column("quantile")
)   
ravi_quants <- (
    data.frame( 
        patient_id = somatic$patient_id, 
        tmb = log(somatic$TMB+1), 
        vals = quantile(somatic$TMB, probs = seq(0,1,1/308))
    ) %>% rownames_to_column("quantile")
)

round_quants <- function(i) floor(as.numeric(as.character(strsplit(i, "%")[[1]][1])))
hmf_quants$round_quantile <- unlist(lapply(hmf_quants$quantile, round_quants))
ravi_quants$round_quantile <- unlist(lapply(ravi_quants$quantile, round_quants))

somatic_go <- (
    ravi_quants 
        %>% select(patient_id, round_quantile, tmb) 
        %>% left_join( hmf_quants , by = "round_quantile")
        %>% transmute(patient_id, tmb = vals)
)

step1 <- rna %>% select(-Name)
step2 <- step1 %>% group_by(Description) %>% summarise(ct = n()) %>% filter(ct == 1) %>% pull(Description)
step3 <- step1 %>% filter(Description %in% step2) %>% column_to_rownames("Description")
step4 <- data.frame(t(step3))

numerizer <- function(ll) as.numeric(gsub(",",".", as.character(ll)))
share <- data.frame(lapply(step4, numerizer))

saveRDS( share, paste0(REF_DIR, "rna_ravi.Rds") )

genes <- unlist(gene_sets)
twist <- data.frame(t(
    column_to_rownames(
        rna 
            %>% filter( Description %in% genes ) 
            %>% select(-Name), 
        "Description")
))

numerizer <- function(ll) as.numeric(gsub(",",".", as.character(ll)))
clean <- data.frame(lapply(twist, numerizer)) %>% mutate_all(~(log(.+1) %>% as.vector))

patient_ids <- unlist(lapply( rownames(twist), function(i) gsub("[.]", "-", as.character(str_split(i, ".T1")[[1]][1] ))))

rna_clean <- data.frame( patient_id = patient_ids, clean)
rna_clean$patient_id <- as.character(rna_clean$patient_id)

rna_clean$tcell  <- apply(rna_clean %>% select( any_of(gene_sets$clusters$tcell)),1,mean)
rna_clean$prolif <- apply(rna_clean %>% select( any_of(gene_sets$clusters$prolif)),1,mean)
rna_clean$tgfb   <- apply(rna_clean %>% select( any_of(gene_sets$clusters$tgfb)),1,mean)

rna_clean$tcell_cluster5  <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$tcell)),1,mean)
rna_clean$prolif_cluster5 <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$prolif)),1,mean)
rna_clean$tgfb_cluster5   <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$tgfb)),1,mean)

rna_clean$tcell_set  <- apply(rna_clean %>% select( any_of(gene_sets$sets1$tcell)),1,mean)
rna_clean$prolif_set <- apply(rna_clean %>% select( any_of(gene_sets$sets1$prolif)),1,mean)
rna_clean$tgfb_set   <- apply(rna_clean %>% select( any_of(gene_sets$sets1$tgfb)),1,mean)

rna_clean$pdl1   <- apply(rna_clean %>% select( CD274 ),1,mean)

rna_go <- rna_clean %>% select(patient_id, tcell, prolif, tgfb, tcell_cluster5, prolif_cluster5, tgfb_cluster5, tcell_set, prolif_set, tgfb_set, pdl1)

ravi_go <- clinical_go %>% left_join( somatic_go, by = "patient_id") %>% left_join( rna_go, by = "patient_id")

saveRDS( ravi_go, paste0( TMP_DIR, "validation-ravi-go.Rds"))
