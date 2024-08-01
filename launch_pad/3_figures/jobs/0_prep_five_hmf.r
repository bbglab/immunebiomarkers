wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

#cpi %>% select(contains("pretreat_comp"))

tmb_features <- c("somatic_TMB", "somatic_TMB_clonal", "somatic_TMB_indel", "somatic_TMB_frameshift", "somatic_TMB_subclonal")
tcell_features <- c("isofox_gene_set_t_cell_effector", "isofox_gene_set_t_cell_gep_18", "isofox_gene_set_cyt", "isofox_CD8A", "isofox_CXCL9")
tgfb_features <- c("isofox_gene_set_Pan_TBRS", "isofox_gene_set_F_TBRS", "isofox_gene_set_mariathan_EMT2", "isofox_TGFB1","isofox_COL4A1")
prolif_features <- c("isofox_gene_set_prolif", "isofox_gene_set_mariathan_Cell_cycle", "isofox_gene_set_mariathan_DNA_replication",  "isofox_TOP2A", "isofox_CDK1")
cor_features <- c(tmb_features, tcell_features, tgfb_features, prolif_features)

#cor_features

ready <- (
    cpi %>% 
    transmute(
        bor = Y_best_response_binary,
        bor_measured = Filter_meta_responseMeasured,
        pfs_event = Survival_pfs_event, 
        pfs_days = Survival_time_to_pfs_event, 
        os_event = Survival_os_event, 
        os_days = Survival_time_to_os_event,
        pan = "Pan-Cancer", 
        tissue, 
        biopsy,
        purity, 
        tmb, 
        tcell, 
        prolif, 
        tgfb,
        pretreat,
        
        ### Add pretreatment features ### 
        pretreat_comp,
        time_since_last_treatment_years = ifelse( pretreat_comp == 3650, NA, pretreat_comp/365),
        overall_pretreat = clinical_pre_treated,
        radio = clinical_meta_hasRadiotherapyPreTreatment,
        chemo = clinical_pre_contains_Chemotherapy, 
        hormonal = clinical_pre_contains_Hormonal, 
        immuno = clinical_pre_contains_Immunotherapy,
        targeted = clinical_pre_contains_Targeted,
        
        ### nice features for correlation plots 
        TMB = somatic_TMB, 
        "TMB clonal" = somatic_TMB_clonal,
        "TMB Indels" = somatic_TMB_indel, 
        "TMB frameshift" = somatic_TMB_frameshift, 
        "TMB sub-clonal" = somatic_TMB_subclonal,
        
        "T-cell Effector" = isofox_gene_set_t_cell_effector, 
        "GEP 18" = isofox_gene_set_t_cell_gep_18, 
        "CYT" = isofox_gene_set_cyt, 
        "CD8A" = isofox_CD8A, 
        "CXCL9" = isofox_CXCL9,
        
        "Pan-TBRS" = isofox_gene_set_Pan_TBRS, 
        "F-TBRS" = isofox_gene_set_F_TBRS, 
        "EMT2" = isofox_gene_set_mariathan_EMT2, 
        "TGFB1" = isofox_TGFB1, 
        "COL4A1" = isofox_COL4A1,
      
        "Proliferation" = isofox_gene_set_prolif, 
        "Cell Cycle" = isofox_gene_set_mariathan_Cell_cycle, 
        "DNA Replication" = isofox_gene_set_mariathan_DNA_replication, 
        "TOP2A" = isofox_TOP2A, 
        "CDK1" = isofox_CDK1
    ) 
)

ready$bor <- ifelse( ready$bor_measured != "Yes", NA, ready$bor)
ready$response_group <- ifelse(ready$bor == 1, "Responders", "Non-Responders")
ready$tissue <- str_to_title(ready$tissue)
ready$tissue <- factor(ready$tissue, levels = c("Skin","Lung","Bladder","Other"))
ready$biopsy <- str_to_title(ready$biopsy)

tcell_tertiles <- quantile(ready$tcell, probs = c(.33,.66), na.rm = TRUE)
ready$pretreat2 <- factor(ifelse(ready$pretreat == 1, "Systemic Pretreatment", "No Systemic Pretreatment"), levels = c("No Systemic Pretreatment", "Systemic Pretreatment"))
ready$tcell2 <- factor( ifelse( ready$tcell > median(ready$tcell, na.rm = TRUE), "T-cell High", "T-cell Low"), levels = c("T-cell Low", "T-cell High"))
ready$tcell3 <- cut(ready$tcell, c(0, tcell_tertiles[1],tcell_tertiles[2], 1000), labels = c("T-cell Low", "T-cell Medium", "T-cell High"))
ready$tgfb2 <- factor( ifelse( ready$tgfb > median(ready$tgfb, na.rm = TRUE), "TGFB High", "TGFB Low"), levels = c("TGFB Low", "TGFB High"))
ready$prolif2 <- factor( ifelse(ready$prolif > median(ready$prolif, na.rm = TRUE), "Proliferation High", "Proliferation Low"), levels = c("Proliferation Low", "Proliferation High"))
ready$tmb2 <- factor( ifelse(exp(ready$tmb)> 11, "TMB High", "TMB Low"), levels = c("TMB Low", "TMB High"))
ready$purity2 <- factor( ifelse(ready$purity > median(ready$purity), "Purity High", "Purity Low"), levels = c("Purity High", "Purity Low"))

ready$multi <- apply(ready %>% select(chemo, hormonal, immuno, targeted), 1, sum)

saveRDS( ready, paste0( TMP_DIR, 'supplement-five.Rds'))
saveRDS( ready %>% select(bor, pfs_event, os_event, tissue, biopsy, tcell), paste0( FIG_FINAL_DIR, '0_prep_five_hmf_summary_table.Rds'))
