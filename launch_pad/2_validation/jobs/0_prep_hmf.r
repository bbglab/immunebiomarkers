wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))
library(tidyverse)

all <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

all$os <- ifelse( all$Survival_os_event == 0, 
                 -all$Survival_time_to_os_event, 
                  all$Survival_time_to_os_event)

hmf_go <- (
    all 
      %>% transmute(
        patient_id = patientIdentifier,
        bor =  Y_best_response_binary,
        os,
        os_event = Survival_os_event,
        os_days = Survival_time_to_os_event,
        age,
        gender = clinical_meta_gender,
        tissue,
        tissue_full = clinical_meta_primaryTumorLocation,
        tmb = somatic_summary_tmbPerMb,
          
        tcell, 
        tgfb, 
        prolif,
          
        tcell_cluster5 = isofox_gene_set_tcell_cluster_05,
        tgfb_cluster5 = isofox_gene_set_tgfb_cluster_05,
        prolif_cluster5 = isofox_gene_set_prolif_cluster_05,
          
        tcell_set = isofox_gene_set_t_cell_effector,
        prolif_set = isofox_gene_set_prolif,
        tgfb_set = isofox_gene_set_Pan_TBRS, 
          
        pdl1 = isofox_CD274,
        pretreat,
        pretreat_comp = clinical_systemic_composite, 
        purity,
        Study = "HMF-CPCT"
    ) 
)

saveRDS( hmf_go, paste0(TMP_DIR, "validation-hmf-go.Rds") )
