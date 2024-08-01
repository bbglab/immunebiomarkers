wd <- dirname(dirname(getwd()))

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
vhio_keep <- c(unique(paste("isofox_",unname(unlist(readRDS(paste0(REF_DIR, "vhio_gene_sets.Rds")))))),'isofox_gene_set_vhio_prolif')

library(tidyverse)
library(data.table)

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

hmf_markers <- fread(paste0(I_DIR, "signals_base.csv")) 

cpi <- (
    hmf_markers
       %>% filter(
           clinical_post_contains_Immunotherapy == TRUE,
           clinical_cpi_mechanism != "CTLA4",
           Survival_time_to_os_event != 0
        )
)
non_cpi <- (
    hmf_markers 
       %>% filter(
          clinical_post_contains_Immunotherapy == FALSE,  
          Survival_time_to_os_event != 0 
       )
)

dim(non_cpi)

feature_types <- unname(sapply( colnames(cpi), function(i) strsplit(i, "_")[[1]][1]))

drivers <- cpi %>% select(contains("driver_"))
drivers_filter <- colnames(drivers)[which(apply( drivers, 2, sum, na.rm = TRUE) < nrow(cpi)/30)]

sigs <- cpi %>% select(contains("sig_"))
sigs_filter <- colnames(sigs)[which(apply( sigs > .02, 2, sum) < nrow(cpi)/20)]

somatic.sig_filter <- gsub('sig_', 'somatic_TMB_', sigs_filter)

isofox_nr_filter <- colnames(cpi %>% select(contains("isofox.nr")))
gene_sets_keep <- colnames(cpi %>% select(contains("isofox_gene_set")))

isofox <- cpi %>% select(contains("isofox_"))
isofox_sd <- apply( isofox, 2, sd, na.rm = TRUE)
isofox_mn <- apply( isofox, 2, mean, na.rm = TRUE)
isofox_filter <- colnames(isofox)[-intersect(which(isofox_sd > .5), which(abs(isofox_mn) > .5))]

somatic.gene <- cpi %>% select(contains("somatic.gene"))
somatic.gene_filter <- colnames(somatic.gene)[apply( somatic.gene > 0, 2, sum) < nrow(cpi)/20]

filters <- unique(c(drivers_filter, sigs_filter, somatic.sig_filter, isofox_nr_filter, isofox_filter, somatic.gene_filter))
filters <- filters[-c(which( filters %in% vhio_keep), which(filters %in% gene_sets_keep))]

cpi_go <- (
    cpi %>% select(-all_of(filters)) 
        %>% mutate(
                tmb = somatic_summary_tmbPerMb, 
                tcell = isofox_gene_set_tcell_cluster, 
                prolif = isofox_gene_set_prolif_cluster,
                tgfb = isofox_gene_set_tgfb_cluster, 
                pretreat = (clinical_meta_hasSystemicPreTreatment2 + clinical_pre_treated)/2,
                pretreat_comp = clinical_systemic_composite,
                tissue = clinical_tumor_location_group,           ### covariate
                tissue_full = clinical_meta_primaryTumorLocation,
                age = clinical_age_at_treatment_start,            ### covariate
                biopsy = clinical_biopsy_site,                    ### covariate
                purity = somatic_summary_purity,                  ### covariate 
                pdl1 = isofox_CD274
        )
)
colnames(cpi_go) <- gsub( "-mutations", ".mb", colnames(cpi_go))

non_cpi_go <- (
    non_cpi
        %>% transmute(
                patient_id = patientIdentifier, 
                Y_best_response_binary, 
                Survival_os_event, 
                Survival_time_to_os_event, 
                Filter_meta_responseMeasured,
                tmb = somatic_summary_tmbPerMb, 
                tcell = isofox_gene_set_tcell_cluster, 
                prolif = isofox_gene_set_prolif_cluster,
                tgfb = isofox_gene_set_tgfb_cluster, 
                pretreat = (clinical_meta_hasSystemicPreTreatment2 + clinical_pre_treated)/2,
                pretreat_comp = clinical_systemic_composite,
                tissue = clinical_tumor_location_group,           ### covariate
                tissue_full = clinical_meta_primaryTumorLocation,
                age = clinical_age_at_treatment_start,            ### covariate
                biopsy = clinical_biopsy_site,                    ### covariate
                purity = somatic_summary_purity,                  ### covariate 
                pdl1 = isofox_CD274
        )
)

print("I ran")

saveRDS( cpi_go, paste0(I_DIR, "cpi_go.Rds"))
saveRDS( cpi_go %>% select( sampleId, contains("clinical"), somatic_summary_tmbPerMb, isofox_CXCL9), paste0(I_DIR, "cpi_check.Rds"))
saveRDS( names(cpi_go[,-seq(15)]), paste0(TMP_DIR, "exhaustive-features-go.Rds"))
saveRDS( non_cpi_go, paste0(TMP_DIR, "validation-hmf-noncpi-go.Rds"))
