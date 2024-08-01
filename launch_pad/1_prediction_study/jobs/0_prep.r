wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(fastDummies)
library(tidyverse)

all <- readRDS(paste0(I_DIR, "cpi_go.Rds")) %>% mutate(tmb_bin = as.numeric(exp(somatic_summary_tmbPerMb)-1 > 10)) 

response_tissue <- c(
                     "Y_best_response_binary", 
                     "Survival_pfs_event",
                     "Survival_time_to_pfs_event", 
                     "Survival_os_event",
                     "Survival_time_to_os_event", 
                     "pfs", 
                     "os", 
                     "tissue",
                     "Filter_meta_responseMeasured"
                    )

clinical_features <- c( "pretreat",
                        "pretreat_comp",
                        "clinical_pre_treated",
                        "clinical_meta_hasSystemicPreTreatment2",
                        "age",
                        "clinical_biopsy_distal_proximal",
                        "clinical_cpi_mechanism3",
                        "hla_lilac_del_hla",
                        "cnv_summary_wholeGenomeDuplication",
                        "sv_summary_svTumorMutationalBurden", 
                        "purity"
                    )

genomic_features <- c( "tmb_bin", "tmb", "tcell", "prolif", "tgfb", "pdl1" )

latent_features <- c(
       "somatic_TMB_clonal",
       "somatic_TMB_vhio",
       "somatic_TMB_exome",
       "isofox_gene_set_t_cell_effector", 
       "isofox_gene_set_prolif", 
       "isofox_gene_set_Pan_TBRS",
       "isofox_gene_set_t_cell_gep_18",
       "isofox_gene_set_mariathan_Cell_cycle", 
       "isofox_gene_set_mariathan_EMT2",
       "isofox_gene_set_vhio_tgfb",
       "isofox_gene_set_vhio_prolif", 
       "isofox_gene_set_vhio_tcell"
)

all$pfs <- ifelse( all$Survival_pfs_event == 0,   
                  -all$Survival_time_to_pfs_event,     
                   all$Survival_time_to_pfs_event)
all$os <- ifelse(  all$Survival_os_event == 0, 
                  -all$Survival_time_to_os_event, 
                   all$Survival_time_to_os_event)

mini <- (
    all %>% column_to_rownames('sampleId')
        %>% select( all_of(response_tissue),   all_of(genomic_features), 
                    all_of(clinical_features), all_of(latent_features) )
        %>% rename( response = Y_best_response_binary)
)

builder <- function( mini ){

    xg_lr <- ( mini 
              %>% select(-pfs,-os, 
                         -Survival_pfs_event, -Survival_time_to_pfs_event, 
                         -Survival_os_event, -Survival_time_to_os_event)
              %>% filter(Filter_meta_responseMeasured == "Yes")
              %>% select(-Filter_meta_responseMeasured)
              %>% drop_na(response))

    xg_pfs <- ( mini 
              %>% select(-response, -os, -Survival_os_event, -Survival_time_to_os_event) 
              %>% rename( event_status = Survival_pfs_event, 
                          time_to_event = Survival_time_to_pfs_event)
              %>% select(-Filter_meta_responseMeasured) 
              %>% drop_na(pfs))

    xg_os <- ( mini 
              %>% select(-response, -pfs, -Survival_pfs_event, -Survival_time_to_pfs_event) 
              %>% rename( event_status = Survival_os_event, 
                          time_to_event = Survival_time_to_os_event)
              %>% select(-Filter_meta_responseMeasured)
              %>% drop_na(os))
    
    ### storage closet ### 
    eval_closet <- list()
    eval_closet[['lr']][['df']] <- xg_lr
    eval_closet[['lr']][['gps']] <- apply(data.frame( xg_lr$response, 
                                                      xg_lr$tissue, 
                                                      is.na(xg_lr$tcell)), 
                                            1, function(i) paste0(i[1],"-",i[2],"-",i[3]))

    eval_closet[['pfs']][['df']] <- xg_pfs
    eval_closet[['pfs']][['gps']] <- apply(data.frame(xg_pfs$event_status, 
                                                            xg_pfs$tissue, 
                                                      is.na(xg_pfs$tcell)), 
                                            1, function(i) paste0(i[1],"-",i[2],"-",i[3]))                                   
                                            
    eval_closet[['os']][['df']] <- xg_os
    eval_closet[['os']][['gps']] <- apply(data.frame(xg_os$event_status, xg_os$tissue, is.na(xg_os$tcell)), 
                                            1, function(i) paste0(i[1],"-",i[2],"-",i[3]))                                                                           
    
    eval_closet
}

saveRDS( builder(mini), paste0(TMP_DIR, "xg-eval-prep.Rds"))
