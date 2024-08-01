wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/validation_help.R"))

library(tidyverse)
library(xgboost)

validation_ready <- readRDS(paste0(TMP_DIR, "validation-hmf-go.Rds"))
hmf_models <- readRDS(paste0(TMP_DIR, "validation-hmf-models.Rds"))

hmf_non_cpi <- (
    readRDS(paste0(TMP_DIR, "validation-hmf-noncpi-go.Rds")) 
        %>% drop_na(tcell)
        %>% select(-Y_best_response_binary, -Survival_os_event, -Survival_time_to_os_event, -Filter_meta_responseMeasured)
        %>% mutate( model_apply = ifelse(tissue!="other",tissue, "all")) 
        %>% mutate_at(vars("tcell","tgfb","prolif", "pretreat", "purity"), scale)
        %>% mutate(
              tcell = as.numeric(tcell),
              prolif = as.numeric(prolif), 
              pretreat = as.numeric(pretreat),
              tgfb = as.numeric(tgfb), 
              pdl1 = as.numeric(pdl1), 
              purity = as.numeric(purity)#,pretreat_comp = as.numeric(pretreat_comp)
        )
)

apply_hmf_mods <- function(df, model, features){
    model_features <- models[[features]]
    mod_lr <- hmf_models[['1']][[features]][[model]]["mod_lr"]$mod_lr
    mod_os <- hmf_models[['1']][[features]][[model]]["mod_os"]$mod_os
    if (nrow(df) > 0) {
        X <- as.matrix(df %>% select(all_of(model_features)))
        get_preds_X(X, df %>% pull(patient_id), mod_lr, mod_os)
    }
}

non_cpi_hmf_preds <- data.frame()
for (j in c("all", "bladder", "lung", "skin")){
    df_go <- hmf_non_cpi %>% filter(model_apply == j)
    preds_i <- apply_hmf_mods( df_go, model = j, features = "five_latent_purity")
    preds_i$model <- "five_latent_purity"
    non_cpi_hmf_preds <- bind_rows(non_cpi_hmf_preds, preds_i)
}

non_cpi <- hmf_non_cpi %>% inner_join(non_cpi_hmf_preds, by = "patient_id") 
non_cpi$Study <- 'HMF'
non_cpi$cpi <- FALSE

non_cpi$pred_os2 <- ifelse( non_cpi$pred_os > 3, 3.1, non_cpi$pred_os)
non_cpi$lr_gp <-  cut( non_cpi$pred_lr,  breaks = c(0,.1,.5,    1) , labels = c("Low", "Medium", "High"))
non_cpi$os_gp <-  cut( non_cpi$pred_os2,  breaks = c(0,.5,1.5,   20) , labels = c("Low", "Medium", "High"))

saveRDS( non_cpi, paste0(TMP_DIR, "validation-non-cpi.Rds") )
