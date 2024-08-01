wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/validation_help.R"))
library(tidyverse)

validation_ready <- readRDS( paste0( TMP_DIR, "validation-go.Rds"  ))
hmf_models <- readRDS( paste0( TMP_DIR, "validation-hmf-models.Rds"))

#head(validation_ready)

apply_hmf_mods <- function (df, build, model, features) 
{
    model_features <- models[[features]]
    mod_lr <- hmf_models[[build]][[features]][[model]]["mod_lr"]$mod_lr
    mod_os <- hmf_models[[build]][[features]][[model]]["mod_os"]$mod_os
    if (nrow(df) > 0) {
        X <- as.matrix(df %>% select(all_of(model_features)))
        get_preds_X(X, df %>% pull(patient_id), mod_lr, mod_os)
    }
}

models <- list( 
    "tmb_bin" = c("tmb_bin"),
    "base_bin" = c("tmb_bin", "pdl1"),
    "tmb" = c("tmb"),
    "base" = c("tmb", "pdl1"),
    "rna" = c("tcell", "prolif", "tgfb"),
    #"no_pretreat" = c("tcell", "prolif", "tgfb", "tmb"),
    "no_tmb" = c("pretreat", "tcell", "prolif", "tgfb"),
    "five_latent" = c("tmb", "tcell", "prolif", "tgfb", "pretreat"),
    "five_latent_purity" = c("tmb", "tcell", "prolif", "tgfb", "pretreat", "purity")
)

preds <- data.frame()
for( i in names(hmf_models)){
    for (j in c("all", "bladder", "lung", "skin", "other")){
        print(j); #flush.console()
        for( k in names(models)){
            if( k != 'five_latent_hmf' ){
                df_go <- validation_ready %>% filter(model_apply == j)
                preds_i <- apply_hmf_mods( df_go, build = i, model = j, features = k )
                preds_i$model <- k
                preds_i$build <- i
                preds <- bind_rows(preds, preds_i)
            }
        }
    }
}
preds <- preds %>% unique() ### remove duplicates

preds_mean <- preds %>% group_by(patient_id, model) %>% summarise_all( ~mean(.x, na.rm = TRUE)) %>% ungroup()

validation_go <- (
    validation_ready 
        %>% inner_join(preds_mean, by = "patient_id")
        %>% mutate(
            tcell = as.numeric(tcell),
            prolif = as.numeric(prolif),
            tgfb = as.numeric(tgfb), 
            pdl1 = as.numeric(pdl1), 
            pretreat = as.numeric(pretreat),
            purity = as.numeric(purity),
            pretreat_comp = as.numeric(pretreat_comp),
            age = as.numeric(age))
)

validation_go$pred_os2 <- ifelse( validation_go$pred_os > 3, 3.1, validation_go$pred_os)
validation_go$lr_gp <-  cut( validation_go$pred_lr,  breaks = c(0,.1,.5,    1) , labels = c("Low", "Medium", "High"))
validation_go$os_gp <-  cut( validation_go$pred_os2,  breaks = c(0,.5,1.5,   20) , labels = c("Low", "Medium", "High"))

saveRDS(validation_go, file = paste0(TMP_DIR,"validation-measure.Rds"))
