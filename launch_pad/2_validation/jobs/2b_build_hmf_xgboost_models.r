wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/validation_help.R"))
library(tidyverse)
library(xgboost)

validation_ready <- readRDS( paste0(TMP_DIR, "validation-hmf-go.Rds") ) 

hmf <- (
    validation_ready 
        %>% filter(Study == "HMF-CPCT") 
        %>% mutate_at(vars("tcell","tgfb","prolif","pdl1", "pretreat", "pretreat_comp", "purity"), scale) 
        %>% mutate(tmb_bin = ifelse(exp(tmb)-1 > 10, 1, 0), model_apply = tissue)
        %>% drop_na(os, bor)
)

builder <- function( df, model_features ){
    files <- list()
    files[['X']] =  as.matrix( df %>% select(all_of(model_features)))
    files[['id']] = df %>% pull(patient_id)
    files[['complete_id']] = df %>% drop_na(tcell) %>% pull(patient_id)
    files[['Y_lr']] = df %>% pull(bor)
    files[['Y_os']] = df %>% pull(os)
    files
}

best_fit <- function(Y, X, hyper_grid = parameter_grid, model, base_model = NULL){
    grid_fit(Y, X, hyper_grid, model, base_model)$best_model
}

set.seed(62220)
K <- 1

mods <- list()

for( i in seq(K) ){
    
    hmf_models <- list()
    for( model in names(models)) {

        print(model); flush.console()

        #### Strore base data structure
        s <- list()
        s$all <- builder( hmf, models[[model]] )
        for (j in unique( hmf %>% pull(model_apply))) { 
            s[[j]] <- builder( hmf %>% filter(model_apply == j), models[[model]] )
        }
        #### Fit overall models 
        s$all$mod_lr <- best_fit( Y = s$all$Y_lr, X = s$all$X, model = "lr", base_model = NULL ) 
        s$all$mod_os <- best_fit( Y = s$all$Y_os, X = s$all$X, model = "os", base_model = NULL )

        #### Tissue specific models guided by overall model
        for (k in unique(hmf$tissue)) {
            s[[k]][['mod_lr']] <- best_fit( s[[k]]$Y_lr, s[[k]]$X, model = "lr", base_model = s$all$mod_lr)
            s[[k]][['mod_os']] <- best_fit( s[[k]]$Y_os, s[[k]]$X, model = "os", base_model = s$all$mod_os) 
            }        
        ### save model
        hmf_models[[model]] <- s
    }    
    ### save as repition model i
    mods[[as.character(i)]] <- hmf_models
}

saveRDS( mods, paste0(TMP_DIR, "validation-hmf-models.Rds"))
