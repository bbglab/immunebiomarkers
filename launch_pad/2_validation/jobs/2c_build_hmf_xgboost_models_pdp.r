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
    %>% mutate_at(vars("tcell","tgfb","prolif", "pretreat", "tmb"), scale)
    %>% mutate(model_apply = tissue)
    %>% drop_na(os, bor, tcell)
    %>% select("patient_id", "bor", "os", "tissue", "tcell","tgfb","prolif", "pretreat", "tmb")
)

best_fit <- function(Y, X, hyper_grid = parameter_grid, model, base_model = NULL){
    grid_fit(Y, X, hyper_grid, model, base_model)$best_model
}

set.seed(62220)
K <- 1

dfs <- list()
for( i in unique(hmf$tissue)){
    dfs[[i]] <- hmf %>% filter(tissue == i)
}

X <- as.matrix(hmf %>% select("tcell","tgfb","prolif", "pretreat", "tmb"))
pan_lr <- grid_fit(as.matrix(hmf$bor), X, parameter_grid, model = "lr")$best_model
pan_os <- grid_fit(as.matrix(hmf$os), X, parameter_grid, model = "os")$best_model

pred_lr <- as.data.frame(predict(pan_lr, X, predcontrib = TRUE)) %>% mutate(mod = "lr", tissue = "pan", patient_id = hmf$patient_id)
pred_os <- as.data.frame(predict(pan_os, X, predcontrib = TRUE)) %>% mutate(mod = "os", tissue = "pan", patient_id = hmf$patient_id)
pan_pred <- rbind(pred_lr, pred_os)

dfs <- list()
for( i in unique(hmf$tissue)){
    dfs[[i]] <- hmf %>% filter(tissue == i)
}

pred_maker <- function(i){
    df <- dfs[[i]]
    Y_bor <- as.matrix( df$bor )
    Y_os <- as.matrix(df$os)
    X <- as.matrix( df %>% select("tcell","tgfb","prolif", "pretreat", "tmb"))
    mod_lr <- grid_fit(Y_bor, X, parameter_grid, model = "lr", base_model = pan_lr)$best_model
    mod_os <- grid_fit(Y_os, X, parameter_grid, model = "os", base_model = pan_os)$best_model
    pred_lr <- as.data.frame(predict(mod_lr, X, predcontrib = TRUE)) %>% mutate(mod = "lr", tissue = i, patient_id = df$patient_id)
    pred_os <- as.data.frame(predict(mod_os, X, predcontrib = TRUE)) %>% mutate(mod = "os", tissue = i, patient_id = df$patient_id)
    rbind(pred_lr, pred_os) 
}

preds <- list()

for( i in unique(hmf$tissue)){
    print(i)
    preds[[i]] <- pred_maker(i)
}

vamonos <- rbind(do.call("rbind", preds), pan_pred)

saveRDS( vamonos, paste0(TMP_DIR, "validation-hmf-preds-pdp.Rds"))
