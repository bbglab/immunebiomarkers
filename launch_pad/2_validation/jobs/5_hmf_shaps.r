wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/shapley_help.R"))

library(tidyverse)
library(xgboost)
library(survival)
library(survcomp)

cpi <- (
    readRDS(paste0(TMP_DIR, "validation-hmf-go.Rds")) 
        %>% drop_na(os, bor)
        %>% mutate_at(vars("tmb", "tcell","tgfb","prolif","pretreat", "purity"), scale) 
)

model_features <- c("tmb", "tcell", "prolif", "tgfb", "pretreat", "purity")

builder <- function( df ){
    files <- list()
    files[['X']] =  as.matrix( df %>% select(all_of(model_features)))
    files[['id']] = df %>% select(patient_id, age, gender, bor)
    files[['complete_id']] = df %>% drop_na(tcell) %>% pull(patient_id)
    files[['Y_lr']] = df %>% pull(bor)
    files[['Y_os']] = df %>% pull(os)
    files
}

store <- list()
store[['all']] <- builder( cpi )
for (i in unique(cpi %>% pull(tissue))) store[[i]] <- builder( cpi %>% filter(tissue == i) )

store[['all']][['mod_lr']] <- grid_fit( store[['all']]$Y_lr, 
                                        store[['all']]$X, 
                                        parameter_grid, 
                                        model = "lr" )$best_model 
store[['all']][['mod_os']] <- grid_fit( store[['all']]$Y_os, 
                                        store[['all']]$X, 
                                        parameter_grid, 
                                        model = "os" )$best_model

for (i in unique(cpi %>% pull(tissue))) {
    print(i)
    flush.console()
    store[[i]][['mod_lr']] <- grid_fit( store[[i]]$Y_lr, 
                                        store[[i]]$X, 
                                        parameter_grid, 
                                        model = "lr", 
                                        base_model = store[['all']][['mod_lr']])$best_model 
    store[[i]][['mod_os']] <- grid_fit( store[[i]]$Y_os, 
                                        store[[i]]$X, 
                                        parameter_grid, 
                                        model = "os",
                                        base_model = store[['all']][['mod_os']]
                                      )$best_model 
}

shap_binder <- function( i, model = "lr" ){

    goods <- store[[i]]
    X <- goods$X
    
    if( model == "lr"){
        mod <- goods$mod_lr
    } else {
        mod <- goods$mod_os
    }    
    
    ### extract Shapleys
    shaps <- data.frame( predict( mod , X, predcontrib = TRUE) )
    shaps$pred <- predict( mod , X, predcontrib = FALSE) 
    shaps$model <- i
    shaps$patient_id <- goods$id$patient_id
    shaps$age <- goods$id$age
    shaps$gender <- goods$id$gender
    shaps$response <- goods$id$response

    shaps %>% select(patient_id, model, everything())
}

for (i in names(store)) {
    print(i)
    flush.console()
    store[[i]][['shaps_lr']] <- shap_binder(i,"lr")
    store[[i]][['shaps_os']] <- shap_binder(i,"os")
}

name_map <- list( "tmb" = "TMB", 
                  "tcell" = "T-cell", 
                  "prolif" = "Proliferation", 
                  "tgfb" = "TGFB", 
                  "pretreat" = "Prior Systemic Therapy",
                  "purity" = "Purity")

chunker <- function( dataset, model, feature){
    
    goods <- store[[dataset]]
    X <- goods$X
    if( model == "lr"){
        shaps <- goods$shaps_lr
    } else {
        shaps <- goods$shaps_os
    }
    
    chunk <- data.frame( scaled_feature = X[,feature], shap_feature = shaps[,feature])
    chunk$feature <- name_map[[feature]]
    chunk$model <- model
    chunk$dataset <- dataset
    chunk$patient_id <- shaps$patient_id
    chunk$age <- shaps$age
    chunk$gender <- shaps$gender
    chunk$pred <- shaps$pred
    chunk$response <- shaps$response

    chunk %>% select( model, dataset, feature, everything())
}

p1 <- data.frame()
for ( model in c("os", "lr") ) {
    print(model)
    for ( dataset in names(store) ) {
        for ( feature in model_features)  {
            p1 <- rbind(p1, chunker(dataset, model, feature))
        }
    }
}

ready <- p1 %>% group_by( dataset ) %>% arrange(scaled_feature) %>% ungroup()
ready$dataset <- str_to_title(ready$dataset)
ready$dataset <- ifelse( ready$dataset == "All", "Pan-Cancer (base-model)", ready$dataset)
ready$dataset <- factor(
    ready$dataset,
    levels = c("Pan-Cancer (base-model)", "Skin", "Lung", "Bladder", "Other")
)
ready$feature <- factor(
    ready$feature,
    levels = c("T-cell", "TMB", "Prior Systemic Therapy", "TGFB", "Proliferation", "Purity")
)

print("worked")

saveRDS( ready, paste0(TMP_DIR, "validation-hmf-shaps.Rds") )
