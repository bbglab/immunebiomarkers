### Functions to complete Cross-validation on HMF data  ###

extract_Y <- function( Y_base, id, i, j ) {
    if( j == "all" ) {
        Y_base %>% filter(patient_id != id) %>% pull(i)
    } else {
        Y_base %>% filter(patient_id != id, tissue == j) %>% pull(i)
    }
}
extract_X <- function( X_base, id, i ) { 
    if( i == "all"){
        as.matrix(X_base %>% filter(patient_id != id) %>% select(-patient_id, -tissue))
    } else {
        as.matrix(X_base %>% filter(patient_id != id, tissue == i) %>% select(-patient_id, -tissue))
    }
}
extract_X_i <- function( X_base, id ) {
    as.matrix(X_base %>% filter(patient_id == id) %>% select(-patient_id, -tissue))
}
get_tissue <- function( X_base, id ) X_base %>% filter(patient_id == id) %>% pull(tissue)

get_preds_X_i <- function( X_base, id, warm_lr, warm_os ){
    
    pred_lr <- predict(  warm_lr, extract_X_i( X_base, id ), predcontrib = FALSE)
    pred_os <- predict(  warm_os, extract_X_i( X_base, id ), predcontrib = FALSE)
    
    shaps_lr <- data.frame(t(predict( warm_lr, extract_X_i( X_base, id ), predcontrib = TRUE)))
    shaps_os <- data.frame(t(predict( warm_os, extract_X_i( X_base, id ), predcontrib = TRUE)))
    
    names(shaps_lr) <- paste0("shap_lr_", names(shaps_lr));
    names(shaps_os) <- paste0("shap_os_", names(shaps_os));
    
    data.frame( "patient_id" = id, "pred_lr" = pred_lr, "pred_os" = pred_os, shaps_lr, shaps_os)

}

run_loo_cv <- function( Y_base, X_base){
    loo_cv <- data.frame()

    for ( id in (X_base %>% pull(patient_id))){ 

        print(id)
        flush.console()

        ### Leave one-out extract ### 
        Y_all_lr <-  extract_Y( Y_base, id, "bor", "all" )
        Y_all_os <-  extract_Y( Y_base, id, "os", "all" )

        Y_tissue_lr <-  extract_Y(  Y_base, id, "bor", get_tissue(X_base, id) )
        Y_tissue_os <-  extract_Y(  Y_base, id, "os",  get_tissue(X_base, id) )

        X_all_i <-  extract_X( X_base, id, "all")
        X_tissue_i <- extract_X( X_base, id, get_tissue( X_base, id ))

        ### Build Models 
        base_lr <-  grid_fit( Y_all_lr,   X_all_i, parameter_grid, model = "lr" )$best_model
        base_os <-  grid_fit( Y_all_os,   X_all_i, parameter_grid, model = "os" )$best_model

        warm_lr <-  grid_fit( Y_tissue_lr,   X_tissue_i, parameter_grid, model = "lr", base_model = base_lr)$best_model
        warm_os <-  grid_fit( Y_tissue_os,   X_tissue_i, parameter_grid, model = "os", base_model = base_os)$best_model

        ### Apply predictions 
        loo_cv <- rbind(loo_cv, get_preds_X_i( X_base, id, warm_lr, warm_os ))
    }
    loo_cv
}


### Functions to apply to external datasets ###
models <- list( 
    "tmb_bin" = c("tmb_bin"),
    "base_bin" = c("tmb_bin", "pdl1"),
    "tmb" = c("tmb"),
    "base" = c("tmb", "pdl1"),
    "rna" = c("tcell", "prolif", "tgfb"),
    "no_pretreat" = c("tmb", "tcell", "tgfb", "prolif"), 
    "no_tmb" = c("pretreat", "tcell", "prolif", "tgfb"),
    "five_latent" = c("tmb", "tcell", "prolif", "tgfb", "pretreat"),
    "five_latent_hmf" = c("tmb", "tcell", "prolif", "tgfb", "pretreat_comp"),
    "five_latent_purity" = c("tmb", "tcell", "prolif", "tgfb", "pretreat", "purity")
)

apply_hmf_mods <- function (df, model, features) {
    model_features <- models[[features]]
    mod_lr <- hmf_models[[features]][[model]]["mod_lr"]$mod_lr
    mod_os <- hmf_models[[features]][[model]]["mod_os"]$mod_os
    if (nrow(df) > 0) {
        X <- as.matrix(df %>% select(all_of(model_features)))
        get_preds_X(X, df %>% pull(patient_id), mod_lr, mod_os)
    }
}
get_preds_X <- function (X, ids, mod_lr, mod_os) {
    pred_lr <- predict(mod_lr, X, predcontrib = FALSE)
    pred_os <- predict(mod_os, X, predcontrib = FALSE)
    shaps_lr <- data.frame(predict(mod_lr, X, predcontrib = TRUE))
    shaps_os <- data.frame(predict(mod_os, X, predcontrib = TRUE))
    names(shaps_lr) <- paste0("shap_lr_", names(shaps_lr))
    names(shaps_os) <- paste0("shap_os_", names(shaps_os))
    data.frame(patient_id = ids, pred_lr = pred_lr, pred_os = pred_os, shaps_lr, shaps_os)
}

### compute performance 
measure <- function( validation_go, i, j, k = "all"){
    
    if( k != "all"){
       tmp <- validation_go %>% filter(Study == i, model == j, tissue == k)    
    } else {
        tmp <- validation_go %>% filter(Study == i, model == j)
    }
    
    if( length(unique(tmp$bor)) > 1 ){
        measure_auc <- roc(tmp$bor, tmp$pred_lr, plot = FALSE, auc = TRUE)$auc
    } else {
        measure_auc <- NA
    }
    measure_c_index <- concordance.index( tmp$pred_os, tmp$os_days, tmp$os_event)$c.index
    data.frame( Study = i, model = j, tissue = k, auc = measure_auc, c_index = measure_c_index)

}
