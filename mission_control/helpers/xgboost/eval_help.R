### 0 - Model Fitting 
eval_map <- list()
eval_map[['lr']] <- c("binary:logistic", "logloss")
eval_map[['pfs']] <- c("survival:cox", "cox-nloglik")
eval_map[['os']] <- c("survival:cox", "cox-nloglik")

xgb_cv <- function( Y, X, params, model, base_model = NULL){ 
    xgb_cv <- xgb.cv(params = params, 
                 data = X, 
                 label = Y, 
                 nrounds = 2000, 
                 nfold = 5,
                 objective = eval_map[[model]][1],  
                 eval_metric = eval_map[[model]][2], 
                 verbose = 0, 
                 early_stopping_rounds = 10,
                 xgb_model = base_model)   
    xgb_cv
}

grid_fit <- function( Y, X, hyper_grid, model, base_model = NULL){
    
    for(i in 1:nrow(hyper_grid)) {
        params <- list(
            eta = hyper_grid$eta[i],
            max_depth = hyper_grid$max_depth[i],
            min_child_weight = hyper_grid$min_child_weight[i],
            subsample = hyper_grid$subsample[i],
            colsample_bytree = hyper_grid$colsample_bytree[i])

        tuner <- xgb_cv( Y, X, params, model, base_model)   
    
        if( model == "lr"){
            errors <- tuner$evaluation_log$test_logloss_mean
        } else if (model %in% c("pfs","os")){
            errors <- tuner$evaluation_log$test_cox_nloglik_mean
        }
        hyper_grid$best_iteration[i] <- tuner$best_iteration
        hyper_grid$error[i] <- min(errors)    
    }
    
    best_fit <- hyper_grid %>% dplyr::arrange(error) %>% head(1)
    params <- list(
        eta = best_fit$eta,
        max_depth = best_fit$max_depth,
        min_child_weight = best_fit$min_child_weight,
        subsample = best_fit$subsample,
        colsample_bytree = best_fit$colsample_bytree
    )
    best <- xgboost(
                params = params, 
                data = X, 
                label = Y, 
                nrounds = best_fit$best_iteration, 
                objective = eval_map[[model]][1], 
                eval_metric = eval_map[[model]][2], 
                verbose = 0,
                xgb_model = base_model)
    list( 'best_model' = best, 'stats' = best_fit )
}    

get_Y <- function( i, j, stores, model) {  ## i represents 'test' or 'train', j represents 'tissue'    
    if( model == "lr"){
        as.matrix(stores[[i]][[j]] %>% select( response ))
    } else if (model == "pfs"){
        as.matrix(stores[[i]][[j]] %>% select( pfs ))
    } else if (model == "os"){
        as.matrix(stores[[i]][[j]] %>% select( os ))
    }
}
get_X <- function( i, j, stores, model_features) {
    as.matrix(stores[[i]][[j]] %>% select( all_of(model_features) ))
}
get_time_to_event <- function( i, j, stores) stores[[i]][[j]] %>% pull( time_to_event )
get_event_status <- function( i, j, stores) stores[[i]][[j]] %>% pull( event_status )

evaluate <- function( eval_ready, model, parameter_grid, model_features, complete, base_tissue = "all", folds = 5){ 
    ## eval_ready from eval-0
    
    for( i in names(eval_ready[[model]])) assign(i, eval_ready[[model]][[i]])  ### extracts df and gps
    
    tissues <- c("all","skin","lung","bladder", "other")
    stores <- list() 
    models <- list() 
    stats <- list()
    preds <- list()
    evals <- list()
    
    #print(folds)
    folds <- createFolds(gps, k = folds, list = FALSE)
    idx <- which(folds == 1)
    df_train <- df[-idx,]; df_test <- df[idx,]
    
    ### store data structures for later eval ### 
    stores[['train']][['all']] <- df_train
    if(complete == TRUE) {
        df_test <- df_test %>% drop_na()
        stores[['test']][['all']] <- df_test
    } else {
        stores[['test']][['all']] <- df_test
    }
    for( i in tissues[-1]) {
        stores[['train']][[i]] <- df_train %>% filter( tissue == i )    
        stores[['test']][[i]] <-  df_test %>% filter( tissue == i ) 
    }
   
    ### fit base model, using base_tissue input
    base_model <- grid_fit(get_Y("train",base_tissue,stores,model), 
                           get_X("train",base_tissue,stores,model_features), 
                           parameter_grid, 
                           model )$best_model
 
    for( i in tissues) { ## models
        simple_model <- grid_fit(get_Y("train",i,stores,model), 
                                 get_X("train",i,stores,model_features), 
                                 parameter_grid, 
                                 model )
        warm_model   <- grid_fit(get_Y("train",i,stores,model), 
                                 get_X("train",i,stores,model_features), 
                                 parameter_grid, 
                                 model, 
                                 base_model = base_model)
    
        models[['simple']][[i]] <- simple_model$best_model
        models[['warm']][[i]] <- warm_model$best_model
        
        stats[['simple']][[i]] <- simple_model$stats
        stats[['warm']][[i]] <- warm_model$stats     
    }
        
    for (i in tissues[-1]){ ## preds
        preds[['simple']][[i]] <- predict(     models[['simple']][[i]],  get_X('test', i, stores, model_features))
        preds[['warm']][[i]] <- predict(       models[['warm']][[i]],    get_X('test', i, stores, model_features))
        preds[[base_tissue]][[i]] <- predict(  base_model,               get_X('test', i, stores, model_features))
    }
    preds[[base_tissue]][['all']] <- predict(  base_model,               get_X('test', 'all', stores, model_features))
    preds[['simple']][['all']] <- c(preds[['simple']][['skin']], 
                                    preds[['simple']][['lung']], 
                                    preds[['simple']][['bladder']], 
                                    preds[['simple']][['other']])
    preds[['warm']][['all']] <- c(preds[['warm']][['skin']], 
                                  preds[['warm']][['lung']], 
                                  preds[['warm']][['bladder']], 
                                  preds[['warm']][['other']])

    if (model == "lr"){ ## evals
        for (i in tissues[-1]){
            Y_test <- as.numeric(get_Y("test", i, stores, model)); #print(Y_test); print(i)
            
            if(sum(Y_test) != 0) {
                evals[['simple']][[i]]      <- roc( Y_test, preds[['simple']][[i]] )$auc
                evals[['warm']][[i]]        <- roc( Y_test, preds[['warm']][[i]]  )$auc
                evals[[base_tissue]][[i]]   <- roc( Y_test, preds[[base_tissue]][[i]]  )$auc
            } else {
                evals[['simple']][[i]]      <- NA
                evals[['warm']][[i]]        <- NA
                evals[[base_tissue]][[i]]   <- NA
            }
        }
        evals[[base_tissue]][['all']] <- roc(as.numeric(get_Y("test","all", stores, model)),preds[[base_tissue]][['all']])$auc
        test_all <- as.numeric(c(get_Y("test","skin", stores, model), 
                                 get_Y("test","lung", stores, model),
                                 get_Y("test","bladder",stores, model),
                                 get_Y("test","other",stores, model)))
        evals[['simple']][['all']] <- roc( test_all, preds[['simple']][['all']] )$auc
        evals[['warm']][['all']] <- roc( test_all, preds[['warm']][['all']] )$auc 
    } else {
        for (i in tissues[-1]){
            evals[['simple']][[i]] <- concordance.index( preds[['simple']][[i]], 
                                                         get_time_to_event("test",i, stores), 
                                                         get_event_status("test",i, stores))$c.index
            evals[['warm']][[i]] <- concordance.index(   preds[['warm']][[i]],   
                                                         get_time_to_event("test",i, stores), 
                                                         get_event_status("test",i, stores))$c.index
            evals[[base_tissue]][[i]] <- concordance.index( preds[[base_tissue]][[i]],   
                                                         get_time_to_event("test",i, stores), 
                                                         get_event_status("test",i, stores))$c.index
     }
       evals[[base_tissue]][['all']] <- concordance.index( preds[[base_tissue]][['all']],   
                                                       get_time_to_event("test",'all', stores), 
                                                       get_event_status("test",'all', stores))$c.index 
       tte_all <- as.numeric(c(get_time_to_event("test","skin",stores),
                               get_time_to_event("test","lung",stores),
                               get_time_to_event("test","bladder",stores),
                               get_time_to_event("test","other",stores)))
       es_all <- as.numeric(c(get_event_status("test","skin",stores),
                              get_event_status("test","lung",stores),
                              get_event_status("test","bladder",stores),
                              get_event_status("test","other",stores)))
       evals[['simple']][['all']] <- concordance.index( preds[['simple']][['all']], tte_all, es_all )$c.index
       evals[['warm']][['all']] <-   concordance.index( preds[['warm']][['all']], tte_all, es_all )$c.index
   }
    list( "evals" = evals, "stats" = stats, "preds" = preds)       
}

evaluate2 <- function( eval_ready, model, parameter_grid, model_features, complete, base_tissue = "all", folds = 5){
	out <- tryCatch(evaluate(eval_ready, model, parameter_grid, model_features, complete, base_tissue, folds), error = function(e) NULL)
	return(out)
}
