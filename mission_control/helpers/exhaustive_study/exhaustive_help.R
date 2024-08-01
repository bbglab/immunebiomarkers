library(tidyverse)

args_map <- list()
args_map[['bor']] <- c("response", "Y_best_response_binary")
args_map[['pfs']] <- c("survival", "Surv(Survival_time_to_pfs_event, Survival_pfs_event)")
args_map[['os']] <- c("survival", "Surv(Survival_time_to_os_event, Survival_os_event)")

create_CPI_data_store <- function( df, tissues){

    response <- list();
    survival <- list(); 
    for (i in c('all','all_adj','all_purity_tissue', 'all_age_biopsy_adj','all_age_biopsy_purity_adj')){
        response[[i]] <- df %>% filter( Filter_meta_responseMeasured == "Yes" ) %>% drop_na(Y_best_response_binary)
        survival[[i]] <- df
    }
    for (i in c('tmb_tcell_clin','tmb_tcell_clin_adj', 'tmb_tcell_clin_age_biopsy_adj', 'tmb_tcell_clin_age_biopsy_purity_adj')){
        response[[i]] <- df %>% drop_na(tcell, Y_best_response_binary) %>% filter( Filter_meta_responseMeasured == "Yes" ) 
        survival[[i]] <- df %>% drop_na(tcell)
    }
    tissues <- c("bladder", "lung", "other", "skin")
    for (i in tissues ) {
        response[[i]] <- response[['all']] %>% filter(tissue == i)
        response[[paste0(i, "_adj")]] <- response[['all']] %>% filter(tissue == i)
        survival[[i]] <- survival[['all']] %>% filter(tissue == i)
        survival[[paste0(i, "_adj")]] <- survival[['all']] %>% filter(tissue == i)                    
    }
    cpi_store = list();
    cpi_store[['os']] = survival
    cpi_store[['pfs']] = survival
    cpi_store[['bor']] = response
    cpi_store

}
scale_the_data <- function( df ){
    col_types <- sapply( df, typeof )
    responses <- colnames(df %>% select(Y_best_response_binary, Y_best_response_time_in_days, Y_relapse, contains("Surv")))
    step1 <- col_types[which(col_types != "character")]
    step2 <- names(step1[-which(names(step1) %in% responses)])
    df %>% mutate_at( vars(all_of(step2)), ~ (scale(.) %>% as.vector))
}
feature_update <- function( feature, dataset ){
    if ( dataset == "tmb_tcell_clin_adj" ) {
        feature <- paste0(feature,"+ pretreat + tcell + tmb + as.factor(tissue)")
    } else if (dataset == "tmb_tcell_clin_age_biopsy_adj"){
        feature <- paste0(feature,"+ pretreat + tcell + tmb + age + as.factor(tissue) + as.factor(biopsy)")
    } else if (dataset == "tmb_tcell_clin_age_biopsy_purity_adj"){
        feature <- paste0(feature,"+pretreat + tcell + tmb + age + purity +as.factor(tissue)+as.factor(biopsy)")
    } else if ( dataset == "residuals" ) {
        feature <- paste0(feature,"+ tgfb + pretreat + tcell + tmb + purity + as.factor(tissue)")
    } else if ( grepl( "recist_response_adj", dataset ) ){
        feature <- paste0(feature,"+ as.factor(Y_best_response_binary) + as.factor(tissue)")
    } else if ( dataset %in% c("all_adj", "responders_adj", "non_responders_adj")){
        feature <- paste0(feature,"+ as.factor(tissue)")
    } else if ( dataset == "all_age_biopsy_adj") {
       paste0(feature,"+ age + as.factor(tissue) + as.factor(biopsy)")
    } else if ( dataset == "all_purity_tissue" ) {
       paste0(feature,"+ purity + as.factor(tissue)")
    } else if ( dataset == "all_age_biopsy_purity_adj") {
       paste0(feature,"+ age + purity + as.factor(tissue) + as.factor(biopsy)")
    } else if ( dataset %in% c("skin_adj", "lung_adj", "bladder_adj", "other_adj")){
       feature <- paste0(feature,"+ pretreat + tcell + tmb + purity") 
    } else if ( dataset %in% c("tmb_tcell_clin")){
       feature <- paste0(feature,"+ pretreat + tcell + tmb") 
    } else {
        feature
    }
}
feature_update_extension <- function(feature, dataset) {
    ifelse (  grepl("ecist_response_adj", dataset),
              paste0(feature, "+ as.factor(Y_best_response_binary) + as.factor(tissue)"), 
              paste0(feature, "+ as.factor(tissue)")
            )
}
get_formula <- function( y, x ){
    tmp <- strsplit(x,"+",fixed = T)[[1]]
    x_reduced <- ifelse(length(tmp) == 1, "1", paste(tmp[-1], collapse='+' ))
    full    <- as.formula(paste(y, x, sep = "~"))
    reduced <- as.formula(paste(y, x_reduced, sep = "~"))
    list( "formula_full" = full, "formula_reduced" = reduced)
}
get_model_stats <- function( y, x, dataset, col_type, model_type ){
    
    formulas <- get_formula(y, x)
    formula_full <- formulas[["formula_full"]]
    
    if (model_type == "response"){
        model_full <- do.call("glm",  list(formula = formula_full, data = as.name(dataset), family="binomial"))
    } else if ( model_type == "survival"){
        model_full <- do.call("coxph",list(formula = formula_full,data = as.name(dataset)))
    }   
    if( col_type == "factor"){     
        est <- NA; se <- NA; ph <- 1
        formula_reduced <- formulas[["formula_reduced"]]
        
        if ( model_type == "response"){
            model_reduced <- do.call("glm", list(formula = formula_reduced, data = as.name(dataset), family = "binomial")) 
        } else if ( model_type == "survival") {
            model_reduced <- do.call("coxph",list(formula = formula_reduced,data = as.name(dataset)))   
            ph <- cox.zph(model_full)$table["GLOBAL", "p"]
        }
        lrt <- anova(model_reduced, model_full)
        df <- as.numeric(lrt[2,3]); 
        
        if (model_type == "response"){
            p_val <- 1-pchisq(as.numeric(lrt[2,4]), df)
        } else if (model_type == "survival") {
            p_val <- as.numeric(lrt[2,4])
        }        
    } else if (col_type == "continuous") {
        ph <- 1
        if ( model_type == "response"){
            coefs <- c(summary(model_full)$coefficients[2,c(1,2,4)])
        } else if ( model_type == "survival"){
            coefs <- c(summary(model_full)$coefficients[1,c(1,3,5)])
            ph <- cox.zph(model_full)$table["GLOBAL", "p"]
        }
        est <- coefs[1]; se <- coefs[2]; p_val <- coefs[3]; df <- 1
    }    
    list( "est" = est, "se" = se, "p_val" = p_val, "df" = df, "ph" = ph)

}
tidy_stats <- function( est, se, p_val, df, feature, dataset, col_type, model_type, ph ){
    output <- data.frame( 
                "est"      = est, 
                "se"       = se, 
                "p_val"    = p_val, 
                "df"       = df, 
                "feature"  = feature, 
                "dataset"  = dataset, 
                "col_type" = col_type, 
                "model_type" = model_type,
                "ph"       = ph
    )
    as_tibble( output )
}
get_stats <- function( response, feature, dataset, col_type, model_type, type = "normal" ){   
    
    model_x <- ifelse(type == "normal",feature_update(feature,dataset),feature_update_extension(feature,dataset)) 
    stats <- get_model_stats( response, model_x, dataset, col_type, model_type )
    
    tidy_stats( est = stats[['est']], se = stats[['se']], p_val = stats[['p_val']], df = stats[['df']],
                model_type = model_type, feature = feature, dataset = dataset, col_type = col_type, ph = stats[['ph']])   
}
get_stats2 <- function( response, feature, dataset, col_type, model_type, type = "normal" ){
    out <- tryCatch( get_stats( response, feature, dataset, col_type, model_type, type ), error = function(e) NULL)
    return(out)
}
get_cors <- function(df,top,i){
    df <- data.frame(df)
    (
        data.frame(cor(df[,c(top,i)],use = "pairwise.complete.obs"))
            %>% tibble::rownames_to_column(var = "feature")
            %>% filter(feature == i)
            %>% select(all_of(c("feature", top)))
    )
}
get_cors2 <- function(df,top,i){
    out <- tryCatch(get_cors(df,top,i), error = function(e) NULL)
    return(out)
}
give_me_pcs <- function( df, match, scale = FALSE ){
    hold <- (
        df
            %>% select("sampleId",contains(match))
            %>% select_if(~sum(!is.na(.)) > 0)
            %>% drop_na()
            %>% mutate_at(vars(-sampleId), ~scale(.x))
            %>% select_if(~sum(!is.na(.)) > 0)
    )
    twirl <- prcomp(as.matrix(hold[,-1]))
    keep <- min( min(which(cumsum(twirl$sdev^2)/sum(twirl$sdev^2) > .95)),20)
    meld <- data.frame(cbind("sampleId" = hold$sampleId, twirl$x[,1:keep]))
    idx <- seq(ncol(meld))[-1]
    meld[idx] <- lapply(meld[idx], function(x) as.numeric(as.character(x)))
    meld %>% rename_at(vars(-sampleId), ~ paste0(paste0(substr(match,1,nchar(match)-1),".pc_"),.x))
}
