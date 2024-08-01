library(tidyverse)

args_map <- list()
args_map[['bor']] <- c("response", "Y_best_response_binary")
args_map[['pfs']] <- c("survival", "Surv(Survival_time_to_pfs_event, Survival_pfs_event)")
args_map[['os']] <- c("survival", "Surv(Survival_time_to_os_event, Survival_os_event)")

### define covariate vectors here 
add_covariates <- function(feature, covariates){ 
    if( covariates == "tissue" ){
        feature <- paste0(feature,"+ as.factor(clinical_tumor_location_group)")
    } else if ( covariates == "purity_tissue" ){
        paste0(feature,"+ purity + as.factor(clinical_tumor_location_group)")
    } else if ( covariates == "age_biopsy_purity_tissue") {
        paste0(feature,"+ age + as.factor(biopsy) + purity + as.factor(clinical_tumor_location_group)")
    } else if ( covariates == "purity") {
        paste0(feature,"+ purity")
    } else if ( covariates == "residuals" ) {
        paste0(feature,"+ tmb + tcell + pretreat + purity + as.factor(clinical_tumor_location_group)")
    } else if ( covariates == "residuals2" ) {
        paste0(feature,"+ tmb + tcell + pretreat + tgfb + purity + as.factor(clinical_tumor_location_group)")
    } else {
	feature
    }
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
        
        lrt <- anova(model_reduced, model_full)
        df <- as.numeric(lrt[2,3]); 
        
        if (model_type == "response"){
            p_val <- 1-pchisq(as.numeric(lrt[2,4]), df)
        } else if (model_type == "survival") {
            p_val <- as.numeric(lrt[2,4])
        }        
    }} else if (col_type == "continuous") {
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
tidy_stats <- function( est, se, p_val, df, feature, dataset, covariates, col_type, model_type, ph ){
    output <- data.frame( 
                "est"      = est, 
                "se"       = se, 
                "p_val"    = p_val, 
                "df"       = df, 
                "feature"  = feature, 
                "dataset"  = dataset, 
                "covariates" = covariates, 
                "col_type" = col_type, 
                "model_type" = model_type,
                "ph"       = ph
    )
    as_tibble( output )
}
get_stats <- function( response, feature, dataset, col_type, model_type, covariates ){
    model_x <- add_covariates(feature, covariates)
    stats <- get_model_stats( response, model_x, dataset, col_type, model_type )
    tidy_stats( 
            est = stats[['est']], 
            se = stats[['se']], 
            p_val = stats[['p_val']], 
            df = stats[['df']],
            model_type = model_type, 
            feature = feature, 
            dataset = dataset, 
            covariates = covariates, 
            col_type = col_type, 
            ph = stats[['ph']]
    )
}
get_stats2 <- function( response, feature, dataset, col_type, model_type, covariates ){
    out <- tryCatch( get_stats( response, feature, dataset, col_type, model_type, covariates ), error = function(e) NULL)
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
