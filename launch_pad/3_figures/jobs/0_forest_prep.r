options(warn=-1)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)
library(survival)

go <- readRDS(paste0(TMP_DIR, "validation-go.Rds")) %>% select(-contains("cluster5"))

full_mod_main <- function( i ){
    if (grepl("Lyon-162", i)){
         "tcell + tgfb + prolif"
    } else if (grepl("Lyon", i)){
         "pretreat + tcell + tgfb + prolif"
    } else if (grepl("DRUP", i) | i == "MARIATHAN") {
        "tmb + tcell + tgfb + prolif"
    } else {
        "tcell + tgfb + prolif + tmb + pretreat"
    }
}
full_mod_clust5 <- function(i){
    if (grepl("Lyon-162", i)){
         "tcell_cluster5 + tgfb_cluster5 + prolif_cluster5"
    } else if (grepl("Lyon", i)){
         "pretreat + tcell_cluster5 + tgfb_cluster5 + prolif_cluster5"
    } else if (grepl("DRUP", i) | i == "MARIATHAN") {
        "tmb + tcell_cluster5 + tgfb_cluster5 + prolif_cluster5"
    } else {
        "tcell_cluster5 + tgfb_cluster5 + prolif_cluster5 + tmb + pretreat"
    }
}
full_mod_sets <- function(i){
    if (grepl("Lyon-162", i)){
         "tcell_set + tgfb_set + prolif_set"
    } else if (grepl("Lyon", i)){
         "pretreat + tcell_set + tgfb_set + prolif_set"
    } else if (grepl("DRUP", i) | i == "MARIATHAN") {
        "tmb + tcell_set + tgfb_set + prolif_set"
    } else {
        "tcell_set + tgfb_set + prolif_set + tmb + pretreat"
    }
}
full_mod <- function(i, sets = "sets"){
    if( sets == "clusters"){
        full_mod_main(i)
    } else if (sets == "sets") {
        full_mod_sets(i)
    } else if (sets == "clusters5") {
        full_mod_clust5(i)
    }
}

model_summary <- function( df, features = "all", i = "HMF-CPCT", model = "glm", sets = "sets1", extra_covariates = NA ){

    ### 0 - filter to correct Study data 
    ready <- df %>% filter(Study_cohort == i)
   
    ### Update feature names by cluster
    if( sets == "sets"){
        features <- ifelse( features %in% c("tcell", "tgfb", "prolif"), 
                            paste0(features, "_set"), 
                            features)
    }

    ### 1 - Construct model formula
    y <- ifelse(model == "glm", "bor", "Surv(os_days, os_event)")
    x <- ifelse( features == "all", full_mod( i, sets ), features)
    
    ### add covariates 
    covariates <- ""
    if( length( unique(ready$tissue) ) > 1) { covariates <- "+ as.factor(study_tissue_cohort)"}
    if( grepl("HMF", i)){ covariates <- paste0( covariates, "+ purity" ) }
    if( !is.na(extra_covariates)) { covariates <- paste0( covariates, extra_covariates) }
    
    formula_go = as.formula( paste( y, x = paste0(x, covariates), sep = "~") )

    ### 2 - Run models
    if( model == "glm"){
        set <- do.call("glm", list( formula = formula_go, family = "binomial", data = as.name("ready")))
    } else {
        set <- do.call("coxph", list( formula = formula_go, data = as.name("ready")))
    } 

    ### 3 - Extract, tag, share output 
    go <- data.frame(summary(set)$coefficients) %>% rownames_to_column("feature")
    go$features <- features
    go$cohort <- i
    go$model <- model
    go$covariates <- covariates
    go$sets <- sets
    go
}

model_summary2 <- function( df, features = "all", i = "HMF-CPCT", model = "glm", sets = "sets1", extra_covariates = NA ){
    out <- tryCatch( model_summary( df, features, i, model, sets, extra_covariates), error = function(e) NULL)
    return(out)
}

studies <- c(
             "HMF-CPCT", 
             "HMF-CPCT-skin", 
             "HMF-CPCT-lung", 
             "HMF-CPCT-bladder", 
             "HMF-CPCT-other",
             "HMF-CPCT-low-purity",
             "HMF-DRUP", 
             "HMF-WIDE", 
             "INSPIRE", 
             "VHIO", 
             "RAVI", 
             "MARIATHAN", 
             "PARKER",
             "Lyon", 
             "External Studies")
features <- c("tcell", "tgfb", "prolif", "tmb", "pretreat", "all" )

rr <- data.frame()
for ( i in studies ){
  for ( j in features){
    for( k in c("clusters", "sets")){
    tmp <- model_summary2( df = go, features = j, i = i, sets = k, model = "glm" )
    rr <- rbind(rr, tmp)
  }
}}

ss <- data.frame()
for ( i in studies ){
    for ( j in features){
        for( k in c("clusters", "sets")){
            tmp <- model_summary2( df = go, features = j, i = i, sets = k, model = "coxph" )
            ss <- rbind(ss, tmp)
    }
}}

alpha <- .05
z_alpha <- qnorm(1-alpha/2)

rr_clean <- (
  rr
    %>% transmute( 
      feature, 
      est = Estimate, 
      ci_low = Estimate - z_alpha*Std..Error, 
      ci_high = Estimate + z_alpha*Std..Error, 
      p_val = Pr...z.., 
      z = z.value,
      features = features,
      cohort, 
      model, 
      sets,
      covariates
))

ss_clean <- (
  ss 
    %>% transmute( 
      feature, 
      est = coef, 
      ci_low = coef - z_alpha*se.coef., 
      ci_high = coef + z_alpha*se.coef., 
      p_val = Pr...z.., 
      z = z,
      features = features,
      cohort, 
      model, 
      sets, 
      covariates
))

base <- bind_rows(rr_clean, ss_clean)

renamer <- function(i){
    if( grepl("tcell", i)){
        "tcell"
    } else if (grepl("prolif", i)){
        "prolif"
    } else if (grepl("tgfb", i)){
        "tgfb"
    } else {
        i
    }
}
base$feature <- unlist(lapply(base$feature, renamer))

feature_map <- list(
    "tmb" = "TMB",
    "tcell" = "T-cell",
    "tgfb" = "TGFB",
    "prolif" = "Proliferation",
    "pretreat" = "Pretreatment",
    "purity" = "Purity"
)
name_mapper <- function( i ){
    if( i %in% names(feature_map)){
        feature_map[[i]]
    } else {
        i
    }
}
base$clean_feature <- unlist(lapply(base$feature, function(i) name_mapper(i)))
base$clean_feature <- factor(base$clean_feature, levels = unlist(feature_map))

study_map <- list(
    "HMF-CPCT" = "HMF CPI Overall",
    "HMF-CPCT-skin" = "Skin",
    "HMF-CPCT-lung" = "Lung",
    "HMF-CPCT-bladder" = "Bladder",
    "HMF-CPCT-other" = "Other",
    "HMF-CPCT-low-purity" = "Low purity",
    "External Studies" = "Validation Overall",
    "PARKER" = "PARKER ICI (Skin)",
    "RAVI" = "RAVI (Lung)",
    "MARIATHAN" = "MARIATHASAN (Bladder)",
    "INSPIRE" = "INSPIRE (mixed)",
    "Lyon" = "Lyon (Lung, HNC)",
    "VHIO" = "VHIO (mixed)",
    "HMF-DRUP" = "HMF-DRUP (mixed)"
)
base$clean_study <- unlist(lapply(base$cohort, function(i) study_map[[as.character(i)]]))
base$clean_study <- factor(base$clean_study, levels = rev(unlist(study_map)))

base$clean_study2 <- paste0(base$clean_study, " ", base$sets)
base$clean_study2 <- factor(base$clean_study2, 
                            levels = c(rbind(paste0( levels(base$clean_study), " clusters"), 
                                             paste0( levels(base$clean_study), " sets"))))

base$clean_model <- ifelse(base$model == "glm", "Best Overall Response", "Overall Survival")

z_cat <- function(i){
    if( is.na(i)){
        "non"
    } else if (abs(i) > 3){
        "strong"
    } else if (abs(i) > 2){
        "moderate"
    } else if(abs(i) > 1) {
        "weak"
    } else {
        "non"
    }
}

base$z_group <- unlist(lapply(base$z, z_cat))

better_or_worse <- function( est, model ){
    if( model == "glm"){
        oo <- ifelse(est > 0, "better", "worse")
    } else if (model == "coxph") {
        oo <- ifelse(est < 0, "better", "worse")
    }
    oo
}

base <- (
  base 
    %>% mutate(tmp_est = ifelse(model == "glm", est, -est))
    %>% mutate( better = ifelse( tmp_est > 0, "yes", "no"))
    %>% mutate( better = ifelse( z_group == "non", "non", better))
)

ready <- (
  base
    %>% filter(feature %in% c("tmb", "tcell", "tgfb", "prolif", "pretreat", "purity"))
    %>% transmute(
     est, 
     ci_low, 
     ci_high, 
     cohort, 
     clean_feature, 
     clean_study, 
     clean_study2,
     clean_model,
     z_group, 
     better, 
     covariates,
     features,
     sets
  )
)

saveRDS( ready, paste0(TMP_DIR, "forest-ready.Rds"))
