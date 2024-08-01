wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/validation_help.R"))

library(tidyverse)
library(survival)
library(pROC)
library(survcomp)

validation_measure <- (
    readRDS(paste0(TMP_DIR, "validation-measure.Rds"))
        %>% select(
          patient_id, 
          build, 
          Study_cohort, 
          model, 
          tissue, 
          complete, 
          bor, 
          os_event, 
          os_days,
          os, 
          pred_lr, 
          pred_os        
        )
        %>% filter(!is.na(bor))
)

measure <- function (df, i, j, k = "all", l = TRUE, q = "1") {
    if (k != "all") {
        tmp <- df %>% filter(Study_cohort == i, model == j, tissue == k, complete == l)
    } else {
        tmp <- df %>% filter(Study_cohort == i, model == j, complete == l)
    }
    
    if (length(unique(tmp$bor)) > 1) {
        #tmp <- tmp %>% filter( !is.na(bor))
        measure_auc <- roc(tmp$bor, tmp$pred_lr, plot = FALSE, auc = TRUE)$auc
    } else {
        measure_auc <- NA
    }
    
    tmp <- tmp %>% drop_na(os_event, os_days)
    
    measure_c_index <- concordance.index(tmp$pred_os, tmp$os_days, tmp$os_event)$c.index
    response_low <- mean(tmp %>% filter(pred_lr < .1, !is.na(bor)) %>% pull(bor))
    n_response_low <- length(tmp %>% filter(pred_lr < .1, !is.na(bor)) %>% pull(bor))
    
    data.frame(
            Study = i, 
            model = j, 
            tissue = k, 
            n = nrow(tmp), 
            complete = l, 
            build = q, 
            auc = measure_auc, 
            c_index = measure_c_index, 
            response_low = response_low, 
            n_response_low = n_response_low
    )
}

lets_see <- data.frame()
for( i in unique(validation_measure $Study_cohort)){
    print(i)
    for ( j in unique(validation_measure$model)){
        for( q in unique(validation_measure$build)){
            for ( l in c(TRUE, FALSE)){
                lets_see <- rbind(lets_see, measure( validation_measure, i, j, k="all", l, q))
            }
        }
    }
}

validation_summary <- (
    lets_see 
        %>% group_by( Study, model, tissue, complete) 
        %>% summarise(auc = mean(auc, na.rm = TRUE), 
                      c_index = mean(c_index, na.rm = TRUE), 
                      non_response = mean(response_low, na.rm = TRUE))
        %>% ungroup()
        %>% rename( features = model)
        %>% gather( model, mn, -Study, -features, -tissue, -complete )
)

name_map <- list(
    "tmb_bin" = "TMB Binary",
    "base_bin" = "TMB Binary + PDL1",
    "tmb" = "TMB", 
    "base" = "TMB + PDL1",
    "rna" = "RNA only",
    "no_tmb" = "Four factors (No TMB)",
    "five_latent" = "Five Factors",
    "fiver_latent_purity" = "Five Factor (purity adjusted)"
)
validation_summary$clean_name <- factor(unlist(lapply(validation_summary$features, function(i) name_map[[i]])), levels = unname(unlist(name_map))) 

validation_summary$model <- ifelse( validation_summary$model == "auc", "Best Overall Response (AUC)", validation_summary$model)
validation_summary$model <- ifelse( validation_summary$model == "c_index", "Overall Survival (C-Index)", validation_summary$model)
validation_summary$model <- ifelse( validation_summary$model == "non_response", "Predicted Non-responder", validation_summary$model)
validation_summary$five <- ifelse( validation_summary$features %in% c("five_latent", "no_tmb"), "yes", "no")

validation_summary$model <- factor(validation_summary$model, levels = c("Best Overall Response (AUC)","Overall Survival (C-Index)","Predicted Non-responder"))

saveRDS(validation_summary, paste0(I_DIR, "validation-summary.Rds"))
