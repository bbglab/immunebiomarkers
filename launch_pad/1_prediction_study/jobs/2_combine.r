wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))

library(tidyverse)
I_DIR <- paste0(TMP_DIR, 'pred_study/')

FIG_FINAL_DIR

results <- data.frame()
for( i in list.files(I_DIR)){
    results <- rbind( results, readRDS(paste0(I_DIR, i)))
}

saveRDS( results, paste0( TMP_DIR, paste0( "xg-eval-results.Rds" ) ))
