wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
library(tidyverse)

results <- data.frame()
for( i in list.files(paste0(TMP_DIR, 'exhaustive_study/'))){
    results <- rbind(results, readRDS( paste0(TMP_DIR, 'exhaustive_study/',i)))
}

saveRDS( results, paste0( TMP_DIR, paste0( "exhaustive-stats.Rds" ) ))
