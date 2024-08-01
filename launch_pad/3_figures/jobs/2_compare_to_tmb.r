library(dplyr)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))

five <- readRDS(paste0( TMP_DIR, "five_factor_responders.rds"))
tmb <- readRDS(paste0( TMP_DIR, "tmb_responders.rds"))

write.csv( rbind( tmb %>% rename(lr_gp = tmb_high) , five ), file = paste0(TMP_DIR, "compare.csv"))

paste0(TMP_DIR, "compare.csv")
