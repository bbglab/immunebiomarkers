options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))

library(tidyverse)

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))

s1 <- ingredients %>% select(dataset, covariates, model, model_type, Type, feature_group, feature, est, p_val, contains("cor"))

write.csv( s1, paste0( FIG_FINAL_DIR , "supplement_exhaustive_results.csv"))
