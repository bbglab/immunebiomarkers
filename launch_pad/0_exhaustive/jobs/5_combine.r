wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/combine_help.R"))
library(tidyverse)

stats <- readRDS(paste0(TMP_DIR,"exhaustive-stats.Rds"))
cors <- readRDS( paste0(TMP_DIR,"exhaustive-cors.Rds") ) 

house_party <- (
    stats 
      %>% left_join(cors, by = "feature")
      %>% drop_na(p_val)
      %>% mutate(log10_p = -log10(p_val))
      %>% relocate(group, dataset, model, model_type, col_type, feature, cor_pretreat, 
                   cor_tmb, cor_tcell, cor_prolif, cor_tgfb, cor_purity, est, se, p_val, log10_p)
)

saveRDS(house_party, paste0(TMP_DIR,"exhaustive-combine.Rds"))
