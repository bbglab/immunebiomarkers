wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
library(tidyverse)

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds"))
features <- readRDS(paste0(TMP_DIR, "exhaustive-features-go.Rds"))

top <- c("pretreat", "tmb", "tcell","prolif", "tgfb", "purity")

cors <- list()
system.time(
for (feature in features){
    cors[[feature]] <- get_cors2(cpi, top, feature)
})

saveRDS(cors, paste0(TMP_DIR,"exhaustive-cors-list.Rds"))

cors_go <- 
dplyr::bind_rows(cors) %>% 
  rename_with(~ paste0("cor_", .x)) %>% 
  rename(feature = cor_feature)

saveRDS(cors_go, paste0(TMP_DIR,"exhaustive-cors.Rds"))
