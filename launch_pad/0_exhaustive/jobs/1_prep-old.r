wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
library(tidyverse)

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

tissues <- c('skin','lung','bladder','other')
cpi_scaled <- create_CPI_data_store( scale_the_data(cpi), tissues )

names(cpi_scaled[['os']])

saveRDS(cpi_scaled, paste0(TMP_DIR, "exhaustive-ready.Rds"))
