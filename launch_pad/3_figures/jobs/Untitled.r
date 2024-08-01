wd <- dirname(dirname(getwd()))
setwd(wd)

suppressMessages(library(tidyverse))

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure2_theme.R"))

go <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds")) %>% filter(feature != "pretreat")

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

som_sig <- (
    go %>% filter(
           !grepl("ciber", feature),
           feature != "clinical_systemic_composite",
           dataset == "all", 
           model == "bor",
           covariates == "age_biopsy_purity_tissue",
           p_val < by_05_fdr,
           grepl("somatic", feature_group)
    )
    %>% pull(feature)
)

cpi$clust_som <- apply(cpi %>% select(all_of(som_sig)), 1, mean)
