wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_prep.R"))
library(tidyverse)
library(stringr)

boom <- readRDS(paste0(TMP_DIR,"exhaustive-combine.Rds"))

boom <- boom %>% filter( log10_p > .01, col_type != "factor", feature != "pretreat_comp")
boom$est <- ifelse( boom$col_type == "factor", 0, boom$est)
boom$cor_pretreat <- ifelse( boom$feature == "clinical_pre_to_post_treatment_time", 0,  boom$cor_pretreat)

boom$feature_group <- unlist(lapply(boom$feature, function(i) strsplit(i, "_")[[1]][1]))
boom$feature_group <- ifelse( boom$feature %in% c("tcell", "tgfb", "prolif", "pdl1"), "isofox", boom$feature_group)
boom$feature_group <- ifelse( boom$feature %in% c("sv", "tmb","purity"), "somatic", boom$feature_group)   
boom$feature_group <- ifelse( boom$feature %in% c("pretreat", "age"), "clinical", boom$feature_group)                                       

boom$p_adj_by <- p.adjust(boom$p_val, method = "BY")
boom$p_adj_bh <- p.adjust(boom$p_val, method = "BH")
boom$by_05_fdr <- boom %>% filter( p_adj_by < .05 ) %>% arrange( log10_p ) %>% head(1) %>% pull(p_val)
boom$bh_05_fdr <- boom %>% filter( p_adj_bh < .05 ) %>% arrange( log10_p ) %>% head(1) %>% pull(p_val)
boom$bf_05_fwe <- .05/nrow(boom)

boom$feature_group <- apply( boom %>% select(feature, feature_group), 1, sort_the_drivers)
boom$feature_group <- ifelse(boom$feature == "hla_lilac_mut_hla", "somatic", boom$feature_group)
boom$plot_est <- ifelse(boom$model %in% c("bor", "relapse", "surv_at_t"), exp(boom$est), 1/exp(boom$est))
boom$big_group <- unlist(lapply( boom$feature_group, big_grouper))
boom$little_group <- unlist(lapply( boom$feature_group, little_grouper))
boom$Direction <- ifelse(boom$plot_est > 1, "Better", "Worse")
boom$Group <- apply( boom %>% select(feature, feature_group), 1, grouper)
boom <- boom %>% filter( plot_est < 5, plot_est > .2)

boom <- boom %>% mutate(big_group = ifelse( big_group %in% c("CNV", "SVs"), "CNV/SVs",big_group)) %>% filter(feature != "clinical_systemic_composite")

grouper <- function( big_group, cor_tcell, cor_prolif, cor_tgfb ){
    if( big_group != "RNA"){ big_group } 
    else {
      if( cor_tcell > .5){  "RNA: T-cell"} 
      else if (cor_tgfb > .5){  "RNA: TGFB"} 
      else if (cor_prolif > .5){ "RNA: Proliferation" } 
      else { "RNA: Remaining" }
    }
}

boom <-
boom %>% 
  rowwise() %>% 
  mutate(discovery_group = grouper( big_group, cor_tcell, cor_prolif, cor_tgfb )) %>% 
  mutate(discovery_group = factor(discovery_group, 
         levels = c("Somatic", "RNA: T-cell", "RNA: TGFB", "RNA: Proliferation", "RNA: Remaining",
                    "CNV/SVs", "Clinical", "HLA" ))) %>% 
  ungroup() %>% 
  mutate(Type = discovery_group)

saveRDS( boom, paste0(TMP_DIR,'exhaustive-plots-base.Rds'))
