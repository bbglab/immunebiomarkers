wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)

ready <- readRDS(paste0(TMP_DIR,'exhaustive-plots-base.Rds'))

k <- .3

(ready %>% filter(model == "bor", covariates == "age_biopsy_purity_tissue", p_val < by_05_fdr, 
                 abs(cor_tcell) < k, abs(cor_tmb) < k, abs(cor_tgfb) < k, abs(cor_prolif) < k, abs(cor_pretreat) < k)
      %>% select(feature, est, p_val, contains("cor"), by_05_fdr))

go <- ready %>% select(dataset, Type, feature, covariates, model, p_val, contains("_05_"))

see <- go %>% mutate(model2 = ifelse( model == "os" & covariates == "residuals", "os_resid", model))

#head(go)

(
 go 
   %>% mutate(model = ifelse( model == "os" & covariates == "residuals", "os_resid", model))
   %>% mutate(covariates = ifelse( model == "os_resid" & covariates == "residuals", "age_biopsy_purity_tissue", as.character(covariates)))
   %>% filter(dataset == "all", covariates == "age_biopsy_purity_tissue")
   %>% group_by(Type, model)
   %>% mutate(ct = n(), 
                 bh = p_val < bh_05_fdr,
                 by = p_val < by_05_fdr,
                 bf = p_val < bf_05_fwe)
   %>% group_by(Type, feature)
   %>% summarise( bh = sum(bh) > 0, by = sum(by) > 0, bf = sum(bf) > 0)
   %>% group_by(Type)
   %>% summarise(ct = n(), 
                 bh = sum(bh), 
                 p_bh = round(100*sum(bh)/ct,1),
                 by = sum(by), 
                 p_by = round(100*sum(by)/ct,1),
                 bf = sum(bf),
                 p_bf = round(100*sum(bf)/ct,1)) 
)

(step1 
    %>% group_by(Type, feature)
    %>% summarise( bh = sum(bh) > 0, by = sum(by) > 0, bf = sum(bf) > 0)
    %>% group_by(Type)
    %>% summarise(ct = n(), 
                  bh = sum(bh), 
                  p_bh = round(100*sum(bh)/ct,1),
                  by = sum(by), 
                  p_by = round(100*sum(by)/ct,1),
                  bf = sum(bf),
                  p_bf = round(100*sum(bf)/ct,1))
)
