wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/validation_help.R"))

library(tidyverse)
library(xgboost)

validation_ready <- readRDS(paste0(TMP_DIR, "validation-hmf-go.Rds"))

set.seed(62220)

hmf <- validation_ready %>% drop_na(os, bor)

#run_loo_cv

hmf_loo_cv <- data.frame()
Y_base <- hmf %>% select(all_of(c("patient_id", "tissue", "bor", "os")))

for( i in c("five_latent_purity")){
    X_base <- hmf %>% select(all_of(c("patient_id", "tissue", models[[i]])))
    hmf_loo_cv <- run_loo_cv(Y_base, X_base)
    hmf_loo_cv$model <- i
}

hmf_loo_cv$pred_os2 <- ifelse( hmf_loo_cv$pred_os > 3, 3.1, hmf_loo_cv$pred_os)
hmf_loo_cv$lr_gp <-  cut( hmf_loo_cv$pred_lr,  breaks = c(0,.1,.5,    1) , labels = c("Low", "Medium", "High"))
hmf_loo_cv$os_gp <-  cut( hmf_loo_cv$pred_os2,  breaks = c(0,.5,1.5,   20) , labels = c("Low", "Medium", "High"))

saveRDS( hmf_loo_cv, paste0(TMP_DIR, "validation-loo-cv.Rds") )
