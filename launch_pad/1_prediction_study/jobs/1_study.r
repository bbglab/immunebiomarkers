wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/eval_help.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/paper_settings.R"))

library(xgboost)
library(tidyverse)
library(survival)
library(pROC)
library(caret)
library(survcomp)

eval_ready <- readRDS(paste0(TMP_DIR, "xg-eval-prep.Rds"))

args <- commandArgs(trailing = TRUE) ##
#args <- list("pfs", "latent_vhio", "true")
model <- args[[1]]
features <- args[[2]]
purity <- args[[3]]
model_features <- unname(unlist(model_features[features]))
if( purity == "true") {
    model_features <- c(model_features, "purity")
}
parameter_grid <- get_parameter_grid(features)

base_tissue <- "all"
N <- 1000
K <- 5

set.seed(622)

evals <- data.frame()

for (j in seq(N)) {
    print(j); flush.console()
    oo <- evaluate2(  eval_ready, 
                      model = model, 
                      parameter_grid = parameter_grid, 
                      model_features = model_features, 
                      complete = TRUE, 
                      base_tissue = base_tissue, 
                      folds = K ); 

    if (!is.null(oo)){
        eval_i <- data.frame(oo$evals) %>% tibble::rownames_to_column(var = "tissue"); 
        eval_i$model <- model; 
        eval_i$rep <- j; 
        eval_i$complete <- TRUE; 
        eval_i$features <- features
        evals <- rbind(evals, eval_i)
    }
}
evals$base_tissue <- base_tissue
evals$purity <- purity

O_DIR <- paste0(TMP_DIR,'pred_study/')
ifelse(!dir.exists(O_DIR), dir.create(O_DIR), FALSE)

saveRDS(evals,paste0( O_DIR, "xg-eval-results-",model,"-",features, "-", purity,".Rds" ))
