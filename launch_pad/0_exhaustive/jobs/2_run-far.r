wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help2.R"))
library(survival)

args <- commandArgs(trailing = TRUE)   ##
#args <- list("os", "all", "response_full")

model <- args[[1]]
model_type <- args_map[[model]][1]
response <- args_map[[model]][2]

dataset <- args[[2]]
covariates <- args[[3]]

boom <- readRDS(paste0(TMP_DIR, "exhaustive-ready.Rds"))[[model_type]]
assign(dataset, boom[[dataset]])

features <- readRDS(paste0(TMP_DIR, "exhaustive-features-go.Rds"))
feature_types <- lapply( boom[[dataset]][,features], typeof )

stats <- data.frame()
for (feature in features){
    col_type <- ifelse( feature_types[feature][[1]] != "character", "continuous", "factor" )
    stats_i <- get_stats2( response, feature, dataset, col_type, model_type, covariates )
    stats <- rbind( stats, stats_i )
}

stats$group <- "cpi"
stats$model <- model

O_DIR <- paste0(TMP_DIR,'exhaustive_study')
ifelse(!dir.exists(O_DIR), dir.create(O_DIR), FALSE)

saveRDS(stats, paste0(O_DIR, '/exhaustive-',model,"-", dataset,"-", covariates, '.Rds'))
