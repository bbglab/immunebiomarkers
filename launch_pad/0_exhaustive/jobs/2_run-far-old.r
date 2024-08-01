wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
library(survival)

#print(sessionInfo())

#args <- commandArgs(trailing = TRUE) ### 
args <- list("os", "all")

model <- args[[1]]
model_type <- args_map[[model]][1]
response <- args_map[[model]][2]

dataset <- args[[2]]

boom <- readRDS(paste0(TMP_DIR, "exhaustive-ready.Rds"))[[model]]
assign(dataset, boom[[dataset]])

features <- readRDS(paste0(TMP_DIR, "exhaustive-features-go.Rds"))
feature_types <- lapply( boom[[dataset]][,features], typeof )

names(boom)

stats <- data.frame()
for (feature in features[1:10]){
    col_type <- ifelse( feature_types[feature][[1]] != "character", "continuous", "factor")
    stats_i <- get_stats2( response, feature, dataset, col_type, model_type )
    stats <- rbind( stats, stats_i )
}

stats

stats$group <- "cpi"
stats$model <- model

#stats

saveRDS(stats, paste0(TMP_DIR,'exhaustive_study/exhaustive-',model,"-", dataset,'.Rds'))
