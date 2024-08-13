wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/extensions_help.R"))

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds")); 
non_cpi <- readRDS(paste0(TMP_DIR, "validation-hmf-noncpi-go.Rds")); 

features <- c( "Y_best_response_binary", "Survival_os_event", "Survival_time_to_os_event", 
                "Filter_meta_responseMeasured", 
                "tissue", "tmb", "tcell", "prolif", "tgfb", "pretreat")

cpi <- cpi %>% select(all_of(features)) %>% mutate(Treatment = "Immunotherapy")
non_cpi <- non_cpi %>% select(all_of(features)) %>% mutate(Treatment = "Non-Immunotherapy")

scale_the_data2 <- function (df) {
    col_types <- sapply(df, typeof)
    responses <- colnames(df %>% select(Y_best_response_binary, contains("Surv")))
    step1 <- col_types[which(col_types != "character")]
    step2 <- names(step1[-which(names(step1) %in% responses)])
    df <- df %>% mutate_at(vars(all_of(step2)), ~(scale(.) %>% as.vector))
    df
}

cpi_scaled <- create_CPI_data_store2(         scale_the_data2(cpi),     cpi %>% pull(Treatment ))
non_cpi_scaled <- create_CPI_data_store2( scale_the_data2(non_cpi), non_cpi %>% pull(Treatment ))

features <- c("tmb", "prolif", "tgfb", "tcell", "pretreat")

print("CPI")
stats_cpi <- data.frame()
for (model in names(cpi_scaled)){
    
    boom <- cpi_scaled[[model]];        
    for(i in names(boom)) assign(i,boom[[i]])
    
    model_type <- args_map[[model]][1]; 
    response <- args_map[[model]][2]
   
    for (feature in features){
            stats_i <- get_stats2( response, feature, "all_adj", "continuous", model_type, type = "extension" )
            stats_cpi <- rbind( stats_cpi, stats_i )
    }
}
stats_cpi$group = "CPI"

print("Non-CPI")
stats_non_cpi <- data.frame()
for (model in names(non_cpi_scaled)){
    
    boom <- non_cpi_scaled[[model]]
    for(i in names(boom)) assign(i,boom[[i]])
    
    model_type <- args_map[[model]][1]
    response <- args_map[[model]][2]
   
    for (feature in features){
            print(feature)
            flush.console()
            stats_i <- get_stats2( response, feature, "all_adj", "continuous", model_type )
            stats_non_cpi <- rbind( stats_non_cpi, stats_i )
    }
}
stats_non_cpi$group = "Non-CPI"

stats <- rbind(stats_cpi, stats_non_cpi)
stats$feature <- as.character(stats$feature)
stats$dataset <- as.character(stats$dataset)
stats <- stats %>% mutate( log10_p = -log10(p_val) )
stats$plot_est <- ifelse(stats$model_type == "response", exp(stats$est), 1/exp(stats$est))
stats$nice_response <- ifelse(stats$model_type == "response", "BOR", "Overall Survival")

stats$clean_feature <- sapply(stats$feature, function(i) ifelse( i %in% names(name_map), name_map[[i]], i))

saveRDS( stats, paste0(TMP_DIR, "supplement-five-extension.Rds"))
