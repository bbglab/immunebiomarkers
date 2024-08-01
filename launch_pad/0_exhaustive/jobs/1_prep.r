wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/exhaustive_study/exhaustive_help.R"))
library(tidyverse)

cpi <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

table(cpi %>% select(contains("clinical")) %>% pull(clinical_cpi_mechanism))

scale_the_data <- function( df ){
    col_types <- sapply( df, typeof )
    responses <- colnames(df %>% select(Y_best_response_binary, Y_best_response_time_in_days, Y_relapse, contains("Surv")))
    step1 <- col_types[which(col_types != "character")]
    step2 <- names(step1[-which(names(step1) %in% responses)])
    df %>% mutate_at( vars(all_of(step2)), ~ (scale(.) %>% as.vector))
}
cpi_scaled <- scale_the_data(cpi)

studies <- list()

### overall
studies[['all']] <- cpi_scaled

### tissues
studies[['skin']]    <- cpi_scaled %>% filter(clinical_tumor_location_group == "skin")
studies[['lung']]    <- cpi_scaled %>% filter(clinical_tumor_location_group == "lung")
studies[['bladder']] <- cpi_scaled %>% filter(clinical_tumor_location_group == "bladder")
studies[['other']]   <- cpi_scaled %>% filter(clinical_tumor_location_group == "other")

get_ready <- function ( studies ) {

    response <- list(); survival <- list()

    for( i in names(studies)){
        response[[i]] <- studies[[i]] %>% filter(Filter_meta_responseMeasured == "Yes") %>% drop_na(Y_best_response_binary)
        survival[[i]] <- studies[[i]]   
    }
    ready = list()
    ready[["survival"]] = survival
    ready[["response"]] = response
    ready
}

go <- get_ready( studies )

saveRDS(go, paste0(TMP_DIR, "exhaustive-ready.Rds"))
