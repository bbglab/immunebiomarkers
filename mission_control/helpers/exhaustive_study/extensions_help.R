create_CPI_data_store2 <- function(df, treatments){
    survival <- list()
    survival[["all"]] <- df
    survival[["non_responders_adj"]] <- df %>% filter(Y_best_response_binary == 0)
    survival[["responders_adj"]] <- df %>% filter(Y_best_response_binary == 1)
    survival[["recist_response_adj"]] <- df %>% drop_na(Y_best_response_binary)
    survival[["all_adj"]] <- survival[["all"]]
    
    for (i in treatments) survival[[i]] <- survival[["all"]] %>% filter(Treatment == i)
    for (i in treatments) survival[[paste0(i,"_recist_response_adj")]] <- survival[["all"]] %>% filter(Treatment == i)

    response <- list()
    response[["all"]] <- df %>% filter(Filter_meta_responseMeasured == "Yes") %>% drop_na(Y_best_response_binary)
    response[["all_adj"]] <- response[["all"]]
    for (i in treatments ) response[[i]] <- response[['all']] %>% filter(Treatment == i)
    
    ### Analysis Store 
    cpi_store = list();
    cpi_store[['os']] = survival
    cpi_store[['bor']] = response
    cpi_store
}
scale_the_data2 <- function( df ){
    col_types <- sapply(df, typeof)
    responses <- colnames(df %>% select(Y_best_response_binary, contains("Surv")))
    step1 <- col_types[which(col_types != "character")]
    step2 <- names(step1[-which(names(step1) %in% responses)])
    df <- df %>% mutate_at(vars(all_of(step2)), ~(scale(.) %>% as.vector))
    df
}
Typer2 <- function (i) 
{
    if (i %in% c("Age", "Gender", "PST", "CPI")) {
        "Clinical"
    }
    else if (i %in% c("Prolif")) {
        "RNA: Proliferation"
    }
    else if (i %in% c("T_cell")) {
        "RNA: T-cell"
    }
    else if (i %in% c("TGFB")) {
        "RNA: TGFB"
    }
    else if (i %in% c("CD276")) {
        "RNA: Other"
    }
    else if (i %in% c("WGD", "Ploidy")) {
        "CNV"
    }
    else if (i == "SV_TMB") {
        "SVs"
    }
    else if (i == "HLA_het") {
        "HLA"
    }
    else if (i == "TMB") {
        "Somatic"
    }
    else {
        "Missed"
    }
}

### Name mapping 
name_map <- list()
name_map[['T_cell']] <- "T-cell"
name_map[['Pre_treated']] <- "Pre-treated"
name_map[['HLA_het']] <- "HLA Heterozygosity"
name_map[['SV_TMB']] <- "SV TMB"

### Add size 
add_size <- function (i) 
{
    if (!grepl("Non", i)) {
        s <- ifelse( grepl("response", i), 415, 483)
    }
    else {
        if (grepl("Chemo", i)) {
            s <- ifelse( grepl("response",i), 886, 1165)
        }
        else if (grepl("Hormon", i)) {
            s <- ifelse( grepl("response",i), 199, 284)
        }
        else if (grepl("Targe", i)) {
            s <- ifelse( grepl("response",i), 315, 388)
        }
        else if (grepl("Multip", i)) {
            s <- ifelse( grepl("response",i), 532, 621)
        }
        else if (grepl("All", i)) {
            s <- ifelse( grepl("response",i), 1957, 2499)
        }
        else if (grepl("Recist_response_adj", i)) {
            s <- ifelse( grepl("response",i), 1957, 2499)
        }
        else {
            s <- NA
        }
    }
    paste0("(n = ", s, ")")
}

namer <- function(i){
    if(grepl("All",i)){
        "Overall"
    } else if (grepl("Chemo",i)){
        "Chemotherapy Only"
    } else if (grepl("Hormonal",i)){
        "Hormonal Only"
    } else if (grepl("Multiple",i)){
        "Multiple Therapies"
    } else if (grepl("Immunotherapy",i)){
        "Immunotherapy Only"
    } else if (grepl("Targeted",i)){
        "Targeted Only"
    } else {
        "Overall"
    }
}
