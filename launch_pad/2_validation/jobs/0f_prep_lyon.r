options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))

library(tidyverse)
library(stringr)
library(survival)
library(GEOquery)

I_DIR <- paste0(E_DIR, "lyon/")

gse159067 = getGEO(filename=paste0(I_DIR, "clinical/GSE159067_series_matrix.txt"))
gse161537 = getGEO(filename=paste0(I_DIR, "clinical/GSE161537_series_matrix.txt"))
gse162519 = getGEO(filename=paste0(I_DIR, "clinical/GSE162519_series_matrix.txt"))
gse162520 = getGEO(filename=paste0(I_DIR, "clinical/GSE162520_series_matrix.txt"))

go159067 <- pData(gse159067@phenoData)
go161537 <- pData(gse161537@phenoData)
go162519 <- pData(gse162519@phenoData)
go162520 <- pData(gse162520@phenoData)

#experimentData(gse162520)

cleaner1 <- function(i){
    abc <- strsplit(as.character(i), " ")[[1]]
    abc[length(abc)]
}
binary_response <- function( recist ){
    recist = as.character(recist)
    if( recist %in% c("CR","PR")){
        1
    } else if ( recist %in% c("SD", "PD")){
        0
    } else {
        NA
    }
}

dim(go159067)

get_id <- function(i) strsplit(strsplit(i,"-")[[1]][2], "_")[[1]][1]
go159067$patient_id <- unlist(lapply( as.character(go159067$description), get_id))

go159067$age <- unlist(lapply( go159067$characteristics_ch1.1, cleaner1))
go159067$sex <- unlist(lapply( go159067$characteristics_ch1.2, cleaner1))
go159067$cpi_line <- unlist(lapply( go159067$characteristics_ch1.6, cleaner1))
go159067$recist <- unlist(lapply( go159067$characteristics_ch1.7, cleaner1))
go159067$bor <- unlist(lapply( go159067$recist, binary_response))
go159067$os_month <- as.numeric(unlist(lapply( go159067$characteristics_ch1.8, cleaner1)))
go159067$os_event <- as.numeric(unlist(lapply( go159067$characteristics_ch1.9, cleaner1)))

ready_159067 <- (
    go159067 %>% transmute(
                 patient_id,
                 bor,
                 os_days = os_month*(365/12),
                 os_event,
                 age, 
                 gender = sex, 
                 tissue = "hnc", 
                 tissue_full = "metastatic",
                 pretreat = ifelse(cpi_line > 1, 1, 0),
                 pretreat_comp = NA, 
                 Study = "Lyon",
                 cohort = "159067"
             )
)

dim(go161537)

go161537$patient_id <- unlist(lapply( go161537$characteristics_ch1, cleaner1))
go161537$age <- unlist(lapply( go161537$characteristics_ch1.2, cleaner1))
go161537$sex <- unlist(lapply( go161537$characteristics_ch1.3, cleaner1))
go161537$cpi_line <- unlist(lapply( go161537$characteristics_ch1.4, cleaner1))
go161537$recist <- unlist(lapply( go161537$characteristics_ch1.5, cleaner1))
go161537$bor<- unlist(lapply( go161537$recist, binary_response))
go161537$os_month <- as.numeric(unlist(lapply( go161537$characteristics_ch1.6, cleaner1)))
go161537$os_event <- as.numeric(unlist(lapply( go161537$characteristics_ch1.7, cleaner1)))

ready_161537 <- (
    go161537 %>% transmute(
                 patient_id,
                 bor,
                 os_days = os_month*(365/12),
                 os_event,
                 age, 
                 gender = sex, 
                 tissue = "lung", 
                 tissue_full = "metastatic",
                 pretreat = ifelse(cpi_line > 1, 1, 0),
                 pretreat_comp = NA, 
                 Study = "Lyon",
                 cohort = "161537"
             )
)

dim(go161537)

go162520$patient_id <- unlist(lapply( go162520$title, cleaner1))
go162520$age <- unlist(lapply( go162520$characteristics_ch1.2, cleaner1))
go162520$sex <- unlist(lapply( go162520$characteristics_ch1.3, cleaner1))
go162520$os_month <- as.numeric(unlist(lapply( go162520$characteristics_ch1.4, cleaner1)))
go162520$os_event <- as.numeric(unlist(lapply( go162520$characteristics_ch1.5, cleaner1)))

ready_162520 <- (
    go162520 %>% transmute(
                 patient_id,
                 bor = NA,
                 os_days = os_month*(365/12),
                 os_event,
                 age, 
                 gender = sex, 
                 tissue = "lung", 
                 tissue_full = "non-advanced",
                 pretreat = NA,
                 pretreat_comp = NA, 
                 Study = "Lyon",
                 cohort = "162520"
             )
)

dim(go162519)

go162519$age <- unlist(lapply( go162519$characteristics_ch1.2, cleaner1))
go162519$sex <- unlist(lapply( go162519$characteristics_ch1.3, cleaner1))
go162519$os_month <- as.numeric(unlist(lapply( go162519$characteristics_ch1.21, cleaner1)))
go162519$os_event <- as.numeric(unlist(lapply( go162519$characteristics_ch1.22, cleaner1)))

ready_162519 <- (
    go162519 %>% transmute(
                 patient_id = title,
                 bor = NA,
                 os_days = os_month*(365/12),
                 os_event,
                 age, 
                 gender = sex, 
                 tissue = "hnc", 
                 tissue_full = "non-advanced",
                 pretreat = NA,
                 pretreat_comp = NA, 
                 Study = "Lyon",
                 cohort = "162519"
             )
)

clinical_go <- bind_rows(ready_159067, ready_161537, ready_162520, ready_162519) %>% filter(os_event %in% c(0,1))

sanity <- function(i) as.numeric(gsub(",", ".", as.character(i)))

rna159067 <- read.table( paste0(I_DIR,"rna/GSE159067_IHN_log2cpm_data.txt"), sep = "\t", header = TRUE)
rna161537 <- read.csv( paste0(I_DIR,"rna/GSE161537_nivobio_log2cpm.csv"),  sep = ';', header = TRUE)
rna162519 <- read.csv( paste0(I_DIR,"rna/GSE162519_GEO_data_LBCC1_log2CPM.csv"), sep = ";", header = TRUE)
rna162520 <- read.csv( paste0(I_DIR,"rna/GSE162520_GEO_data_TUMADOR_log2cpm.csv"), sep = ";", header = TRUE)

share <- list()
share[['rna159067_log2']] <- data.frame(t(rna159067 %>% column_to_rownames("ID_REF")))
share[['rna161537_log2']] <- data.frame(t(rna161537 %>% column_to_rownames("Patient_number"))) %>% mutate_all(sanity)
share[['rna162519_log2']] <- data.frame(t(rna162519 %>% column_to_rownames("X"))) %>% mutate_all(sanity)
share[['rna162519_log2']] <- data.frame(t(rna162520 %>% column_to_rownames("X"))) %>% mutate_all(sanity)
saveRDS( share, paste0(REF_DIR, "rna_lyon.Rds") )

rna_ready_159067 <- (
data.frame(
    t( rna159067
         %>% filter(ID_REF %in% unlist(gene_sets)) 
         %>% column_to_rownames("ID_REF") 
      )) %>% rownames_to_column("patient_id")
    %>% select_if( names(.) %in% c("patient_id", unlist(gene_sets)))
)

get_clean_id <- function(i) strsplit(strsplit(i,"\\.")[[1]][2], "_")[[1]][1]
rna_ready_159067$patient_id <- unlist(lapply(rna_ready_159067$patient_id, get_clean_id))

rna_ready_159067$prolif <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters$prolif), 1, mean)
rna_ready_159067$tcell <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters$tcell), 1, mean)
rna_ready_159067$tgfb <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters$tgfb), 1, mean)

rna_ready_159067$prolif_cluster5 <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters5$prolif), 1, mean)
rna_ready_159067$tcell_cluster5 <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters5$tcell), 1, mean)
rna_ready_159067$tgfb_cluster5 <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$clusters5$tgfb), 1, mean)

rna_ready_159067$prolif_set <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$sets1$prolif), 1, mean)
rna_ready_159067$tcell_set <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$sets1$tcell), 1, mean)
rna_ready_159067$tgfb_set <- apply( rna_ready_159067 %>% select_if(names(.) %in% gene_sets$sets1$tgfb), 1, mean)

rna_ready_161537 <- (
data.frame(
    t( rna161537
         %>% filter(Patient_number %in% unlist(gene_sets)) 
         %>% column_to_rownames("Patient_number") 
      )) %>% rownames_to_column("patient_id")
    %>% select_if( names(.) %in% c("patient_id", unlist(gene_sets)))
)

rna_ready_161537 <- rna_ready_161537 %>% mutate_at(vars(-patient_id), sanity)

rna_ready_161537$prolif <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters$prolif), 1, mean)
rna_ready_161537$tcell <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters$tcell), 1, mean)
rna_ready_161537$tgfb <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters$tgfb), 1, mean)

rna_ready_161537$prolif_cluster5 <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters5$prolif), 1, mean)
rna_ready_161537$tcell_cluster5 <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters5$tcell), 1, mean)
rna_ready_161537$tgfb_cluster5 <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$clusters5$tgfb), 1, mean)

rna_ready_161537$prolif_set <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$sets1$prolif), 1, mean)
rna_ready_161537$tcell_set <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$sets1$tcell), 1, mean)
rna_ready_161537$tgfb_set <- apply( rna_ready_161537 %>% select_if(names(.) %in% gene_sets$sets1$tgfb), 1, mean)

rna_ready_161537$patient_id <- unlist(lapply(rna_ready_161537$patient_id, function(i) strsplit(i, "X")[[1]][2]))

rna_ready_162519 <- (
data.frame(
    t( rna162519
         %>% filter(X %in% unlist(gene_sets)) 
         %>% column_to_rownames("X") 
      )) %>% rownames_to_column("patient_id")
    %>% select_if( names(.) %in% c("patient_id", unlist(gene_sets)))
)

rna_ready_162519 <- rna_ready_162519 %>% mutate_at(vars(-patient_id), sanity)

rna_ready_162519$prolif <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters$prolif), 1, mean)
rna_ready_162519$tcell <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters$tcell), 1, mean)
rna_ready_162519$tgfb <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters$tgfb), 1, mean)

rna_ready_162519$prolif_cluster5 <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters5$prolif), 1, mean)
rna_ready_162519$tcell_cluster5 <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters5$tcell), 1, mean)
rna_ready_162519$tgfb_cluster5 <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$clusters5$tgfb), 1, mean)

rna_ready_162519$prolif_set <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$sets1$prolif), 1, mean)
rna_ready_162519$tcell_set <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$sets1$tcell), 1, mean)
rna_ready_162519$tgfb_set <- apply( rna_ready_162519 %>% select_if(names(.) %in% gene_sets$sets1$tgfb), 1, mean)

rna_ready_162520 <- (
data.frame(
    t( rna162520
         %>% filter(X %in% unlist(gene_sets)) 
         %>% column_to_rownames("X") 
      )) %>% rownames_to_column("patient_id")
    %>% select_if( names(.) %in% c("patient_id", unlist(gene_sets)))
)

rna_ready_162520 <- rna_ready_162520 %>% mutate_at(vars(-patient_id), sanity)

rna_ready_162520$prolif <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters$prolif), 1, mean)
rna_ready_162520$tcell <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters$tcell), 1, mean)
rna_ready_162520$tgfb <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters$tgfb), 1, mean)

rna_ready_162520$prolif_cluster5 <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters5$prolif), 1, mean)
rna_ready_162520$tcell_cluster5 <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters5$tcell), 1, mean)
rna_ready_162520$tgfb_cluster5 <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$clusters5$tgfb), 1, mean)

rna_ready_162520$prolif_set <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$sets1$prolif), 1, mean)
rna_ready_162520$tcell_set <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$sets1$tcell), 1, mean)
rna_ready_162520$tgfb_set <- apply( rna_ready_162520 %>% select_if(names(.) %in% gene_sets$sets1$tgfb), 1, mean)

rna_ready <- bind_rows(rna_ready_159067, rna_ready_161537, rna_ready_162519, rna_ready_162520)

#### To be added

ok <- clinical_go %>% inner_join(rna_ready, by = "patient_id")
ok$os <- ifelse( ok$os_event == 0, -ok$os_days, ok$os_days)

lyon_go <- (
    ok 
        %>% transmute(
          patient_id, 
          bor, 
          os = ifelse( os_event == 0, -os_days, os_days), 
          os_event,
          os_days, 
          age, 
          gender, 
          tissue, 
          tissue_full,
          tmb = NA,
          tcell,
          prolif, 
          tgfb,
          tcell_cluster5,
          prolif_cluster5, 
          tgfb_cluster5,
          tcell_set, 
          prolif_set, 
          tgfb_set, 
          pdl1 = CD274,
          pretreat, ## these cohorts are early stage..
          pretreat_comp = NA, 
          purity = NA, 
          Study,
          cohort,
          extra = NA, 
          extra2 = NA
    )
)

saveRDS( lyon_go, paste0( TMP_DIR, "validation-lyon-go.Rds"))
