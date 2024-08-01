options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))

library(dplyr)
library(stringr)
library(tidyverse)

I_DIR <- paste0(E_DIR, "/mariathan/")

oo <- readRDS( paste0(I_DIR,"cds_extract.RData"))
ff <- readRDS( paste0(I_DIR, "fmone_extract.RData"))

clinical <- oo$pData_df

binary_response <- function( recist ){
    if( recist %in% c("CR","PR")){
        1
    } else if ( recist %in% c("SD", "PD")){
        0
    } else {
        NA
    }
}

clinical <- tibble::rownames_to_column(clinical, "patient_id")
clinical$bor <- unlist(lapply( as.character(clinical$"Best Confirmed Overall Response"), binary_response))
clinical$ecog <- as.numeric(clinical$"Baseline ECOG Score" > 0)
#clinical$pretreat <- ifelse( clinical$"Sample collected pre-platinum" == "N", 1, 0 )
clinical$os_days <- clinical$os * 365/12
clinical$os_event <- ifelse(clinical$censOS == 1, 1, 0)
clinical$os <- ifelse(clinical$os_event == 1, clinical$os_days, -clinical$os_days)
clinical$tmb <- clinical$"FMOne mutation burden per MB"
clinical$ip <- clinical$"Immune phenotype"
clinical$plat <- clinical$"Received platinum"
clinical$sample_age <- clinical$"Sample age"

clinical_go <- (
    clinical
        %>% transmute(
          patient_id, 
          bor,
          os, 
          os_event, 
          os_days, 
          age = NA, 
          gender = Sex,
          tissue = "bladder", 
          tissue_full = Tissue,
          tmb = log(tmb+1),
          pretreat = NA,
          pretreat_comp = NA,
          purity = ip,
          Study = "MARIATHAN",
          extra = ecog, 
          extra2 = sample_age
        )
)

df_tmp <- as.data.frame(t(oo$cds))

share <- data.frame(
    t(data.frame(t(df_tmp)) 
        %>% rownames_to_column("entrez_id")
        %>% inner_join( oo$fData %>% select(entrez_id, Symbol), by = "entrez_id")
        %>% filter(! Symbol %in% c("", "CSNK1E")) %>% drop_na(Symbol)
        %>% select(-entrez_id)
        %>% column_to_rownames("Symbol")
    )    
)
saveRDS( share, paste0(REF_DIR, "rna_mariathasan.Rds"))

genes <- unlist(gene_sets)
symbols <- oo$fData[which(oo$fData$Symbol %in% genes),]
df_tmp<- as.data.frame(oo$cds)[which(rownames(as.data.frame(oo$cds))%in% symbols$entrez_id),]
rna <- as.data.frame(t(df_tmp))

colnames(rna) <- oo$fData %>% filter(entrez_id %in% colnames(rna)) %>% pull(symbol)
rna <- rna %>% mutate_at(colnames(rna), ~(log(.+1) %>% as.vector))
rna <- tibble::rownames_to_column(rna, "patient_id")

names(gene_sets)

rna$tcell <-  apply(rna %>% select( any_of(gene_sets$clusters$tcell)),1,mean)  
rna$prolif <- apply(rna %>% select( any_of(gene_sets$clusters$prolif)),1,mean) 
rna$tgfb <-   apply(rna %>% select( any_of(gene_sets$clusters$tgfb)),1,mean) 

rna$tcell_cluster5 <-  apply(rna %>% select( any_of(gene_sets$clusters5$tcell)),1,mean)  
rna$prolif_cluster5 <- apply(rna %>% select( any_of(gene_sets$clusters5$prolif)),1,mean) 
rna$tgfb_cluster5 <-   apply(rna %>% select( any_of(gene_sets$clusters5$tgfb)),1,mean) 

rna$tcell_set <-  apply(rna %>% select( any_of(gene_sets$sets1$tcell)),1,mean)  
rna$prolif_set <- apply(rna %>% select( any_of(gene_sets$sets1$prolif)),1,mean) 
rna$tgfb_set <-   apply(rna %>% select( any_of(gene_sets$sets1$tgfb)),1,mean) 

rna$pdl1 <-   apply(rna %>% select( CD274),1,mean) 

rna_go <- rna %>% select( patient_id, tcell, prolif, tgfb, tcell_cluster5, prolif_cluster5, tgfb_cluster5, tcell_set, prolif_set, tgfb_set, pdl1)

mariathan_go <- clinical_go %>% left_join(rna_go, by = "patient_id")

saveRDS( mariathan_go, paste0( TMP_DIR, "validation-mariathan-go.Rds"))
