options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))

library(tidyverse)
library(stringr)

I_DIR <- paste0(E_DIR, "parker_melanoma/")
clinical <- read.csv( paste0(I_DIR, "Clinical/Subjects-CancerCell-MORRISON1-metadata.tsv"), sep = "\t")
somatic <- read.csv( paste0(I_DIR, "WES/WES-CancerCell-MORRISON1-metadata.tsv"), sep = "\t") 
rna <- read.csv( paste0(I_DIR, "RNASeq/data/RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv"), sep = "\t")
#rna <- read.csv( paste0(I_DIR, "RNASeq/data/RNA-CancerCell-MORRISON1-no_batch_correction-logcpm-all_samples.tsv"), sep = "\t")
rna_link <- read.csv( paste0(I_DIR,"RNASeq/RNA-CancerCell-MORRISON1-metadata.tsv"), sep = "\t") #%>% transmute(sample.id, patient_id = subject.id)

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

clinical$bor <- unlist(lapply( as.character(clinical$"bor"), binary_response))

clinical_go <- (
    clinical
        %>% filter(! treatment.regimen.name %in% c("CTLA4", "PD1-to-CTLA4"),
                   !cohort %in% c("064", "va") ### stated not pd1 in paper
                  )
        %>% transmute(
          patient_id = subject.id, 
          bor,
          os_days = os,   
          os = NA, 
          os_event = NA, 
          age = subject.age, 
          gender = subject.sex,
          tissue = "skin", 
          tissue_full = sample.tumor.type,
          pretreat = ifelse(previous.treatment == "naive", 0, 1),
          pretreat_comp = NA,
          Study = "PARKER",
          cohort = cohort,
          extra = treatment.regimen.name, 
          extra2 = meddra.disease.preferred.name
        )
)

somatic_go <- (
    somatic
        %>% filter(! treatment.regimen.name %in% c("CTLA4", "PD1-to-CTLA4", "CTLA4-to-PD1"), #, "PD1-to-CTLA4", "CTLA4-to-PD1"
                  grepl("pre",timepoint.id)
                   )
        %>% drop_na(tmb)
        %>% transmute(
            patient_id = subject.id,
            tmb = log(tmb + 1), 
            purity
        )
)

share <- data.frame(t(rna %>% column_to_rownames("gene.hgnc.symbol")))
saveRDS( share, paste0(TMP_DIR, "rna_parker.Rds") )

rna_link <- rna_link %>% transmute(sample.id, patient_id = subject.id)

rna_clean <- (
    rownames_to_column(
        data.frame(
            t(column_to_rownames(data.frame(rna %>% filter( gene.hgnc.symbol %in% unlist(gene_sets) )), "gene.hgnc.symbol")
        )
    ), "sample.id")
)

rna_clean$tcell  <- apply(rna_clean %>% select( any_of(gene_sets$clusters$tcell)),1,mean)
rna_clean$prolif <- apply(rna_clean %>% select( any_of(gene_sets$clusters$prolif)),1,mean)
rna_clean$tgfb   <- apply(rna_clean %>% select( any_of(gene_sets$clusters$tgfb)),1,mean)

rna_clean$tcell_cluster5  <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$tcell)),1,mean)
rna_clean$prolif_cluster5 <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$prolif)),1,mean)
rna_clean$tgfb_cluster5   <- apply(rna_clean %>% select( any_of(gene_sets$clusters5$tgfb)),1,mean)

rna_clean$tcell_set  <- apply(rna_clean %>% select( any_of(gene_sets$sets1$tcell)),1,mean)
rna_clean$prolif_set <- apply(rna_clean %>% select( any_of(gene_sets$sets1$prolif)),1,mean)
rna_clean$tgfb_set   <- apply(rna_clean %>% select( any_of(gene_sets$sets1$tgfb)),1,mean)

rna_clean$pdl1   <- apply(rna_clean %>% select( CD274 ),1,mean)

rna_go <- rna_clean %>% select(sample.id, tcell, prolif, tgfb, tcell_cluster5, prolif_cluster5, tgfb_cluster5, tcell_set, prolif_set, tgfb_set, pdl1)

parker_go <- (
    clinical_go 
        %>% left_join( somatic_go, by , by = "patient_id") 
        %>% left_join( rna_link, by = "patient_id") 
        %>% left_join( rna_go, by = "sample.id")
        %>% filter( !grepl("_on_", sample.id),  !grepl("_ON_", sample.id), !grepl("EDT", sample.id)) 
)

saveRDS( parker_go, paste0( TMP_DIR, "validation-parker-go.Rds"))
