wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/validation/settings.R"))
library(tidyverse)

load(paste0( E_DIR, "inspire/inspire-genomics/data/RData/PanCanMut.data.RData"))

rna <- read.table( paste0(E_DIR, "inspire/inspire-genomics/data/Source Data/SourceData_Fig4/gene-expression-matrix-TPM-final.tsv")) 

pt.dat$bor <- pt.dat[,"PR/CR"]
pt.dat$clinical_recist <- pt.dat[,"Best overall response"]
pt.dat$benefit <- pt.dat[,"Clinical Benefit"]

clinical <- (pt.dat 
     %>% rownames_to_column("patient_id")
     %>% transmute( 
         patient_id, 
         gender = SEX, 
         age = AGE, 
         tissue = COHORT, 
         tissue_full = COHORT,
         pretreat = PRIOR_SYSTEMIC_THERAPY, 
         bor,
         pfs_event = PFS_EVENT, 
         pfs_days = PFS, 
         os_event = OS_EVENT, 
         os_days = OS,
         clinical_recist,
         benefit, 
         biopsy_location = NA,
         mechanism = "pd"
     )
)

mut <- suv.dat %>% transmute(
     patient_id = pt,
     tmb_raw = ns.per.Mb, 
     tmb = log(ns.per.Mb + 1),
     OS, PFS, OS_EVENT, PFS_EVENT, response
 )

clin_mut <- clinical %>% left_join(mut, by = "patient_id")

tmp <- data.frame(t(rna %>% rownames_to_column("gene")))
rna_send <- tmp[-1,]
names(rna_send) <- unname(unlist(lapply(tmp[1,], as.character)))
rna_send <- log(data.frame(lapply(rna_send, as.numeric)))
saveRDS( rna_send, paste0( TMP_DIR, "rna_validation_inspire.Rds"))

express <- (
    rna
        %>% rownames_to_column("gene")
        %>% filter( gene %in% unlist(gene_sets))
)

rna <- data.frame(t(express)[-1,])
colnames(rna) <- t(express)[1,]
patients <- rna %>% rownames_to_column("patient_id") %>% pull(patient_id)
rna <- log(data.frame(lapply(rna, as.numeric)))
rna$patient_id <- patients 

rna$id <- unlist(lapply( strsplit(rna$patient_id, "\\."), function(i) paste( i[c(2,3)], collapse = "-") ))
rna_mns <- rna %>% select(-patient_id) %>% group_by(id) %>% summarise_all( mean)

clin_mut$id <- unlist(lapply( strsplit(clin_mut$patient_id, "-"), function(i) paste( i[c(2,3)], collapse = "-") ))
together <- clin_mut %>% inner_join( rna_mns, by = "id") %>% select(-contains("patient"))

names(gene_sets)

together$tcell <- apply( together %>% select( any_of( gene_sets$clusters$tcell ) ), 1, mean, na.rm = TRUE)
together$tgfb <- apply( together %>% select( any_of( gene_sets$clusters$tgfb ) ), 1, mean, na.rm = TRUE)
together$prolif <- apply( together %>% select( any_of( gene_sets$clusters$prolif) ), 1, mean, na.rm = TRUE)

together$tcell_cluster5 <- apply( together %>% select( any_of( gene_sets$clusters5$tcell ) ), 1, mean, na.rm = TRUE)
together$tgfb_cluster5 <- apply( together %>% select( any_of( gene_sets$clusters5$tgfb ) ), 1, mean, na.rm = TRUE)
together$prolif_cluster5 <- apply( together %>% select( any_of( gene_sets$clusters5$prolif) ), 1, mean, na.rm = TRUE)

together$tcell_set <- apply( together %>% select( any_of( gene_sets$sets1$tcell ) ), 1, mean, na.rm = TRUE)
together$tgfb_set <- apply( together %>% select( any_of( gene_sets$sets1$tgfb ) ), 1, mean, na.rm = TRUE)
together$prolif_set <- apply( together %>% select( any_of( gene_sets$sets1$prolif) ), 1, mean, na.rm = TRUE)

inspire_go <- (
    together 
        %>% transmute(
          patient_id = id, 
          bor, 
          os = ifelse( os_event == 0, -os_days, os_days), 
          os_event,
          os_days, 
          age, 
          gender, 
          tissue, 
          tissue_full = tissue,
          tmb,
          tcell, 
          prolif, 
          tgfb,
          tcell_cluster5,
          tgfb_cluster5,
          prolif_cluster5,  
          tcell_set,
          prolif_set, 
          tgfb_set,
          pdl1 = CD274,
          pretreat,
          pretreat_comp = NA, 
          purity = NA,
          Study = "INSPIRE"
    )
)

saveRDS( inspire_go, paste0( TMP_DIR, "validation-inspire-go.Rds"))
