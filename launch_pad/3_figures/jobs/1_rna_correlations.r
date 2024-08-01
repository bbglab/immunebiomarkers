#wd <- dirname(dirname(getwd()))
#source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)

REF_DIR <- "/workspace/projects/immune_biomarkers/repo/immune_biomarkers/ref/"
setwd(REF_DIR)

hmf <- readRDS('rna_hmf.Rds')   ### log-transformed
inspire <- readRDS('rna_inspire.Rds') ### some transformation
mariathasan <- readRDS('rna_mariathasan.Rds')  ### raw
ravi <- readRDS('rna_ravi.Rds')  ### ravi
lyon <- readRDS('rna_lyon.Rds')  ### log2 transformations
parker <- readRDS('rna_parker.Rds') ### raw

gene_sets = list(
    'gene_set_prolif' = c('BUB1','CCNB2','CDK1','CDKN3','FOXM1','KIAA0101','MAD2L1','MELK','MKI67','TOP2A'),
    'gene_set_t_cell_effector' = c('CD8A','CD27','IFNG','GZMA','GZMB','PRF1','EOMES','CXCL9','CXCL10','CXCL11','CD274','CTLA4','FOXP3','TIGIT','IDO1','PSMB8','PSMB9','TAP1','TAP2'),
    'gene_set_Pan_TBRS' = c('TGFB1','TGFBR2','ACTA2','COL4A1','TAGLN','SH3PXD2A')
)

#library(corrplot)
#corrplot(cor(mariathasan %>% select(all_of(unlist(gene_sets)))))
