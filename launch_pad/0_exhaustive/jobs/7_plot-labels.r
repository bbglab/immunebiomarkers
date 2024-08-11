wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_prep.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
library(tidyverse)
library(stringr)

boom <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))

overall <- list(
    'somatic_TMB' = 'TMB',
    'clinical_pre_contains_Chemotherapy' = 'Prior Chemotherapy',
    'clinical_meta_hasRadiotherapyPreTreatment' = 'Prior Radiotherapy',
    'clinical_meta_hasSystemicPreTreatment2' = "Prior systemic therapy",
    'clinical_pre_treated' = "Prior therapy",
    'clinical_systemic_composite' = "Prior Systemic Therapy (Composite)",
    'clinical_pre_to_post_treatment_time' = "Time since Prior Systemic Therapy",
    'isofox_gene_set_prolif' = 'Proliferation',
    'isofox_gene_set_Pan_TBRS' = 'TGFB',
    'isofox_CD276' = 'CD276' ,
    'isofox_gene_set_t_cell_effector' = 'T-cell effector gene set',
    'sv_summary_svTumorMutationalBurden' = 'SV TMB',
    'cnv_scna' = 'SCNA',
    'cnv_summary_wholeGenomeDuplication' = 'WGD',
    'cnv_summary_diploidProportion' = 'Ploidy',
    'hla_HLA_all_LOH' = 'HLA LOH' 
)

tmb <- list(
    "somatic_TMB" = "TMB",
    "somatic_TMB_clonal" = "TMB clonal",
    "somatic_TMB_subclonal" = "TMB subclonal",
    "somatic_TMB_frameshift" = "TMB frameshift")

tmb_genes <- list()
somatic_genes = c('BRCA2', 'MSH2', 'POLD1')
idx <- paste0("somatic.gene_",somatic_genes,".mb")
tmb_genes[idx] <- somatic_genes

sigs <- list( "sig_SBS7a" = "Signature 7a")

t_cells <- list()
genes <- c('CXCL9', 'CD274','CTLA4','TIGIT')
idx <- paste0('isofox_', genes)
t_cells[idx] <- genes

prolif <- list()
genes <- gene_sets[['gene_set_prolif']]
idx <- paste0('isofox_', genes)
prolif[idx] <- genes

tgfb <- list()
genes <- gene_sets[['gene_set_Pan_TBRS']]
idx <- paste0('isofox_', genes)
tgfb[idx] <- genes

cibersort <- list(
    'cibersort_T.cells.CD8' = "CD8 T-cells",
    'cibersort_Dendritic.cells.resting' = "Dendritic Cells Resting",
    'cibersort_Macrophages.M2' = "Macrophages M2",
    'cibersort_Macrophages.M0' = "Macrophages M0",
    'cibersort_T.cells.gamma.delta' = "Gamma T-cells")

pre_set_labels <- c( overall, tmb, tmb_genes, t_cells, prolif, tgfb, cibersort, sigs )
labeller <- function(i) pre_set_labels[[i]]

annotate <- 
boom %>% 
  filter( feature %in% names(pre_set_labels)) %>% 
  rowwise() %>% 
  mutate(clean_label = labeller(feature)) %>%
  ungroup()

saveRDS( annotate, paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
