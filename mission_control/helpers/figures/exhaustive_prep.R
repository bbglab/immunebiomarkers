library(ggplot2)

overall <- list(
    'somatic_TMB' = 'TMB',
    'clinical_age_at_treatment_start' = 'Age',
    'clinical_meta_gender' = 'Gender',
    'clinical_pre_treated' = 'Pre-treated',
    'isofox_gene_set_prolif' = 'Proliferation', 
    'isofox_gene_set_Pan_TBRS' = 'TGFB',
    'isofox_CD276' = 'CD276' ,
    'isofox_gene_set_t_cell_effector' = 'T-cell',
    'sv_summary_svTumorMutationalBurden' = 'SV TMB',
    'cnv_scna' = 'SCNA',
    'cnv_summary_wholeGenomeDuplication' = 'WGD', 
    'cnv_summary_diploidProportion' = 'Ploidy',
    'hla_HLA_all_LOH' = 'HLA LOH' ,
    'hla_HLA_all_heterozygous' = 'HLA Heterozygosity'
)
tmb <- list(                
    "somatic_TMB" = "TMB",
    "somatic_TMB_clonal" = "TMB clonal",
    "somatic_TMB_subclonal" = "TMB subclonal",
    "somatic_TMB_frameshift" = "TMB frameshift",
    "somatic_TMB_indel" = "TMB indel",
    "somatic_TMB_sbs" = "TMB sbs",
    "somatic_TMB_dbs" = "TMB dbs",
    "somatic_TMB_mbs" = "TMB mbs",
    "somatic_summary_msIndelsPerMb" = "TMB MS-indels",
    "somatic_summary_msStatus" = "MSI status"
)
somatic_genes <- c( 'BCLAF1','SERPINB3','SERPINB4','B2M','JAK1','KRAS','TP53','PTEN','RTK','STK11','BAP1','BRCA1','BRCA2','ATM','POLE','ERCC2','FANCA','MSH2','MLH1','POLD1','MSH6','BRAF','TP53')

cibersort <- list(
    'cibersort_T.cells.CD8' = "CD8 T-cells",
    'cibersort_Macrophages.M1' = "Macrophages M1",
    'cibersort_Dendritic.cells.resting' = "Dendritic Cells Resting",
    'cibersort_Macrophages.M2' = "Macrophages M2",
    'cibersort_Mast.cells.activated' = "Mast Cells Activated",
    'cibersort_T.cells.gamma.delta' = "Gamma T-cells",
    'cibersort_NK.cells.activated' = "NK Cell Activated"
)
gene_sets = list(
    'gene_set_cyt' = c('GZMA', 'PRF1'),
    'gene_set_t_cell_gep_6' = c('IFNG','STAT1','CXCL9','CXCL10','IDO1', 'HLA-DRA'),
    'gene_set_t_cell_gep_10' = c('GZMA', 'PRF1', 'IFNG','STAT1','CCR5','CXCL9','CXCL10','CXCL11','IDO1','HLA-DRA'),
    'gene_set_t_cell_gep_18' =  c('CD3D','IDO1','CIITA','CD3E','CCL5','GZMK','CD2','HLA-DRA','CXCL13','IL2RG','NKG7','HLA-E','CXCR6','LAG3','TAGAP','CXCL10','STAT1','GZMA','GZMB'),
    'gene_set_prolif' = c('BUB1','CCNB2','CDK1','CDKN3','FOXM1','KIAA0101','MAD2L1','MELK','MKI67','TOP2A'),
    'gene_set_tim3' = c('HAVCR2', 'LAG3', 'PDCD1', 'CTLA4', 'C10orf54', 'BTLA', 'FOXP3'),
    'gene_set_t_cell_effector' = c('CD8A','CD27','IFNG','GZMA','GZMB','PRF1','EOMES','CXCL9','CXCL10','CXCL11','CD274','CTLA4','FOXP3','TIGIT','IDO1','PSMB8','PSMB9','TAP1','TAP2'),
    'gene_set_myeloid_inflammation' = c('CXCL1','CXCL2','CXCL3','IL8','IL6','PTGS2'),
    'gene_set_stroma_emt_shortened' = c('FLNA','EMP3','CALD1','FN1','FOXC2','LOX','FBN1','TNC'),
    'gene_set_Pan_TBRS' = c('TGFB1','TGFBR2','ACTA2','COL4A1','TAGLN','SH3PXD2A'),
    'gene_set_impres' = c('PDCD1','CD27','CTLA4','CD40','CD86','CD28','CD80','TNFRSF14','TNFSF4','TNFRSF9','C10orf54','HAVCR2','CD200','CD276','CD274'), ### based on pairwise comparisons
    'gene_set_12_chemokine' = c('CCL2', 'CCL3', 'CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13'), ### usa PCA to combine (1st PC)
    'gene_set_immune_checkpoint_genes' = c('CD8A','PDCD1','CD274','PDCD1LG2','CTLA4','CD80','CD86','BTLA','TNFRSF14','LAG3'), ### from herv-3 paper
    'gene_set_cd8_t_effector' = c('CD8A','GZMA','GZMB','IFNG', 'CXCL9','CXCL10','PRF1','TBX21'),
    'gene_set_infiltrate'  = c('CD247', 'CD2', 'CD3E', 'GZMH', 'NKG7', 'PRF1','GZMK')
)
sort_the_drivers <- function(i){
  if( grepl("driver", i[1])){
	ifelse(grepl("MUTATION",i[1]), "driver_somatic", "driver_cnv")
  } else {
	i[2]
  }
}
big_grouper <- function(i){
    if ( i %in% c( "somatic.gene.pc" , "somatic.pc" , "somatic.gene" , "somatic", "sig", "driver_somatic") ){
        "Somatic"
    } else if ( i %in% c("hla")){
        "HLA"
    } else if ( i %in% c("sv") ) {
        "SVs"
    } else if ( i %in% c( "isofox", "cibersort", "isofox.pc" ) ){
        "RNA"
    } else if ( i %in% c( "cnv", "cnv.region","driver_cnv") ){
        "CNV"
    } else {
        "Clinical"
    } 
}
little_grouper <- function(i){
    if ( i %in% c( "somatic.gene.pc" , "somatic.pc" , "somatic.gene" , "somatic") ){
        "Somatic"
    } else if ( i %in% c( "isofox", "isofox.pc" ) ){
        "isofox"
    } else if ( i %in% c( "cnv", "cnv.region" ) ){
        "CNV"
    } else {
        i
    } 
}
cafe_con_leche <- function (cafe, k) {
    cortado <- list(cafe)[[1]]
    big_group <- cortado[2]
    if (big_group %in% c("Clinical", "HLA")) {
        big_group
    } else if ( big_group == "Somatic" ){
        "Somatic: Mutation"
    } else if ( big_group == "CNV" ){
        "Somatic: CNV"
    } else if ( big_group == "SVs" ){
        "Somatic: SV"
    } else {
        tmp <- list(`RNA: T-cell` = as.numeric(cortado[3]), 
                    `RNA: Proliferation` = as.numeric(cortado[4]), 
                    `RNA: TGFB` = as.numeric(cortado[5]))
        max_cor <- max(unlist(tmp))
        ifelse(max_cor > k, names(which(tmp == max_cor))[1], "RNA: Remaining")
    }
}
grouper <- function(k){
    j <- k[1]
    i <- k[2]   
    if( grepl("gene_set", j)){
        "Gene Set"
    } else if (i == "somatic") {
        "TMB Variant"
    } else if (i == "somatic.gene"){
        "Gene Mutations"
    } else if (i == "driver_somatic"){
        "Driver"
    } else if (i == "sig"){
        "Signature Exposure"	
    } else if (i == "isofox"){
        "Gene"
    } else {
        str_to_title(i)
    }
}
Typer2 <- function(i){  
    if( i %in% c("Age", "Gender", "Pre_treated")){
        "Clinical"
    } else if (i %in% c("Proliferation")){
        "RNA: Proliferation"
    } else if (i %in% c("T_cell")) {
        "RNA: T-cell"
    } else if (i %in% c("TGFB", "CD276")) {
        "RNA: CD276/TGFB"
    } else if (i %in% c("WGD", "Ploidy")){
        "CNV"
    } else if (i == "SV_TMB"){
        "SVs"
    } else if (i == "HLA_LOH"){
        "HLA"
    } else if (i == "TMB"){
        "Somatic"
    } else {
        "Missed"
    }
}