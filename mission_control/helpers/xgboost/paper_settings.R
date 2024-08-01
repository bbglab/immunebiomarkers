### MODELS ##

tmb_bin <- c("tmb_bin")
tmb <- c("tmb")
base_model <- c("tmb", "pdl1")
rna_only <- c( "tcell", "prolif", "tgfb")
molecular_only <- c( "tmb", "tcell", "prolif", "tgfb")
main_model <- c( "tmb", "tcell", "prolif", "tgfb", "pretreat")

extra_features <- c( "age",
                     "clinical_biopsy_distal_proximal",
                     "clinical_cpi_mechanism3",
                     "hla_lilac_del_hla",
                     "cnv_summary_wholeGenomeDuplication",
                     "sv_summary_svTumorMutationalBurden")

vhio_model <- c('somatic_TMB_vhio', 
                'isofox_gene_set_vhio_tgfb', 
                'isofox_gene_set_vhio_prolif', 
                'isofox_gene_set_vhio_tcell', 
                'pretreat')

latent1 <- c( "somatic_TMB_clonal", 
              "isofox_gene_set_t_cell_effector", 
              "isofox_gene_set_prolif", 
              "isofox_gene_set_Pan_TBRS",
              "clinical_pre_treated")

latent2 <- c( "somatic_TMB_exome", 
              "isofox_gene_set_vhio_tcell", 
              "isofox_gene_set_mariathan_Cell_cycle", 
              "isofox_gene_set_mariathan_EMT2",
              "clinical_meta_hasSystemicPreTreatment2")

model_features <- list(
    "tmb_bin" = tmb_bin,
    "tmb" = tmb,
    "base" = base_model,
    "rna_only" = rna_only, 
    "molecular_only" = molecular_only,
    "no_tmb" = main_model[-1],
    "five_latent" =  main_model,
    "five_latent_interaction" = main_model, 
    "five_latent_complex" = main_model,
    "full_mod" = c(main_model, extra_features),
    "latent_vhio" = vhio_model,
    "latent1" = latent1,
    "latent2" = latent2
)   

### PARAMETER GRIDS ###
parameter_grid <- expand.grid(
    eta = c(.05),
    max_depth = c(1),
    min_child_weight = c(5),
    subsample = c(.75), 
    colsample_bytree = c(1),
    optimal_trees = 0,                
    error = 0
)

parameter_grid_int <- expand.grid(
    eta = c(.05),
    max_depth = c(1,2),
    min_child_weight = c(5),
    subsample = c(.75),
    colsample_bytree = c(1),
    optimal_trees = 0,
    error = 0
)

parameter_grid_comp <- expand.grid(
    eta = c(.05, .1),
    max_depth = c(1,2),
    min_child_weight = c(5,10),
    subsample = c(.75,1),
    colsample_bytree = c(1),
    optimal_trees = 0,
    error = 0
)

get_parameter_grid <- function(i){
      if( grepl("complex",i)){
          parameter_grid_comp
      } else if (grepl("interaction",i)){
          parameter_grid_int
      } else {
          parameter_grid
      }
}
