wd <- dirname(dirname(getwd()))
setwd(wd)

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure2_theme.R"))

suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
labels <- readRDS(paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
threshold <- -log10(ingredients$by_05_fdr)[1]

geneset_names = data.table::fread(paste0(REF_DIR, "CPI_genesets.txt"))
colnames(geneset_names)[1] = "feature"

geneset_names %>% filter(feature == "mariathan_Cell_cycle")

#fill_map <- list("Clinical"="#9E7788","CNV/SVs"="#B3DE69","HLA"="#FFFF99","RNA"="#FC913A","Somatic"="#80B1D3")
fill_map <- list("Clinical"="#9E7788","CNV/SVs"="#B3DE69","HLA"="#FFFF99","RNA"="grey","Somatic"="#80B1D3")
fill2_map <- list( "Cluster 1 Gene Set" = "#8DD3C7", "Cluster 2 Gene Set" = "#BEBADA","Gene Set" = "black","RNA" = "grey")
size_map <- list( "TRUE" = 4, "FALSE" = 2)
alpha_map <- list( "TRUE" = 1, "FALSE" = .4)
color_map <- list( "TRUE" = "black", "FALSE" = "white")
color_map2 <- list( "In Cluster" = "black", "Not in Cluster" = "white")

names_map <- list("isofox_gene_set_mariathan_Cell_cycle" = "Cell Cycle",
                  "isofox_gene_set_HALLMARK_G2M_CHECKPOINT" = "G2M Checkpoint",
                  "isofox_gene_set_prolif" = "Proliferation Potential",
                  "isofox_gene_set_HALLMARK_ANGIOGENESIS" = "Angiogenesis",
                  "isofox_gene_set_KEGG_ECM_RECEPTOR_INTERACTION" = "ECM receptor interaction", 
                  "isofox_gene_set_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
                  "isofox_gene_set_KEGG_FOCAL_ADHESION" = "Focal Adhesion",
                  "isofox_gene_set_Pan_TBRS" = "TGFB Response signature")
namer <- function(i) { if( i %in% names(names_map)){ names_map[[i]]} else {""}}

typer <- function( big_group, Group, cor_tgfb, cor_prolif ){
    if( Group == "Gene Set"){
        if( cor_tgfb > .8 ){ "Cluster 2 Gene Set"} 
        else if ( cor_prolif > .8 ) {"Cluster 1 Gene Set"}
        else {"Gene Set"}
    } else {
        big_group
    }
}

base_data <- (
    ingredients 
      %>% filter(
           !grepl("ciber", feature),
           feature != "clinical_systemic_composite",
           dataset == "all", 
           model == "os",
           covariates == "residuals",
           !grepl("cluster", feature),
           !grepl("vhio", feature),
           !grepl("rand", feature),
           ! feature %in% c("isofox_gene_set_tcell_cluster", "tcell", "tmb", "tgfb","prolif", "isofox_gene_set_prolif_cluster", "isofox_gene_set_tgfb_cluster"),
           cor_tmb < .3, cor_tcell < .3, cor_pretreat < .3
      )
      %>% mutate(in_cluster = ifelse(Group != "Gene Set" & log10_p > threshold, "In Cluster", "Not in Cluster"))
      %>% rowwise()
      %>% mutate( gene_set_type = typer(big_group, Group, cor_tgfb, cor_prolif))
      %>% mutate( highlight = namer(feature) )
      %>% mutate( big_group = factor(big_group, levels = c("RNA", "Somatic", "Clinical", "CNV/SVs", "HLA")))
)

plots_go <- function( settings ){(
      ggplot( settings$data, settings$aes) + 
        geom_point(shape = 21, stroke = 0.3) + 
        scale_x_continuous(n.breaks = 3)  +
        scale_y_continuous(breaks = c(seq(0,8,2)), limits = c(0,8)) + 
        scale_size_manual(values = unlist(size_map))  +
        scale_alpha_manual(values = unlist(alpha_map)) + 
        scale_color_manual(values = unlist(color_map) ) + 
        scale_fill_manual(values = unlist(fill_map), limits = force) + 
        geom_hline( yintercept = threshold, linetype="dashed", color = "black" ) + 
        xlab( settings$x ) + ylab("-Log10 p-value") +  ggtitle( settings$title) +
        settings$theme
)}

guides_main <- guides(alpha = "none", size = "none", color = "none", fill = guide_legend(override.aes = list(size=5)))

aes_main <- aes(x = plot_est, y = log10_p, fill = big_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_tgfb <- aes(x = cor_tgfb, y = log10_p, fill = gene_set_type, size = cor_tgfb > .5 , alpha =  cor_tgfb > .5, color = log10_p > threshold ) #+ guides_main
aes_prolif <- aes(x = cor_prolif, y = log10_p, fill = gene_set_type, size = cor_prolif > .5, alpha =  cor_prolif > .5, color = log10_p > threshold ) #+ guides_main

plts <- list( 
    "main" = list("data" = base_data %>% filter(Group != "Gene Set"), 
                  "theme" = base_theme,
                  "x" = "1 / OS Hazard Estimate",
                  "title" = "Overall Survival (OS) Residuals vs All Features",
                  "aes" = aes_main,
                  "labels" = FALSE),
    "remaining" = list("data" = base_data %>% filter(cor_tgfb < .3, cor_prolif < .3, cor_tmb < .3, cor_tcell < .3, cor_pretreat < .3), 
                  "theme" = base_theme,
                  "x" = "1 / OS Hazard Estimate",
                  "title" = "Overall Survival (OS) Residuals vs Remaining Features",
                  "aes" = aes_main,
                  "labels" = FALSE),
    "clust1" = list("data" = base_data %>% filter(big_group %in% c("RNA")),
                    "theme" = base_theme, 
                    "x" = "Correlation to Cluster 1 Mean",
                    "title" = "OS Residuals vs RNA Features",
                    "aes" = aes_prolif,
                    "labels" = FALSE),
    "clust2" = list("data" = base_data %>% filter(big_group %in% c("RNA")),
                    "theme" = base_theme, 
                    "x" = "Correlation to Cluster 2 Mean",
                    "title" = "OS Residuals vs RNA Features",
                    "aes" = aes_tgfb,
                    "labels" = FALSE))

fig_3a <- plots_go(plts$main) + theme(legend.position = c(0.8, 0.85)) + guides_main  + geom_text(x = 1, y = threshold + .25, label = "BY Threshold", size = 4, color = "black")

fig_3d <- plots_go(plts$remaining) + theme(legend.position = c(0.8, 0.85)) + guides_main  + geom_text(x = 1, y = threshold + .25, label = "BY Threshold", size = 4, color = "black")

clust1 <- (
  plots_go(plts$clust1) 
    + scale_fill_manual(values = unlist(fill2_map), limits = force) 
    + geom_label_repel(data = base_data %>% filter(gene_set_type == "Cluster 1 Gene Set"), inherit.aes = FALSE,
                     aes(x = cor_prolif, y = log10_p, label = highlight), size = 4, nudge_y = 1, fill = fill2_map$`Cluster 1 Gene Set`)
    #+ theme(legend.position = c(0.15, 0.85)) + guides_main
)

clust2 <- (
  plots_go(plts$clust2) 
    + scale_fill_manual(values = unlist(fill2_map), limits = force) 
    + geom_label_repel(data = base_data %>% filter(gene_set_type == "Cluster 2 Gene Set"), inherit.aes = FALSE,
                     aes(x = cor_tgfb, y = log10_p, label = highlight), size = 3, nudge_y = 1, fill = fill2_map$`Cluster 2 Gene Set`)
    + theme(axis.title.y = element_blank(), axis.text.y = element_blank() )
    + theme(legend.position = c(0.2, 0.85)) + guides_main
)

fig_3c <- as_ggplot(arrangeGrob( clust1, clust2, ncol = 2))

saveRDS( list( "a" = fig_3a, "c" = fig_3c, "d" = fig_3d), file = paste0(FIG_DIR, "fig_3acd.Rds"))
