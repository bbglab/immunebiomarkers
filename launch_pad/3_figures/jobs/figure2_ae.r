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

fill_map <- list("Clinical"="#9E7788","CNV/SVs"="#B3DE69","HLA"="#FFFF99","RNA"="#FC913A","Somatic"="#80B1D3")
size_map <- list( "TRUE" = 4, "FALSE" = 2)
alpha_map <- list( "TRUE" = 1, "FALSE" = .4)
color_map <- list( "TRUE" = "black", "FALSE" = "white")
color_map2 <- list( "In Cluster" = "black", "Not in Cluster" = "white")

base_data <- (
    ingredients 
      %>% filter(
           !grepl("ciber", feature),
           feature != "clinical_systemic_composite",
           dataset == "all", 
           model == "bor",
           covariates == "age_biopsy_purity_tissue"
      )
      %>% mutate(in_cluster = ifelse(Group != "Gene Set" & log10_p > threshold, "In Cluster", "Not in Cluster"))
)

plot_ready <- (
  base_data 
    %>% select(feature, big_group, cor_pretreat, cor_tcell, cor_tmb, plot_est, log10_p, Group, in_cluster)
    %>% gather(plot, val, -feature, -big_group, -log10_p, -Group, -in_cluster)
    %>% drop_na()
    %>% mutate( big_group = factor(big_group, levels = c("Somatic", "RNA", "Clinical", "CNV/SVs", "HLA")))
)

labels_ready <- (
  labels
    %>% filter(clean_label %in% c("TMB", "Prior therapy", "T-cell effector gene set"), dataset == "all", model == "bor", covariates == "age_biopsy_purity_tissue")
    %>% select(feature, big_group, cor_pretreat, cor_tcell, cor_tmb, plot_est, log10_p, clean_label)
    %>% gather(plot, val, -feature, -big_group, -log10_p, -clean_label)
)

plots_go <- function( settings ){
    base <- (
      ggplot( settings$data, settings$aes) + 
        geom_point(shape = 21, stroke = 0.3) + 
        scale_y_continuous(breaks = c(seq(0,8,2)), limits = c(0,8)) + 
        scale_x_continuous(n.breaks = 3)  +
        scale_size_manual(values = unlist(size_map))  +
        scale_alpha_manual(values = unlist(alpha_map)) + 
        scale_color_manual(values = unlist(color_map) ) + 
        scale_fill_manual(values = unlist(fill_map), limits = force) + 
        geom_hline( yintercept = threshold, linetype="dashed", color = "black" ) + 
        xlab( settings$x ) + ylab("-Log10 p-value") +  ggtitle( settings$title) +
        settings$theme
    )
    if( !is.data.frame(settings$labels) ){
        base 
    } else {
        base + geom_label_repel(data = settings$labels, aes(x = val, y = log10_p, label = clean_label), 
                                size = 5, nudge_y = 1,inherit.aes = FALSE)
    }
}

guides_main <- guides(alpha = "none", size = "none", color = "none", fill = guide_legend(override.aes = list(size=5)))
guides_3 <- guides(alpha = "none", size = "none", fill = guide_legend(override.aes = list(size=5), order = 1), color = guide_legend(override.aes = list(size=5), order = 2))

aes_main <- aes(x = val, y = log10_p, fill = big_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_clust <- aes(x = val, y = log10_p, fill = big_group, size = val > .5, alpha = log10_p > threshold, color = in_cluster ) 
aes_clust2 <- aes(x = val, y = log10_p, fill = big_group, size = val > .5, alpha = log10_p > threshold, color = log10_p > 0 ) #+ guides_main

plts <- list( 
    "main" = list("data" = plot_ready %>% filter(plot == "plot_est", Group != "Gene Set"), 
                  "theme" = base_theme,
                  "x" = "Response Odds Ratio Estimate",
                  "title" = "Best Overall Response (BOR) vs All Features",
                  "aes" = aes_main,
                  "labels" = FALSE),
    "clust1" = list("data" = plot_ready %>% filter(plot == "cor_tmb", big_group %in% c("Somatic", "CNV/SVs")),
                    "theme" = base_theme, 
                    "x" = "Correlation to Cluster 1 Mean",
                    "title" = "BOR vs Somatic and SV features",
                    "aes" = aes_clust,
                    "labels" = labels_ready %>% filter(plot == "cor_tmb", clean_label == "TMB")),
    "clust2" = list("data" = plot_ready %>% filter(plot == "cor_pretreat",  big_group %in% c("Clinical", "HLA")),
                    "theme" = base_theme,
                    "x" = "Correlation to Cluster 2 Mean",
                    "title" = "BOR vs Clinical and HLA features", 
                    "aes" = aes_clust2,
                    "labels" = labels_ready %>% filter(plot == "cor_pretreat", clean_label == "Prior therapy")),
    "clust3" = list("data" = plot_ready %>% filter(plot == "cor_tcell", big_group == "RNA"),
                    "theme" =  base_theme,
                    "x" = "Correlation to Cluster 3 Mean",
                    "title" = "BOR vs RNA Expression features",
                    "aes" = aes_clust,
                   "labels" = labels_ready %>% filter(plot == "cor_tcell", clean_label == "T-cell effector gene set"))
)

fig_2a <- plots_go(plts$main) + theme(legend.position = c(0.85, 0.25)) + guides_main  + geom_text(x = 1, y = threshold + .25, label = "BY Threshold", size = 4, color = "black")# + geom_vline( xintercept = 1, linetype="dashed", color = "black" ) 

clust1 <- plots_go(plts$clust1) + theme(legend.position = c(0.17, 0.92)) + guides_main + geom_text(x = -.3, y = threshold + .25, label = "BY Threshold", size = 4, color = "black")
clust2 <- plots_go(plts$clust2) + theme(legend.position = c(0.17, 0.92)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank() ) + guides_main
clust3 <- plots_go(plts$clust3) + theme(legend.position = c(0.17, 0.87)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + guides_3 + scale_color_manual(values = unlist(color_map2) )

fig_2e <- as_ggplot(arrangeGrob( clust1, clust2, clust3, ncol = 3))

saveRDS( list( "a" = fig_2a, "e" = fig_2e), file = paste0(FIG_DIR, "fig_2ae.Rds"))
