wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/general.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_plots.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
source(paste0(wd,"/mission_control/helpers/figures/themes.R"))

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
labels <- readRDS(paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
labels$size <- labels$size
threshold <- -log10(ingredients$by_05_fdr)[1]

#o_dir <- paste0(FIG_DIR ,"supplement_note/")
#annote <- function( i, lab, size = 20 ) annotate_figure( i, fig.lab = lab, fig.lab.size = size, fig.lab.face = "bold")

fill_map <- list("Clinical"="#9E7788","CNV/SVs"="#B3DE69","HLA"="#FFFF99","RNA"="#FC913A","Somatic"="#80B1D3")
size_map <- list( "TRUE" = 4, "FALSE" = 2)
alpha_map <- list( "TRUE" = 1, "FALSE" = .4)
color_map <- list( "TRUE" = "black", "FALSE" = "white")
color_map2 <- list( "In Cluster" = "black", "Not in Cluster" = "white")

base_data <- ingredients %>% filter(dataset == "all")
aes_main <- aes(x = plot_est, y = log10_p, fill = Type_plus_main, size = Type_plus_main, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main

#table(base_data$Type, base_data$big_group)

plots_go <- function( settings ){(
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
}

table(base_data$covariates)

settings <- list(
    "bor" = list( "none" = )
    "pfs" = 
    "os" = 
    "residuals" = 
)

base <- annote(latent_plot(type = "base", k = 1, "All features"), "A")
latent <- annote(latent_plot(type = "latent", k = .7, "Features with > .7 correlation to latent factors"), "B")
non_latent <- annote(latent_plot(type = "non_latent", k = .3, "Features with < .3 correlation to latent factors"), "C")
together <- arrangeGrob( base, latent, non_latent, ncol = 1)
go <- get_dressed(together, "HMF Exhaustive Analysis - Latent factors", size = 22, vjust = 2)

options(repr.plot.width = 17, repr.plot.height = 17, repr.plot.res = 100)
plot(go)
ggsave( paste0(o_dir, "sn_exhaustive_latent.png"), go, width = 17, height = 17, dpi = 300)
