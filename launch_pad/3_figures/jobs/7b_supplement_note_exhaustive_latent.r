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

o_dir <- paste0(FIG_DIR ,"supplement_note/")
annote <- function( i, lab, size = 20 ) annotate_figure( i, fig.lab = lab, fig.lab.size = size, fig.lab.face = "bold")

type <- "base"

latent_plot <- function( type = "base", k, title ){

    if(type != "base"){
        labels <- labels %>% filter(clean_label == "dont_show")
    }
    
    if(type == "base"){
        ready <- ingredients
    } else if (type == "latent"){
        ready <- ingredients %>% filter( abs(cor_tmb) > k | abs(cor_tcell) > k | abs(cor_tgfb) > k | abs(cor_prolif) > k | abs(cor_pretreat) > k)
    } else {
        ready <- ingredients %>% filter( abs(cor_tmb) < k, abs(cor_tcell) < k, abs(cor_tgfb) < k, abs(cor_prolif) < k, abs(cor_pretreat) < k)
    }

    m0 <- (
        super_plot(
            ready , 
            labels %>% filter(clean_label %in% c("Prior Systemic Therapy", "TMB", "T-cell")), 
            type = "main",
            i = "bor", 
            j = "all_adj", 
            x_axis = "plot_est", 
            y_axis = "log10_p", 
            subset = "overall", 
            k = .5, 
            simple = "False", 
            nudge_y = 1,
            repel = TRUE, 
            push = 1,
            pull = 1, 
            seed = 622
    ) + exhaustive_theme_left + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.04)) 
      + scale_x_continuous(breaks = c(1,2), limits = c(.4,2.5))     
      + geom_vline(xintercept = 1, linetype="dashed", color = "black", size=.1)
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1) 
    )

    m1 <- (
        super_plot(
            ready, 
            labels %>% filter(clean_label %in% c("Prior Systemic Therapy", "TMB", "T-cell")), 
            type = "main",
            i = "pfs", 
            j = "all_adj", 
            x_axis = "plot_est", 
            y_axis = "log10_p", 
            subset = "overall", 
            k = .5, 
            simple = "False", 
            nudge_y = 1,
            repel = TRUE, 
            push = 1,
            pull = 1, 
            seed = 622
    ) + exhaustive_theme_right + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.04))
      + scale_x_continuous(breaks = c(.5,1,1.5), limits = c(.67,1.6))
      + geom_vline(xintercept = 1, linetype="dashed", color = "black", size=.1)  
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1)  
    )

    m2 <- (
        super_plot(
            ready,
            labels %>% filter(clean_label %in% c("Prior Systemic Therapy", "TMB", "T-cell")), 
            type = "main",
            i = "os", 
            j = "all_adj", 
            x_axis = "plot_est", 
            y_axis = "log10_p", 
            subset = "overall", 
            k = .5, 
            simple = "False", 
            nudge_y = 1,
            repel = TRUE, 
            push = 1,
            pull = 1, 
            seed = 622
    ) + exhaustive_theme_right + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.04)) 
      + scale_x_continuous(breaks = c(.5,1,1.5), limits = c(.65,1.9))   
      + geom_vline(xintercept = 1, linetype="dashed", color = "black", size=.1)  
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1) 
    )

    m3 <- (
        super_plot(
            ready %>% filter(Type != "Somatic: Mutation", cor_tcell < .5, cor_pretreat < .5), 
            labels %>% filter(clean_label %in% c("TGFB", "Proliferation")), 
            type = "main",
            i = "os", 
            j = "tmb_tcell_clin_adj", 
            x_axis = "plot_est", 
            y_axis = "log10_p", 
            subset = "overall", 
            k = .5, 
            simple = "False", 
            nudge_y = 1.2,
            repel = TRUE
    ) + exhaustive_theme_right_legend + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + ggtitle("OS Residuals vs Features") 
      + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.1))   
      + scale_x_continuous(breaks = c(.75,1,1.25), limits = c(.6,1.4))  
      + geom_vline(xintercept = 1, linetype="dashed", color = "black", size=.1)  
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.2) 
      + guides( color = "none", alpha = "none", size = "none", fill = guide_legend(override.aes = list(size=3)))  
      + annotate(geom="text", x=1.2, y=4.1, label="BY threshold",color="black")  
    ) 
    leg <- get_legend(m3)

    if( type == "non_latent"){
        m0 = m0 + ggtitle("BOR vs non-latent features") + scale_color_manual( values = c("white"))
        m1 = m1 + ggtitle("PFS vs non-latent features") + scale_color_manual( values = c("white"))
        m2 = m2 + ggtitle("OS vs non-latent features") + scale_color_manual( values = c("white"))
        m3 = m3 + ggtitle("OS Residuals vs non-latent features") + scale_color_manual( values = c("white"))
    } else if (type == "latent"){
        m0 = m0 + ggtitle("BOR vs latent features")
        m1 = m1 + ggtitle("PFS vs latent features")
        m2 = m2 + ggtitle("OS vs latent features")
        m3 = m3 + ggtitle("OS Residuals vs latent features")
    } else {
        m0 = m0 + ggtitle("BOR vs features")
        m1 = m1 + ggtitle("PFS vs features")
        m2 = m2 + ggtitle("OS vs features")
        m3 = m3 + ggtitle("OS Residuals vs features")
    }
    basic <- arrangeGrob(m0,m1,m2, m3 + guides(fill = "none"), leg, layout_matrix = rbind( c( rep(1,10), rep(2,9), rep(3,9),rep(4,9),rep(5,4))))  
    get_dressed(basic, title, vjust = 0, mar = 0, size = 18, hjust = .45)
}   

base <- annote(latent_plot(type = "base", k = 1, "All features"), "A")
latent <- annote(latent_plot(type = "latent", k = .7, "Features with > .7 correlation to latent factors"), "B")
non_latent <- annote(latent_plot(type = "non_latent", k = .3, "Features with < .3 correlation to latent factors"), "C")
together <- arrangeGrob( base, latent, non_latent, ncol = 1)
go <- get_dressed(together, "HMF Exhaustive Analysis - Latent factors", size = 22, vjust = 2)

options(repr.plot.width = 17, repr.plot.height = 17, repr.plot.res = 100)
plot(go)
ggsave( paste0(o_dir, "sn_exhaustive_latent.png"), go, width = 17, height = 17, dpi = 300)
