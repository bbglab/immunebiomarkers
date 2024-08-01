wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
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
threshold <- -log10(ingredients$by_05_fdr)[1]

o_dir <- paste0(FIG_DIR ,"supplement_note/")
annote <- function( i, lab, size = 16 ) annotate_figure( i, fig.lab = lab, fig.lab.size = size, fig.lab.face = "bold")

make_cor_plots <- function( model = "bor", title){

    s1 <- (
        super_plot(
            ingredients, # %>% filter( abs(cor_tmb) < k, abs(cor_t_cell) < k, abs(cor_tgfb) < k, abs(cor_prolif) < k, abs(cor_pre_treat) < k), 
            labels %>% filter(! clean_label %in% c("TMB clonal", "TMB subclonal", "BRCA2", "MSH2")), 
            type = "cor",
            i = model, 
            j = "all", 
            x_axis = "cor_tmb", 
            y_axis = "log10_p", 
            subset = "somatic", 
            k = .5, 
            simple = "False", 
            repel = TRUE, 
            nudge_y = .8, 
            nudge_x = 0,
            push = 10
    ) + exhaustive_theme_left + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + ggtitle(paste0( toupper(model)," vs Somatic Features"))
      + scale_x_continuous(n.breaks = 3, limits = c(-.7,1.15))  
      + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=.1)  
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1)  
    )
    
    c1 <- (
        super_plot(
            ingredients %>% filter( Type %in% c("Clinical", "HLA"), ! grepl("composite", feature)), 
            labels %>% filter( Type %in% c("Clinical", "HLA"), ! grepl("composite", feature), !clean_label %in% c("pre_treated","HLA LOH", "Prior Radiotherapy"), clean_label != "Prior Chemotherapy"), 
            type = "cor",
            i = model, 
            j = "all_adj", 
            x_axis = "cor_pretreat", 
            y_axis = "log10_p", 
            subset = "overall", 
            k = .5, 
            simple = "False",
            repel = TRUE, 
            nudge_y = .8, 
            nudge_x = 0,
            push = 10
    ) + exhaustive_theme_right + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
      + ggtitle(paste0( toupper(model), " vs Clinical + HLA Features"))
      + scale_x_continuous(breaks = c(0,1), limits = c(-.4,1.05))  
      + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=.1)  
      + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1)     
      + guides( fill = "none", alpha = "none", size = "none", shape = "none", color = guide_legend(override.aes = list(size=3)))   
    )
    if(model == "bor"){
        c1 <- c1 + annotate("text", x=-.2, y=4, label= "BY Threshold")
    }


    leg <- get_legend(c1)

    t1 <- (
        super_plot(
            ingredients, # %>% filter( abs(cor_tmb) < k, abs(cor_t_cell) < k, abs(cor_tgfb) < k, abs(cor_prolif) < k, abs(cor_pre_treat) < k), 
            labels %>% filter(!clean_label %in% c("TGFB", "Proliferation", "CXCL9", "Macrophages M0","Macrophages M2", "CD274")), 
            type = "cor",
            i = model, 
            j = "all_adj", 
            x_axis = "cor_tcell", 
            y_axis = "log10_p", 
            subset = "rna", 
            k = .5, 
            simple = "FALSE", 
            repel = TRUE, 
            nudge_y = .8, 
            nudge_x = 0,
            push = 10
        ) + exhaustive_theme_right + theme(plot.margin = unit(c(.4,0,1,0), "cm"))
          + ggtitle(paste0( toupper(model), " vs RNA Features"))
          + scale_x_continuous(n.breaks = 3, limits = c(-.8,1.15))
          + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=.1)  
          + geom_vline(xintercept = .5, linetype="dashed", color = "black", size=.1)  
          + geom_hline(yintercept = threshold, linetype="dashed", color = "black", size=.1)          
    )
    
    ### Add Y-scales
    if(model == "bor"){
        s1 <- s1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9))
        t1 <- t1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9))
        c1 <- c1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9)) + guides(color = "none")
    } else if (model == "os"){
        s1 <- s1 + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.1))
        t1 <- t1 + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.1)) 
        c1 <- c1 + scale_y_continuous(breaks = c(3,6,9,12), lim = c(0,12.1)) + guides(color = "none")    
    } else {
        s1 <- s1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9.8) )
        t1 <- t1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9.8))
        c1 <- c1 + scale_y_continuous(breaks = c(3,6,9), lim = c(0,9.8)) + guides(color = "none")    
    }
    
    basic <- combine_exhaustive1(s1,t1,c1)
    get_dressed(basic, title, size = 18, mar = 0, vjust = 0)
}

bor <- annote(make_cor_plots("bor", "BOR vs all features"), "A")
pfs <- annote(make_cor_plots("pfs", "PFS vs all features"), "B")
os <- annote(make_cor_plots("os", "OS vs all features"), "C")
together <- arrangeGrob(bor, pfs, os)
ready <- get_dressed(together, title = "HMF Exhaustive Analyses - Correlations with TMB, T-cell, Prior systemic therapy", size = 18, vjust = 3)

options(repr.plot.width = 12, repr.plot.height = 15, repr.plot.res = 400)
plot(ready)
ggsave( paste0(o_dir, "sn_exhaustive_correlations.png"), ready, width = 12, height = 15, dpi = 300)
