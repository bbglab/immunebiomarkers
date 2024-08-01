wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure2_theme.R"))

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds")) %>% mutate(Direction = factor(ifelse(Direction == "Worse", "Worse Outcome", "Better Outcome"), levels = c("Worse Outcome", "Better Outcome")))
labels <- readRDS(paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
labels$size <- labels$size
threshold <- -log10(ingredients$by_05_fdr)[1]

size_map <- list( "TRUE" = 4, "FALSE" = 2)
alpha_map <- list( "TRUE" = 1, "FALSE" = .4)
color_map <- list( "TRUE" = "black", "FALSE" = "white")
shape_map <- list(  "Better Outcome" = 24, "Worse Outcome" = 25)
fill_map <- list("Clinical"="#9E7788",
                 "CNV/SVs"="#B3DE69",
                 "HLA"="#FFFF99",
                 "RNA: Remaining"= "#E0BA92", #"#FC913A",
                 "Somatic"="#80B1D3",
                 "RNA: T-cell" = "#FB8072",
                 "RNA: TGFB" = "#BEBADA", 
                 "RNA: Proliferation" = "#8DD3C7")

aes_tmb <- aes(x = cor_tmb, y = log10_p, shape = Direction, fill = discovery_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_tcell <- aes(x = cor_tcell, y = log10_p, shape = Direction, fill = discovery_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_pretreat <- aes(x = cor_pretreat, y = log10_p, shape = Direction,fill = discovery_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_tgfb <- aes(x = cor_tgfb, y = log10_p, shape = Direction,fill = discovery_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main
aes_prolif <- aes(x = cor_prolif, y = log10_p, shape = Direction, fill = discovery_group, size = log10_p > threshold, alpha = log10_p > threshold, color = log10_p > threshold ) #+ guides_main

guides_main <- guides(alpha = "none", size = "none", color = "none", shape = guide_legend(override.aes = list(size=7, ncol = 2)), fill = guide_legend(override.aes = list(shape=21, size=7, ncol=2)))

all <- ingredients %>% filter(dataset == "all")

max_y <- max(all %>% filter(covariates == "age_biopsy_purity_tissue") %>% pull(log10_p))

plots_go <- function( df = all, type = "tmb",  xlab = "X", title = "X"){
    if( type == "tmb"){ gg <- ggplot(df, aes_tmb) } 
    else if (type == "tcell"){gg <- ggplot(df, aes_tcell) } 
    else if (type == "pre"){ gg <- ggplot(df, aes_pretreat)}
    else if (type == "tgfb"){ gg <- ggplot(df, aes_tgfb)}
    else if (type == "prolif"){ gg <- ggplot(df, aes_prolif)}
    (gg + 
        geom_point(stroke = 0.3) + 
        scale_y_continuous(breaks = c(seq(0,max_y,2)), limits = c(0,max_y)) + 
        scale_x_continuous(n.breaks = 3)  +
        scale_size_manual(values = unlist(size_map))  +
        scale_alpha_manual(values = unlist(alpha_map)) + 
        scale_color_manual(values = unlist(color_map) ) + 
        scale_shape_manual(values = unlist(shape_map)) + 
        scale_fill_manual(values = unlist(fill_map), limits = force) + 
        geom_hline( yintercept = threshold, linetype="dashed", color = "black" ) + 
        xlab( xlab ) + 
        ylab("-Log10 p-value") +  
        ggtitle( title) +
        base_theme
    )
}

remove_y  <- theme(axis.text.y=element_blank(), axis.title.y=element_blank())

make_layer <- function( df = all, type = "tmb", overall_title = "X", aes_go = aes_tmb, xlab = "Correlation TMB" ){
    p1 <- plots_go(df %>% filter(model == "bor", covariates == "age_biopsy_purity_tissue"), type = type, xlab = xlab, title = paste0("BOR vs All features"))
    p2 <- plots_go(df %>% filter(model == "pfs", covariates == "age_biopsy_purity_tissue"), type = type, xlab = xlab, title = paste0("PFS vs All features")) + remove_y
    p3 <- plots_go(df %>% filter(model == "os", covariates == "age_biopsy_purity_tissue"), type = type, xlab = xlab, title = paste0("OS vs All features")) + remove_y
    p4 <- plots_go(df %>% filter(model == "os", covariates == "residuals"), type = type, xlab = xlab, title = paste0("OS Residuals vs All features")) + remove_y
    legend <- get_legend(p3 + theme(legend.position = "right") + guides_main)
    
    if( type %in% c("tmb", "tcell","pre")){
        as_ggplot(arrangeGrob(p1,p2,p3, ncol = 3)) + theme(plot.title = element_text(size = 18, vjust = 0, hjust = .4)) + ggtitle(overall_title) + theme(plot.margin = margin(.1,.1,.1,.1, "cm"))
    } else {
        as_ggplot(arrangeGrob(p1,p3,p4, ncol = 3)) + theme(plot.title = element_text(size = 18, vjust = 0, hjust = .4)) + ggtitle(overall_title) + theme(plot.margin = margin(.1,.1,.1,.1, "cm"))      
    }
}

options(repr.plot.width = 14, repr.plot.height = 4.5)

annote <- function( i, lab ) annotate_figure( i, fig.lab = lab, fig.lab.size = 18, fig.lab.face = "bold")

tmb <- annote(make_layer(df = all, overall_title = "Signficance vs Correlation with TMB Cluster", type = "tmb", xlab = "Correlation TMB"), "A")
pre <- annote(make_layer(df = all, overall_title = "Signficance vs Correlation with Pretreatment Cluster", type = "pre", xlab = "Correlation Pretreatment"), "B")
tcell <- annote(make_layer(df = all, overall_title = "Signficance vs Correlation with T-cell Cluster", type = "tcell", xlab = "Correlation T-cell"), "C")
tgfb <- annote(make_layer(df = all, overall_title = "Signficance vs Correlation with TGFB Cluster", type = "tgfb", xlab = "Correlation TGFB"), "D")
prolif <- annote(make_layer(df = all, overall_title = "Signficance vs Correlation with Proliferation Cluster", type = "prolif", xlab = "Correlation Proliferation"), "E")

legend <- get_legend(plots_go(all %>% filter(model == "os", covariates == "age_biopsy_purity_tissue"), title = "OS vs All features") + remove_y + theme(legend.position = "right") + guides_main+ theme(plot.margin = margin(.1,.1,.1,.1, "cm")))

options(repr.plot.width = 19, repr.plot.height = 14)

lets_go <- as_ggplot(arrangeGrob( tmb, legend, pre, tcell,  tgfb, prolif, ncol = 2)) + ggtitle("HMF CPI: Exhaustive Analysis by Correlation to 5 Latent Factors")

lets_go + theme(plot.margin = margin(.5,.5,.5,.5, "cm"), plot.title = element_text(vjust = 2, hjust = .4, size = 18))

ggsave(file = paste0(FIG_FINAL_DIR, "7c_supplement_note_exhaustive_correlation_go.png"), width =19, height = 14)

ggsave(file = paste0(FIG_FINAL_DIR, "7c_supplement_note_exhaustive_correlation_go.pdf"), width =19, height = 14)
