wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/pred_study.R"))
source(paste0(wd,"/mission_control/helpers/figures/general.R"))

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(survminer)
library(ggpubr)

evals <- readRDS(paste0(TMP_DIR,"xg-eval-results.Rds"))

data_long <- (
    gather(evals, type, eval_metric, simple:all, factor_key=TRUE) 
       %>% drop_na(eval_metric)
       %>% mutate( tissue = ifelse(tissue == "all", "Pan-Cancer", str_to_title(tissue) ) )
)
data_long$tissue <- factor(
    data_long$tissue,
    levels = c("Pan-Cancer", "Skin", "Lung", "Bladder", "Other")
)
data_long$fit <- paste0(data_long$type,"-", data_long$features)
data_long <- rbind(data_long, data_long %>% filter(fit == "warm-five_latent") %>% mutate(fit = "hybrid"))

results <- (
    data_long
         %>% filter(complete) 
         %>% group_by(model, fit, tissue, type, purity) 
         %>% summarise(mn = mean(eval_metric, na.rm = TRUE))
    )

#table(results$fit)
#results %>% filter(fit == "warm-five_latent_complex")

name_map <- list(
    "warm-tmb" = 'TMB\n',
    "warm-base" = 'Base: TMB + PDL1\n',
    "warm-rna_only" = 'RNA only: T-cell, TGFB, Proliferation\n',
    "warm-no_tmb" = 'No TMB: T-cell + TGFB + Proliferation + Pretreat\n', 
    "warm-five_latent" = 'Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n', 
    "hybrid" = 'Hybrid Models\n', 
    "warm-full_mod" = 'Full Model: Five Latent + Age + Distal Biopsy + HLA LOH \n+ CPI Mechanism + WGD + SVB',
    "simple-five_latent" = 'Simple Tissue Models\n',
    "all-five_latent" = 'Pan-Cancer Model\n',
    "warm-latent1" = 'Five Latent: Alernative Measures 1\n',
    "warm-latent2" = 'Five Latent: Alernative Measures 2\n',
    "warm-five_latent_interaction" = 'Five Latent: Allow Interaction\n',
    "warm-five_latent_complex" = 'Five Latent: Allow Interaction, Full Grid Search\n'
)

results <- results %>% filter( fit %in% names(name_map))
data_long <- data_long %>% filter( fit %in% names(name_map))

results$fit2 <- factor(unlist(lapply(results$fit, function(i) name_map[[i]])), levels = unname(unlist(name_map)))
data_long$fit2 <- factor(unlist(lapply(data_long$fit, function(i) name_map[[i]])), levels = unname(unlist(name_map)))

studies <- list( 
    "figure_5b" = c('Base: TMB + PDL1\n', 
                    'Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n'),
    "figure_s6a" = c('TMB\n', 
                     'Base: TMB + PDL1\n', 
                     'RNA only: T-cell, TGFB, Proliferation\n', 
                     'No TMB: T-cell + TGFB + Proliferation + Pretreat\n',
                     'Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n', 
                     'Full Model: Five Latent + Age + Distal Biopsy + HLA LOH \n+ CPI Mechanism + WGD + SVB'),
    "figure_s6b" = c('Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n', 
                     'Five Latent: Alernative Measures 1\n', 
                     'Five Latent: Alernative Measures 2\n'),
    "figure_sn2_1" = c('Hybrid Models\n', 
                       'Simple Tissue Models\n',
                       'Pan-Cancer Model\n'),
    "figure_sn7" = c('Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n',
                     'Five Latent: Allow Interaction\n',
                     'Five Latent: Allow Interaction, Full Grid Search\n'
                    )
)

color_maps <- list(
    'Base: TMB + PDL1\n' = '#A6CEE3',
    'Five Latent: TMB + T-cell + TGFB + Proliferation + Pretreat\n' = '#1F78B4', 
    'TMB\n' = '#B2DF8A',
    'RNA only: T-cell, TGFB, Proliferation\n' = 'grey',
    'No TMB: T-cell + TGFB + Proliferation + Pretreat\n' = 'orange',
    'Full Model: Five Latent + Age + Distal Biopsy + HLA LOH \n+ CPI Mechanism + WGD + SVB' = '#FB9A99',
    'Hybrid Models\n' = '#1F78B4',
    'Simple Tissue Models\n' = '#E31A1C',
    'Pan-Cancer Model\n' = '#FDBF6F',
    'Five Latent: Alernative Measures 1\n' = '#6A3D9A',
    'Five Latent: Alernative Measures 2\n' = '#CAB2D6',
    'Five Latent: Allow Interaction\n' = 'dark grey',
    'Five Latent: Allow Interaction, Full Grid Search\n' = 'white'    
)

get_cols <- function(fits) unlist(lapply( fits, function(i) color_maps[i]))

title_maps <- list(
    "t1" = list("lr" = "Test AUC", "pfs" = "Test Concordance Index", "os" = "Test Concordance Index"),
    "t2" = list("lr" = "HMF Prediction Study - Response", "pfs" = "HMF Prediction Study - Progression Free Survival", "os" = "HMF Prediction Study - Overall Survival")
)

j <- "true"

ready_data <- list()
for( i in names(studies)){
    ready_data$raw[[i]] <- data_long %>% filter(fit2 %in% studies[[i]], purity == j)
    ready_data$bars[[i]] <- results %>% filter(fit2 %in% studies[[i]], purity == j)
}

#ready_data$bars[[i]]

plts <- list()
for ( i in names(studies) ){
    cols <- get_cols(studies[[i]])
    for( j in c("lr", "pfs", "os")){ 
        plts[[i]][[j]] <- get_dressed( 
            pred_plots( ready_data$raw[[i]], ready_data$bars[[i]], mod = j, cols ), title = title_maps$t2[[j]]
        )
        tmp <- (
            barplots(ready_data$bars[[i]], ylab = "", theme = theme_bars_main, size = 4, cols) + 
            theme(legend.position = "right", 
                  legend.text = element_text(size=20), 
                  legend.title = element_blank(),
                  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) 
        )
        if( i == "figure_5b"){
            plts[[i]][['legend']] <- get_legend(tmp + theme(legend.direction="horizontal", legend.key.height = unit(1, 'cm'), legend.key.width = unit(1.5, 'cm')))
        } else {
            plts[[i]][['legend']] <- get_legend(tmp)
        }
    }
}

o_simple <- function( plts, i ) {
    arrangeGrob(plts[[i]]$lr, plts[[i]]$os, plts[[i]]$legend, layout_matrix = rbind(c(1,2), c(1,2), c(1,2), c(1,2), c(1,2), c(1,2), c(3,3)))
}
o_full <- function( plts, i ) {
    arrangeGrob(plts[[i]]$lr, plts[[i]]$legend, plts[[i]]$pfs, plts[[i]]$os,  ncol = 2)
}

prediction_plots <- list(
    "figure_5b" = as_ggplot(o_simple( plts, 'figure_5b')),
    "figure_s6a" = as_ggplot(o_full( plts, 'figure_s6a')),
    "figure_s6b" = as_ggplot(o_full( plts, 'figure_s6b')),
    "figure_sn2_1" = as_ggplot(o_full( plts, 'figure_sn2_1')),
    "figure_sn7" = as_ggplot(o_full( plts, 'figure_sn7'))
)

saveRDS( list( "individual" = plts, "combined" = prediction_plots), file = paste0(FIG_DIR, "pred_study_plots.Rds"))
