wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)
library(gridExtra)
library(ggpubr)

pred_plts <- readRDS(paste0(FIG_DIR, "pred_study_plots.Rds"))

get_dressed <- function(gg, title = "", hjust = .5, size = 27, vjust = 1, mar = c(0,0,.2,.2)){
    gg + ggtitle(title) + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(mar, "cm"))
}
annote <- function( i, lab, size = 25 ) annotate_figure( i, fig.lab = lab, fig.lab.size = size, fig.lab.face = "bold")

a <- get_dressed( annote(arrangeGrob(pred_plts$combined$figure_sn2_1), "A"), title = "Modelling Strategies: Comparison of performance", mar = rep(1,4))
b <- get_dressed( annote(pred_plts$combined$figure_s6a, "B"), title = "Hybrid Models: Evaluation over different features", mar = rep(1,4))
c <- get_dressed( annote(pred_plts$combined$figure_sn7, " "), title = "Evaluation of different parameter settings for XGBoost models", mar = rep(1,4))
go <- as_ggplot(arrangeGrob(a, b))

options(repr.plot.width = 20, repr.plot.height= 23, resolution = 200)
go
ggsave( paste0(FIG_FINAL_DIR, "0a_supplement_main_sim_modelling.png"), go, width = 20, height = 23)

options(repr.plot.width = 20, repr.plot.height= 11.5, resolution = 200)
c
ggsave( paste0(FIG_FINAL_DIR, "0a_supplement_main_sim_modelling_7b.png"), c, width = 20, height = 11.5)
