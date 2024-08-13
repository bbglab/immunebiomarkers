wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

library(tidyverse)
library(gridExtra)
library(ggpubr)

ef <- readRDS(paste0(FIG_DIR, "figure5_ef.Rds"))

get_dressed <- function(go, hjust = .5, size = 25, vjust = 4, m = 1){
    go + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(rep(m,4), "cm"))
}
annote <- function( i, lab ) annotate_figure( i, fig.lab = lab, fig.lab.size = 35, fig.lab.face = "bold")
left_title <- theme(plot.title = element_text(hjust = 0))

fig_e <- annote(get_dressed( ef$e_main,  "Overall survival by response strata", m = .5, size = 30 ), "A")
fig_f <- annote(get_dressed( as_ggplot(ef$low),  "Overall survival by response strata", m = .5, size = 30 ), "B")
lower_panel <- as_ggplot(arrangeGrob( fig_e, fig_f, fig_g, ncol = 3))
lower <- get_dressed(lower_panel, m = 1, size = 40) + ggtitle("HMF response vs survival and shapley example")

figure5_supp <- as_ggplot(arrangeGrob( lower, ncol = 1))

options(repr.plot.width = 30, repr.plot.height= 11, resolution = 200)
figure5_supp
ggsave(file = paste0(FIG_FINAL_DIR, "figure5_supp.png"), width = 30, height = 11)
