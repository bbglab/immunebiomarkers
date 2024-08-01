wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

library(tidyverse)
library(gridExtra)
library(ggpubr)

a <- readRDS(paste0(FIG_DIR, "6_sm_stratify_external_a.Rds"))
b <- readRDS(paste0(FIG_DIR, "6_sm_stratify_external_b.Rds"))

get_dressed <- function(grob, hjust = .5, size = 25, vjust = 4, m = 1){
    as_ggplot(grob) + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(rep(m,4), "cm"))
}
annote <- function( i, lab ) annotate_figure( i, fig.lab = lab, fig.lab.size = 35, fig.lab.face = "bold")
left_title <- theme(plot.title = element_text(hjust = 0))

fig_a <- annote(a, "A") + theme(plot.margin = unit(c(1,1,1,1), "cm"))
fig_b <- annote(b, "B") + theme(plot.margin = unit(c(1,7,1,1), "cm"))

supp_strata <- arrangeGrob(fig_a, fig_b, layout_matrix = cbind(c(1,1,1,2,2,2,2)))

paste0(FIG_FINAL_DIR, "6_supplement_main_stratify_external.pdf")

options(repr.plot.width = 16, repr.plot.height= 20, resolution = 200)
get_dressed( supp_strata , size = 31) + ggtitle("Validation cohorts with complete data on five factors (n = 179)")
ggsave(file = paste0(FIG_FINAL_DIR, "6_supplement_main_stratify_external.png"), width = 16, height = 20)
ggsave(file = paste0(FIG_FINAL_DIR, "6_supplement_main_stratify_external.pdf"), width = 16, height = 20)
