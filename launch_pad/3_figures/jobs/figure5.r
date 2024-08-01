wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

library(tidyverse)
library(gridExtra)
library(ggpubr)

b <- readRDS(paste0(FIG_DIR, "figure5_b.Rds"))
c <- readRDS(paste0(FIG_DIR, "figure5_c.Rds"))
d <- readRDS(paste0(FIG_DIR, "figure5_d.Rds"))

get_dressed <- function(go, hjust = .5, size = 25, vjust = 4, m = 1){
    go + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(rep(m,4), "cm"))
}
annote <- function( i, lab ) annotate_figure( i, fig.lab = lab, fig.lab.size = 35, fig.lab.face = "bold")
left_title <- theme(plot.title = element_text(hjust = 0))

fig_b <- annote(get_dressed(as_ggplot(b$abc) + ggtitle("Patients Stratified by Response Groups"), size = 35), "B")
fig_c <- annote(get_dressed( as_ggplot(c$response) , size = 35) + ggtitle("Surivival by Response Groups"), "C")
fig_d <- annote(get_dressed( d , size = 35) + ggtitle("Patients Stratified by TMB High or Low (TMB per MB > 10)"), "D")

layout <- cbind(c(rep(1,5), rep(2,4)), c(rep(1,5), rep(2,4)), c(rep(3,5), rep(4,4)), c(rep(3,5), rep(4,4)), c(rep(3,5), rep(4,4)))

go <- as_ggplot(arrangeGrob(ggplot() + geom_blank(), fig_c, fig_b, fig_d, layout_matrix = layout))

go_go <- (get_dressed(go, m = 1, size = 40, vjust = 2) )

paste0(FIG_FINAL_DIR, "figure5_new.png")

options(repr.plot.width = 30, repr.plot.height= 20, resolution = 200)
go_go
ggsave(file = paste0(FIG_FINAL_DIR, "figure5_new.png"), width = 30, height = 20)
ggsave(file = paste0(FIG_FINAL_DIR, "figure5_new.pdf"), width = 30, height = 20)
