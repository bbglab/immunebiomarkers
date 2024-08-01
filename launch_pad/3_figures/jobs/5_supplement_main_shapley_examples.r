wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

library(tidyverse)
library(gridExtra)
library(ggpubr)

annote <- function( i, lab, size = 25 ) annotate_figure( i, fig.lab = lab, fig.lab.size = size, fig.lab.face = "bold")

go <- readRDS(paste0(FIG_DIR, "figure5_g.Rds"))
main_p <- annote(go$main + theme(plot.title = element_text(size = 35), axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.margin = unit(c(.5,.5,1.5,.5), "cm")), "A")
exs <- go$exs

labs <- seq(4); 
j <- 1; 
clean_exs <- list()
for( i in names(exs) ){
    clean_exs[[i]] <- annote(exs[[i]] + theme(plot.margin = unit(c(7,2,.5,2), "cm")), lab = as.character(labs[j])) 
    j <- j+1
}

top <- arrangeGrob(main_p, layout_matrix = rbind(c(2,1,1,1,3)))

options(repr.plot.width=24, repr.plot.height=36, repr.plot.res = 200)
top <- arrangeGrob(main_p, layout_matrix = rbind(c(2,1,1,1,3)))
bottom <- arrangeGrob( clean_exs[[1]], clean_exs[[2]], clean_exs[[3]], 
            clean_exs[[4]], 
            layout_matrix = rbind(seq(2), seq(2)+2))
go <- as_ggplot(arrangeGrob(top, bottom, layout_matrix = cbind(c(1,1,2,2,2))))

options(repr.plot.width=20, repr.plot.height=25, repr.plot.res = 200)
go
ggsave( paste0(FIG_FINAL_DIR, "5_supplement_main_shapley_examples.png"), width = 20, height = 25)
