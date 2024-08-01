wd <- dirname(dirname(getwd()))
setwd(wd)

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure2_theme.R"))

library(tidyverse)

ready <- readRDS(paste0(I_DIR, "cpi_go.Rds"))

fig_2c <- (
  ggplot(ready , aes(x = isofox_gene_set_t_cell_effector, y = tcell)) + 
    geom_point(color = "white", fill = "grey20", shape = 21, size = 3) + 
    ggpubr::stat_cor(aes(label = after_stat(r.label)), size = 6) + 
    ylab("T-cell effector gene set") +
    xlab("Mean expression Cluster 3 genes") + 
    ggtitle("T-cell effector gene set vs Cluster 3") +
    base_theme
)

saveRDS(fig_2c, file = paste0(FIG_DIR, "fig_2c.Rds"))
