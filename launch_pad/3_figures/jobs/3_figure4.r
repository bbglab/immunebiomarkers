wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)
library(ggpubr)
library(gridExtra)

base_theme <- theme_classic(base_size = 21) 

theme <- (
    base_theme + 
    theme(
        axis.title.y=element_blank(), 
        plot.title = element_text(hjust = 0.5, size =25),#, size = 21),
        legend.position="none"
    )
)

forest_plots <- readRDS(paste0(FIG_DIR, "forest-plots.Rds"))

FIG_DIR

names(forest_plots)

silhouette <- readRDS(paste0(FIG_DIR, "silhouette_values_sqrt_new_clustering_data.Rds"))
sil_values<-(ggplot(silhouette$silhouette_results, aes(x=factor(gene_name,levels=silhouette$gene_levels),y=sil_width,fill=cluster_annotated)) 
    + geom_bar(stat = "identity") 
    + facet_grid(cols = vars(dataset), rows = vars(factor(cluster_annotated, levels = c("tcell", "proliferation","tgfb"))),scales="free_y")
    + scale_fill_manual(values = c("#8DD3C7", "#FB8072", "#BEBADA"))
    + theme_bw()
    + coord_flip()
    + xlab("")
    + ylab("silhuette score")
    + theme
    + theme(axis.text.y=element_text(size=7), legend.position = "none")
)

annote <- function( i, lab ) annotate_figure( i, fig.lab = lab, fig.lab.size = 20, fig.lab.face = "bold")

plot_a <- annote(forest_plots$multi  + theme(plot.margin = unit(c(0,0,2,0), "cm")), "A")

plot_b <- annote( sil_values + ggtitle("Expression Silhouette Values") + theme(plot.margin = unit(c(0,1,0,1), "cm")), "B")

options(repr.plot.width = 20, repr.plot.height= 24, resolution = 200)

together <- arrangeGrob( plot_a, plot_b)

plot( as_ggplot(together))

ggsave( paste0(FIG_DIR, "fig4_new.png"), width = 20, height = 24, dpi = 300)
