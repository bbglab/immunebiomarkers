library(ggpubr)
library(ggplot2)

base_theme <- theme_classic(base_size = 21) 

color_map <- list(
    "Prior Systemic Therapy (Composite)" = '#FFFFB3',
    "TGFB" = "#BEBADA",
    "Proliferation" = "#8DD3C7",
    "T-cell" = '#FB8072',
    "TMB" = '#80B1D3'
)
get_dressed <- function(grob, title, hjust = .5, size = 25, vjust = 1, mar = .5){
    as_ggplot(grob) + ggtitle(title) + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(rep(mar,4), "cm"))
}

