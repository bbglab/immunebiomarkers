library(ggplot2)
base_theme <- (
    theme_classic(base_size = 18) + 
    theme(
            plot.margin = unit(c(.4,0,1.6,0), "cm"), 
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.title.x=element_text(size=14),
            axis.text.x = element_text(size = 14),
            legend.position="none",
            legend.background = element_blank(), 
            legend.box.background = element_blank(), 
            legend.title=element_blank()
        )
)
