library(ggplot2)
base_theme <- theme_classic(base_size = 20) 
theme1 <- (
    base_theme 
        + theme(plot.margin = unit(c(0,0,0,0), "cm"), 
                axis.text = element_text(size = 16), 
                axis.title = element_text(size = 18), 
                plot.title = element_text(hjust = .5, size = 20), 
                legend.position = "none")
)        
theme2 <- theme1 + theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 12))
