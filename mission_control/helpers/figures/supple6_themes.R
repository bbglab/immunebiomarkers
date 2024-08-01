library(ggplot2)
library(RColorBrewer)

### themes 
base_theme <- (
    theme_classic(base_size = 25)  + 
    theme( legend.position = "none", 
           plot.margin = unit(c(.3,0,.2,1.5), "cm"), 
           plot.title = element_text(hjust = 0.5, size = 30),
           axis.title.y= element_blank(),
           axis.ticks.y= element_blank(), 
           axis.text.y = element_blank(),
           axis.line.y = element_blank(),
           axis.text.x = element_text(size = 25)
         )
)

### bar plots 
theme_bar_main <- base_theme
theme_bar_tissue <- base_theme + theme(axis.title.x = element_text(size = 13), axis.text.x = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 25))

### compute plots 
theme_compute_main <- base_theme + theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,1.2,1), "cm"))
theme_compute_tissue <- base_theme + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 25), plot.margin = unit(c(0,0,1.2,1), "cm"))

### dist plots 
theme_dist <- (
    base_theme + 
    theme(axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = unit(c(1,0,0,1.5), "cm")
         )
)

layout <- rbind(c(1,1), c(1,1), c(1,1), c(2,2), c(2,2))
layout2 <- rbind(c(1,1), c(1,1), c(1,1), c(1,1), c(2,2), c(3,3), c(3,3))

### colors 
my_palette<-brewer.pal(6,"RdYlGn")[c(1,3,6)]
