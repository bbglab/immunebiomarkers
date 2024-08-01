library(ggplot2)
library(gridExtra)

### Base theme for project ### 
base_theme <- theme_classic(base_size = 16)

### Exhaustive Plots 
exhaustive_theme_left <- (
    base_theme + 
    theme(
            plot.margin = unit(c(.4,0,1.6,0), "cm"), 
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.title.x=element_text(size=16),
            axis.text.x = element_text(size = 16),
            legend.position="none"
     ) 
)
exhaustive_theme_right <- (
    exhaustive_theme_left + 
    theme(
        axis.title.y=element_blank(), 
        axis.text.y=element_blank()
    )
)
exhaustive_theme_right_legend <- exhaustive_theme_right + theme(legend.position = "right")
combine_exhaustive0 <- function( p0, p1 ) {arrangeGrob(p0,p1,layout_matrix = rbind( c(rep(1,10), rep(2,9))))}
combine_exhaustive1 <- function( p0, p1, p2) {arrangeGrob(p0,p1,p2, layout_matrix = rbind( c(rep(1,10), rep(2,9), rep(3,9))))}
combine_exhaustive2 <- function( p0, p1, p2, p3) {arrangeGrob(p0,p1,p2,p3, layout_matrix = rbind( c(rep(1,10), rep(2,9), rep(3,9), rep(4,9))))}


### Making grobs back into ggplots - get dressed! 
get_dressed <- function(grob, title, hjust = .5, size = 25, vjust = 1, mar = .5){
    as_ggplot(grob) + ggtitle(title) + theme(plot.title = element_text(hjust = hjust, size = size, vjust = vjust), plot.margin = unit(c(mar,mar,mar,mar), "cm"))
}
