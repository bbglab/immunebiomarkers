library(ggplot2)
library(gridExtra)

### THEMES
base_theme <- theme_classic(base_size = 21)

theme_violins_main <- (
    base_theme + 
    theme( axis.ticks.x = element_blank(),axis.text.x=element_blank(),
           axis.title.x=element_blank(),legend.position="none",
           plot.margin = unit(c(0,0,0,0), "cm")) 
)
theme_violins_tissue <- (
    theme_violins_main + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank())
)

### Main and tissues 
theme_bars_main <- (
    theme_violins_main + 
    theme( strip.background = element_blank(), strip.text.x = element_blank())
)
theme_bars_tissue <- (
    theme_bars_main + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank())
)

### plot layouts 
layout_main <- rbind( c(1,2), c(1,3))
layout_individual <- rbind( c(1,1), c(1,1), c(1,1), c(1,1), c(2,2))


### PLOT FUNCTIONS 
facet_plot <- function( i_data, ylab, k = 4, theme, cols){
    gg <- (
        ggplot( i_data, aes(x = as.factor(fit2), y=eval_metric, fill = fit2)) 
            + geom_violin(width=.75)
            + geom_boxplot(width=.35, color="black", alpha=.75, fill = "white") 
            + scale_fill_manual(values=cols)
            + labs( y = ylab)
            + geom_hline(yintercept=.5,linetype="dashed", color = "grey")
            + theme
            + scale_y_continuous(breaks = seq(0,1,.1)) 
    )
    gg + facet_wrap( vars(tissue), ncol = 2)
}

barplots <- function( df, ylab = "", theme = theme_bars_main, size = 4, cols ){
    p<- (
        ggplot( df, aes(x=fit2, y=mn, fill=fit2)) 
            + geom_bar(stat="identity", position = "dodge", color = "black")
            + scale_fill_manual(values=cols)
            + facet_wrap( ~ tissue, ncol = 2)
            + labs(y = ylab, x = "")
            + theme #+ theme(legend.position = "bottom")
            + geom_text(aes(label=round(mn,2)),color="black", size=size, vjust=-.2)
            + scale_y_continuous(breaks = seq(0,1,.5), lim = c(0,1.4)) 
    )
    p  
}

pred_plots <- function( raw, bars, mod, cols) {
    
    v1 <- facet_plot( raw %>% filter( tissue == "Pan-Cancer", model == mod), ylab = title_maps$t1[[mod]], theme = theme_violins_main, cols = cols ) 
    v2 <- facet_plot( raw %>% filter( tissue %in% c("Skin", "Lung"), model == mod), ylab = "", k = 2, theme = theme_violins_tissue, cols = cols) 
    v3 <- facet_plot( raw %>% filter( tissue %in% c("Bladder", "Other"), model == mod), ylab = "", k = 2, theme = theme_violins_tissue, cols = cols)

    b1 <- barplots( bars %>% filter(tissue == "Pan-Cancer", model == mod), ylab = "", theme = theme_bars_main, size = 6, cols = cols)
    b2 <- barplots( bars %>% filter(tissue %in% c("Skin", "Lung"), model == mod), ylab = "", theme = theme_bars_tissue, size = 4, cols = cols) 
    b3 <- barplots( bars %>% filter(tissue %in% c("Bladder", "Other"), model == mod), ylab = "", theme = theme_bars_tissue, size = 4, cols = cols) 

    l1 <- arrangeGrob(v1, b1 , layout_matrix = layout_individual)
    l2 <- arrangeGrob(v2, b2 , layout_matrix = layout_individual)
    l3 <- arrangeGrob(v3, b2 , layout_matrix = layout_individual)
    l_o <- arrangeGrob(l1, l2, l3, layout_matrix = layout_main)
    l_o
}
