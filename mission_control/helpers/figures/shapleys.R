
### THEME 
base_theme <- theme_classic(base_size = 21) 

theme_pdp <- (
    base_theme + 
    theme(axis.title.x=element_blank(), 
          legend.position="none",
          plot.margin = unit(c(0,0,0,0), "cm"))
)

theme_pdp_main <- theme_pdp
theme_pdp_tissue <- (
    theme_pdp + 
    theme(axis.text.x=element_blank(), 
          axis.title.y=element_blank(), 
          axis.text.y=element_blank())
)

theme_bars <- (
    theme_pdp + 
    theme ( axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            strip.background = element_blank(), 
            strip.text.x = element_blank())
)

theme_bars_main <- theme_bars 
theme_bars_tissue <- (
    theme_bars + 
    theme(axis.text.x=element_blank(), 
          axis.title.y=element_blank(), 
          axis.text.y=element_blank())
)

### plot layouts 
layout_main <- rbind( c(1,2), c(1,3))
layout_individual <- rbind( c(1,1), c(1,1), c(1,1), c(1,1), c(2,2))


### PLOT FUNCTIONS 
pdp_plot <- function( df, ylab, k = 2, theme){
    gg <- (
    df %>% 
      ggplot(aes(scaled_feature, y = shap_feature, color = feature, group = feature)) + 
      scale_color_manual(values=unlist(color_map)) + 
      geom_borderline(size = k, bordercolour = "black") + 
      facet_wrap( ~ dataset, ncol = 2) + 
      labs(x = "Scaled Feature Values") + 
      theme +
      scale_x_continuous( breaks = seq(-3,3,3), limits = c(-3,3))
    ) 
    if( grepl("BOR",ylab)) {
      gg <- gg + scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.56,1.56) ) + labs(y = ylab)  
    } else {
      gg <- gg + scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.39,1.39)) + labs(y = ylab)    
    }
    gg
}

importance_plot <- function( df , theme, size ) {
    p <- (
        ggplot(df, aes(x=feature, y = vals, fill = feature))     
            + geom_bar(stat="identity", position = "dodge", color = "black")
            + scale_fill_manual(values=unlist(color_map))
            + facet_wrap( ~ dataset, ncol = 2) 
            + labs(y = "", x = "")  
            + theme
            + geom_text(aes(label=round(vals,1)),color="black", size=size, vjust=-.2)
    )
    if( "Pan-Cancer (base-model)" %in% df$dataset){
        p + scale_y_continuous(breaks = 0, lim = c(0,max(3*df$vals)))  
    } else {
        p + scale_y_continuous(breaks = 0, lim = c(0,max(1.5*df$vals)))  
    }
}
