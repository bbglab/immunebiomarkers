
response_plot <- function( i_data, y, ylab, xlab, facet = "tcell"){
    vals <- i_data %>% pull(.data[[y]])
    gg <- (
        ggplot( i_data, aes(as.factor(bor2), y=.data[[y]], fill = as.factor(Responder))) 
            + geom_violin(width=.75)
            + geom_boxplot(width=.35, color="black", alpha=.9, fill = "white") 
            + labs( y = ylab , x = xlab)
            + ylim(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)+.3)
            + theme1
            + stat_compare_means(size = 5)
    )
    if( facet == "tissue" ){
        gg + facet_wrap( vars(tissue), ncol = 4) + ggtitle("by Tissue") + scale_x_discrete(guide = guide_axis(n.dodge=2))
    } else if( facet == "pan") {
        gg + facet_wrap( vars(pan)) + ggtitle("All samples") + scale_x_discrete(guide = guide_axis(n.dodge=2))
    } else if ( facet == "tcell"){
        gg + facet_wrap( vars(tcell2))
    } 
}

together <- function( i_data, y, ylab, xlab, plt = "main"){
    a <- response_plot( i_data, y, ylab, xlab, facet = "pan") 
    b <- response_plot( i_data %>% filter(tissue != "Skin"), y, ylab, xlab, facet = "tissue")
    c <- response_plot( i_data, y, ylab, xlab, facet = "tcell")
    
    if( plt == "main" ){
        arrangeGrob(a + theme1, b + theme2, layout_matrix = rbind(c(1,1,2,2,2,2), c(1,1,2,2,2,2)))
    } else if (plt == "tgfb") {
        arrangeGrob(a + theme1, c + theme2, ncol = 2)    
    }
}

ggsurvplot_facet2 <- function(pval.size = 5, ...)
{
  newcall <- bquote(
    p <- p + geom_text(data = pvals.df, aes(x = pval.x, y = pval.y, 
    label = pval.txt), size = .(pval.size), hjust = 0)
    )

    body(ggsurvplot_facet)[[20]][[3]][[8]] <- newcall
    ggsurvplot_facet(...)
}
