
i <- 5; k <- .3; l <- .4

base_theme <- theme_classic(base_size = 21+i) 

theme1 <- (
    base_theme + 
    theme(plot.margin = unit( rep(k,4), "cm"), 
          legend.position = "none",
           axis.title.x=element_blank(),
           plot.title = element_text(size = 18+i, hjust = .5),
           axis.text = element_text(size = 18+i)
         )
) 

theme2 <-(
    theme1 + 
    theme(axis.title.y=element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 13+i))
)
surv_theme1 <- theme1 + theme(axis.title.x = element_text(size = 21+i))
surv_theme2 <- (
    surv_theme1 + 
    theme(axis.title.y=element_blank(), 
          axis.text.y = element_blank()
          )
)
