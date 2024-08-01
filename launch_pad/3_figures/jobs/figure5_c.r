wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

library(survival)
library(survminer)
library(survcomp)

base_theme <- theme_classic(base_size = 31) 
theme <- (
    base_theme + 
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position="none", 
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 15))
)
remove <- theme(axis.title.y=element_blank(), 
                axis.text.y=element_blank())
remove_x <- theme(axis.title.x=element_blank())

fig5_data <- readRDS(paste0(TMP_DIR, "figure5_data.Rds"))

fits_lr <- list()
fits_os <- list()
for( i in names(fig5_data)){
    fits_lr[[i]] <- survfit(Surv(os_days/(365/12), os_event) ~ lr_gp, data = fig5_data[[i]])
}
for( i in names(fig5_data)){
    fits_os[[i]] <- survfit(Surv(os_days/(365/12), os_event) ~ os_gp, data = fig5_data[[i]])
}

ggsurvplot_facet2 <- function(pval.size = 5, ...)
{
  newcall <- bquote(
    p <- p + geom_text(data = pvals.df, aes(x = pval.x, y = pval.y, label = pval.txt), size = .(pval.size), hjust = 0)
    )

    body(ggsurvplot_facet)[[20]][[3]][[8]] <- newcall
    ggsurvplot_facet(...)
}
cols <- brewer.pal(6,"RdYlGn")[c(1,3,6)]
cols2 <- brewer.pal(6,"RdYlGn")[c(6,1,3)]

go_surv <- function(pts, i) {
    (
        ggsurvplot_facet2(pts[[i]], 
                          fig5_data[[i]], 
                          ggtheme = theme, 
                          xlab = "Months", 
                          ylab = "Survival Probability",
                          facet.by = "tissue", 
                          pval = TRUE, 
                          palette = cols2,
                          short.panel.labs=T,
                          pval.coord = c(15,.94), 
                          pval.size = 7) + theme(legend.position="none")
    )
}

plts_lr <- list()
plts_os <- list()
for( i in names(fig5_data)){
    plts_lr[[i]] <- go_surv(fits_lr, i)
    plts_os[[i]] <- go_surv(fits_os, i) #+ scale_x_continuous(breaks = c(0,10,20,50), limits = c(0,55))
}

plots_lr <- arrangeGrob(plts_lr$'HMF-CPCT', plts_lr$Melanoma + remove, plts_lr$Lung + remove, 
                         plts_lr$Bladder + remove, plts_lr$Other + remove, 
                         layout_matrix = rbind( c(1,1,2,3), c(1,1,4,5)))

plots_os <- arrangeGrob(plts_os$'HMF-CPCT' +remove_x, plts_os$Melanoma + remove, plts_os$Lung + remove, 
                         plts_os$Bladder + remove, plts_os$Other + remove, 
                         layout_matrix = rbind( c(1,1,2,3), c(1,1,4,5)))

options(repr.plot.width=11, repr.plot.height=7, resolution = 500)
plot(plots_lr)

saveRDS( list( "response" = plots_lr , "survival" = plots_os), file = paste0(FIG_DIR, "figure5_c.Rds") )
