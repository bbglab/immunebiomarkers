wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))

library(survival)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

go <- readRDS(paste0(TMP_DIR, "forest-ready.Rds"))  %>% filter(clean_feature != "Purity")

base_theme <- theme_classic(base_size = 21) 

theme <- (
    base_theme + 
    theme(
        axis.title.y=element_blank(), 
        plot.title = element_text(hjust = 0.5, size =25),#, size = 21),
        legend.position="none"
    )
)

my_palette <- brewer.pal(6,"RdYlGn")[c(1,6)]

color_map <- list(
    "yes" = my_palette[2],
    "no" = my_palette[1],
    "non" = "grey"
)
alpha_map <- list(
    "strong" = 1,
    "moderate" = 1,
    "weak" = .4,
    "non" = .8
)
feature_color_map <- list(
    "TMB" = '#80B1D3',
    "T-cell" = '#FB8072',
    "TGFB" = "#BEBADA",
    "Proliferation" = "#8DD3C7",
    "Pretreatment" = "#9E7788"
)

plot_maker <- function( df, title = "HMF CPCT Study", study = "clean_study2" ){
    
    if( study == "clean_study"){
        df <- df %>% filter( sets == "clusters")
    }
    
    ggplot( df,  
        aes_string(x = study, 
            y= "est", 
            ymin= "ci_low", 
            ymax= "ci_high", 
            linetype = "sets",
            color = "better", 
            alpha = "z_group")
      ) +
    geom_pointrange(lwd = 1.2) + 
    geom_hline(yintercept=0, lty=2, col = "grey") + 
    coord_flip() +  
    facet_grid( rows = vars(clean_feature), cols = vars(clean_model), scales = "free") + 
    scale_color_manual(values = unlist(color_map)) +
    scale_alpha_manual(values = unlist(alpha_map)) +
    scale_y_continuous(n.breaks = 4 )  +
    ggtitle(title) + 
    ylab("Coefficient Estimate (95% CI)") + 
    theme  
}

save_the_forests <- list()

for (i in c("multi", "single")){
    ready <- if( i == "multi"){ go %>% filter(features == "all" )} else {go %>% filter(features != "all" )}
    for( j in c("clean_study", "clean_study2")){
        save_the_forests[[i]][[j]][['cpct']] <- (
            plot_maker( ready %>% filter( grepl("CPCT", cohort) ),  title = "HMF CPI", study = j))
        save_the_forests[[i]][[j]][['external']] <- (
            plot_maker( ready %>% filter( ! grepl("CPCT", cohort) ),  title = "Validation Cohorts", study = j))
    }
}

saveRDS( save_the_forests, paste0(FIG_DIR, "figure4-and-supplement-forest-plots.Rds"))

color_strips <- function(gg){
  g <- ggplot_gtable(ggplot_build(gg))
  strips <- which(grepl('strip', g$layout$name))
  pal <- c("white", "white",unlist(feature_color_map))
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i] 
  } 
  as_ggplot(g)
}

layout_matrix = rbind(c(rep(1,11), rep(2,12)))

plot_together <- function( p1, p2, layout = layout_matrix, title ){
    step1 <- as_ggplot(
        arrangeGrob(color_strips(p1) + theme(plot.margin = margin(0,1,0,1, "cm")), 
                    color_strips(p2) + theme(plot.margin = margin(0,1,0,1, "cm"))
    , layout_matrix = layout_matrix))
    
    ready <- ( step1 
        + theme(plot.margin = margin(0,0,0,0, "cm"), plot.title = element_text(hjust = 0.5, size =25)) 
        + ggtitle(title)
    )
    ready
}

options(repr.plot.width = 20, repr.plot.height= 12, resolution = 200)
figure4 <- plot_together( save_the_forests$multi$clean_study$cpct, save_the_forests$multi$clean_study$external,
               title = "Multivariate Model Estimates")
figure4
ggsave( paste0(FIG_FINAL_DIR, "figure4.png"), width = 20, height = 12)

#FIG_FINAL_DIR

options(repr.plot.width = 20, repr.plot.height= 20, resolution = 200)
plot_together( save_the_forests$multi$clean_study2$cpct, 
               save_the_forests$multi$clean_study2$external, 
               title = "Multivariate Model Estimates")
ggsave( paste0(FIG_FINAL_DIR, "7_forest_plot_sets_clusters.png"), width = 20, height = 20)

options(repr.plot.width = 20, repr.plot.height= 12, resolution = 200)
forest_uni <- plot_together( save_the_forests$single$clean_study$cpct, 
               save_the_forests$single$clean_study$external, 
               title = "Univariate Model Estimates")
saveRDS( forest_uni, paste0(FIG_DIR, "supplement_forest_univariate.Rds"))
#ggsave( paste0(FIG_FINAL_DIR, "7_forest_plot_supp_univariate.png"), width = 20, height = 12)
