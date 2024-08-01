library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/xgboost/shapley_example.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure5_themes.R"))

validation_ready <- readRDS(paste0(TMP_DIR, "validation-hmf-go.Rds"))
hmf_loo <- readRDS(paste0(TMP_DIR, "validation-loo-cv.Rds")) #%>% drop_na(tcell)
start <- validation_ready %>% left_join(hmf_loo, by = "patient_id") %>% drop_na(pred_lr) %>% drop_na(tcell)

color_map <- list(
    "Prior Systemic Therapy" = '#FFFFB3',
    "TGFB" = "#BEBADA",
    "Proliferation" = "#8DD3C7",
    "T-cell" = '#FB8072',
    "TMB" = '#80B1D3'
)
base_theme <- theme_classic(base_size = 22) 

theme_shaps <- (
   base_theme + 
   theme(axis.title.y=element_blank(), 
         legend.position = "none",
         axis.text.y=element_blank(), 
         axis.ticks.y=element_blank(),
         plot.title = element_text(hjust = 0.5, size = 22),
         axis.text.x = element_text(size = 22), 
         plot.margin = unit(c(.5,.5,.5,.5), "cm")
        )
)
theme_patient <- (
    theme_shaps + 
    theme(plot.margin = unit(c(0,.5,0,.5), "cm"),
          axis.title.x=element_blank())
)

theme_main <- (
    base_theme + 
    theme( legend.position = "none", 
           plot.title = element_text(hjust = 0.5, size = 16),
           plot.margin = unit(rep(.2,4), "cm"),
           axis.ticks.x=element_blank()#,axis.title.x=element_blank()
)
)

my_palette<-brewer.pal(6,"RdYlGn")[c(1,3,6)]

lr <- (
    start
        %>% rename( pred = pred_lr, dataset = tissue, tissue = tissue_full )
        %>% select(
                patient_id, 
                dataset,
                tissue,
                age, 
                gender, 
                contains("shap_lr"),
                #-shap_lr_pretreat,
                -shap_lr_BIAS, 
                #-shap_lr_pdl1,
                -shap_lr_purity, 
                pred, 
                bor,
                os
        )
        %>% rename( "T-cell" = shap_lr_tcell, "TMB" = shap_lr_tmb, "TGFB" = shap_lr_tgfb, 
                    "Proliferation" = shap_lr_prolif, "Prior Systemic Therapy" = shap_lr_pretreat)
)
lr$model <- "lr"
ready_lr <- (
    lr 
        %>% gather("feature", 
                   "shap_feature", 
                   "TMB", 
                   'T-cell', 
                   "Proliferation", 
                   'TGFB', 
                   'Prior Systemic Therapy'
                  )
)

os <- (
    start
        %>% rename( pred = pred_os, dataset = tissue, tissue = tissue_full )
        %>% select(
                patient_id, 
                dataset,
                tissue,
                age, 
                gender, 
                contains("shap_os"), 
                #-shap_os_pretreat,
                -shap_os_BIAS, 
                #-shap_os_pdl1,
                -shap_os_purity,
                pred, 
                bor,
                os
        )
        %>% rename( "T-cell" = shap_os_tcell, "TMB" = shap_os_tmb, "TGFB" = shap_os_tgfb, 
                    "Proliferation" = shap_os_prolif, "Prior Systemic Therapy" = shap_os_pretreat)
)
os$model <- "os"
ready_os <- (
    os 
        %>% gather("feature", 
                   "shap_feature", 
                   "TMB", 
                   'T-cell', 
                   "Proliferation", 
                   'TGFB', 
                   'Prior Systemic Therapy'
                  )
)

ready <- rbind(ready_lr, ready_os)

ready$feature <- factor(
    ready$feature,
    levels = rev(c("T-cell", "TMB", "Prior Systemic Therapy", "TGFB", "Proliferation"))
)

ready$col <- "Medium"
ready$col <- ifelse(ready$pred < .1 & ready$model == "lr", "Low", ready$col)
ready$col <- ifelse(ready$pred > .5 & ready$model == "lr", "High", ready$col)

ready$col2 <- "Medium"
ready$col2 <- ifelse(ready$pred < .5 & ready$model == "os", "Low", ready$col2)
ready$col2 <- ifelse(ready$pred > 1.5 & ready$model == "os", "High", ready$col2)

ready$more <- "all"

get_examples <- function( lr_low, lr_high, os_low, os_high){
    (ready 
    %>% filter(model == "os", 
               patient_id %in% (
                   ready 
                       %>% filter(model == "lr", feature == "TMB", pred > lr_low, pred < lr_high) 
                       %>% pull(patient_id)
                ))
    %>% filter( pred > os_low, pred < os_high, feature == "TMB" )
    %>% arrange(pred)
    )
}

ll <- get_examples( lr_low = 0,lr_high = .1, os_low = 0, os_high = .5)
lh <- get_examples( lr_low = 0,lr_high = .1, os_low = 3, os_high = 6)
hh <- get_examples( lr_low = .65,lr_high = 1, os_low = 0, os_high = .5)
mm <- get_examples( lr_low = .2,lr_high = .3, os_low = .7, os_high = 1)

#mm

examples <- list( "ll" = "CPCT02140093", "lh" = "CPCT02010333", "hh" = "CPCT02020853", "mm" = "CPCT02020894")

pps <- list()
for( i in names(examples)){
    pps[[i]] <- as_ggplot(patient_plots( ready, examples[[i]] )) + theme(plot.margin = unit(c(8,0,0,0), "cm"))
}

start$highlight <- ifelse( start$patient_id %in% unname(unlist(examples)), "yes", "no")

main <- (
    ggplot(data=start, aes(x=pred_lr, y=pred_os, fill = lr_gp, shape = os_gp, size = highlight, alpha = highlight)) 
    + geom_point(color = "black") 
    + theme_main
    + scale_fill_manual(values = c('Low' = my_palette[1], 'Medium' = my_palette[2],'High' = my_palette[3]))
    + scale_shape_manual(values = c(21,21,21))
    + scale_size_manual(values = c(3,13))
    + scale_alpha_manual(values = c(.4,1))
    + labs( x = "Probability of Response", y = "OS Hazard", title = "OS Hazard vs Probability of Response")
    + geom_vline(xintercept = .1,linetype = "dotted")
    + geom_vline(xintercept = .5,linetype = "dotted")
    + geom_hline(yintercept = .5,linetype = "dotted")
    + geom_hline(yintercept = 1.5, linetype = "dotted")
    + scale_x_continuous(breaks=c(.1,.5,1), lim = c(0,1.05), labels = scales::percent_format(scale = 100))
    + scale_y_continuous(breaks=c(.5,1.5,3, 4.5,6,10), lim = c(0,6.5))
)

saveRDS( list("main" = main, "exs" = pps), file = paste0(FIG_DIR, "figure5_g.Rds"))
