wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/general.R"))
source(paste0(wd,"/mission_control/helpers/figures/shapleys.R"))

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(ggborderline)

#print(sessionInfo())

ready <- readRDS(paste0(TMP_DIR, "validation-hmf-shaps.Rds")) %>% filter(feature != "Purity")

ready <- readRDS(paste0(TMP_DIR, "validation-hmf-shaps.Rds"))

importance <- (ready
    %>% group_by(model, dataset, feature) 
    %>% summarise(vals = sd(shap_feature))
    %>% ungroup()
)
importance$feature <- factor(
    importance$feature,
    levels = c("T-cell", "TMB", "Prior Systemic Therapy", "TGFB", "Proliferation")
)

p1 <- pdp_plot(ready  %>% filter(model == "lr", dataset == "Pan-Cancer (base-model)"), k = 2, ylab = "BOR Shapley Values", theme = theme_pdp_main)
p2 <- pdp_plot(ready  %>% filter(model == "lr", dataset %in% c("Skin", "Lung")), k = 1, ylab = "BOR Shapley Values", theme = theme_pdp_tissue)
p3 <- pdp_plot(ready  %>% filter(model == "lr", dataset %in% c("Bladder", "Other")), k = 1, ylab = "BOR Shapley Values", theme = theme_pdp_tissue)
p4 <- pdp_plot(ready  %>% filter(model == "os", dataset == "Pan-Cancer (base-model)"), k = 2, ylab = "OS Hazard Shapley Values", theme = theme_pdp_main)
p5 <- pdp_plot(ready  %>% filter(model == "os", dataset %in% c("Skin", "Lung")), k = 1, ylab = "OS Hazard Shapley Values", theme = theme_pdp_tissue)
p6 <- pdp_plot(ready  %>% filter(model == "os", dataset %in% c("Bladder", "Other")), k = 1, ylab = "OS Hazard Shapley Values", theme = theme_pdp_tissue)

i1 <- importance_plot(importance %>% filter(model == "lr", dataset == "Pan-Cancer (base-model)"), theme = theme_bars_main, size = 5)
i2 <- importance_plot(importance %>% filter(model == "lr", dataset %in% c("Skin","Lung")), theme = theme_bars_tissue, size = 4)
i3 <- importance_plot(importance %>% filter(model == "lr", dataset %in% c("Bladder","Other")), theme = theme_bars_tissue, size = 4)
i4 <- importance_plot(importance %>% filter(model == "os", dataset == "Pan-Cancer (base-model)"), theme = theme_bars_main, size = 5)
i5 <- importance_plot(importance %>% filter(model == "os", dataset %in% c("Skin","Lung")), theme = theme_bars_tissue, size = 4)
i6 <- importance_plot(importance %>% filter(model == "os", dataset %in% c("Bladder","Other")), theme = theme_bars_tissue, size = 4)

o1 <- arrangeGrob(p1 , i1, layout_matrix = layout_individual) 
o2 <- arrangeGrob(p2 , i2, layout_matrix = layout_individual) 
o3 <- arrangeGrob(p3 , i3, layout_matrix = layout_individual) 
l_o <- arrangeGrob(o1, o2, o3, layout_matrix = layout_main)

ll_ready <- get_dressed(l_o, title = "HMF Shapley Values - Response", size = 21)

o4 <- arrangeGrob(p4 , i4, layout_matrix = layout_individual) 
o5 <- arrangeGrob(p5 , i5, layout_matrix = layout_individual) 
o6 <- arrangeGrob(p6 , i6, layout_matrix = layout_individual) 
o_o <- arrangeGrob(o4, o5, o6, layout_matrix = layout_main)

oo_ready <- get_dressed(o_o, title = "HMF Shapley Values - Overall Survival", size = 21)

together <- arrangeGrob(ll_ready, oo_ready, ncol = 2)

options(repr.plot.width = 14)

as_ggplot(together)

saveRDS( list("ll_ready" = ll_ready, "oo_ready" = oo_ready, "together" = together), 
              paste0(FIG_DIR, "shapley_pdp_plots.Rds"))
