wd <- dirname(dirname(getwd()))
setwd(wd)

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/figure2_theme.R"))

suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
labels <- readRDS(paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
threshold <- -log10(ingredients$by_05_fdr)[1]

geneset_names = data.table::fread(paste0(REF_DIR, "CPI_genesets.txt"))
colnames(geneset_names)[1] = "feature"

go <- (
    ingredients 
      %>% filter(feature != "clinical_systemic_composite") 
      %>% filter(!grepl("ciber", feature))
      %>% filter(Group == "Gene Set",dataset == "all", model == "bor",covariates == "age_biopsy_purity_tissue")
      %>% mutate(selected = ifelse(cor_tcell > .8, TRUE, FALSE))
      %>% mutate(feature = sub("isofox_", "", feature))
      %>% left_join(geneset_names, by = "feature")
      %>% filter( !grepl("cluster", feature), !grepl("ariathan", feature), !grepl("vhio", feature), !grepl("rand", feature))
)

fig_2d <- (
  ggplot(go, aes(x = cor_tcell, y = log10_p, col = selected)) +
    geom_point(size = 3) +
    xlab("Gene set correlation with expression cluster") +
    ylab("-Log10 (p-value)") +
    ggtitle("Cluster 3 represents T-cell infiltration") +
    scale_color_manual(values=c("darkgrey","#FB8072")) +
    geom_vline(xintercept = 0.8, linetype="dashed", color = "black", linewidth =.1) +
    geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth =.2) +
    geom_text_repel(aes(label = if_else((cor_tcell > 0.8 & log10_p > threshold), Short_name, "")), size=2.7, col="grey20" ,max.overlaps=10)+
    scale_y_continuous(breaks = c(2,4,6,8), lim = c(0,8)) + 
    scale_x_continuous(breaks = c(0,.5,.8,1)) +  
    base_theme
)

saveRDS(fig_2d, file = paste0(FIG_DIR, "fig_2d.Rds"))
