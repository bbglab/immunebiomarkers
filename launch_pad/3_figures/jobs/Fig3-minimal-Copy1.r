wd <- dirname(dirname(getwd()))
setwd(wd)

source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_plots.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
source(paste0(wd,"/mission_control/helpers/figures/themes.R"))

.libPaths("/home/jusset/anaconda3/envs/biomarker_r/lib/R/library")

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggh4x)
library(lemon)

ready <- readRDS(paste0(I_DIR, "cpi_go.Rds"))
results <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
cutoff <- results$by_05_fdr[1]

one_set_result<- results %>% filter(model == "os" & dataset == "all" & covariates == "age_biopsy_purity_tissue") 
os_resid_set_result<- results %>% filter(model == "os" & dataset == "all" & covariates == "residuals") 

# drop features are duplicated features
drop_features <- c("tcell", "tmb", "prolif","tgfb","pretreat","pretreat_comp","tissue","tissue_full","age","biopsy","purity","pdl1")
#There 2 significant features not expression - check!!! For now they are excluded
weird_features <- c("cnv.region_loh_chr11.q23.3","somatic.gene_MATN2.mb")
top_os_resid <- results %>%
    filter(model == "os", dataset == "all", covariates == "residuals") %>%
    filter(!feature %in% drop_features) %>%
    filter(!feature %in% weird_features) %>%
    filter(!grepl("gene_set_", feature)) %>%
    filter(!grepl("cibersort", feature)) %>%
    filter(p_val < cutoff) %>%
    pull(feature)

os_resid_cor_go <- ready %>% select(all_of(top_os_resid))  %>% drop_na() 
#cor_os_resid_with_pval <- rcorr(as.matrix(os_resid_cor_go))
#cor_matrix_bor <- ifelse((cor_bor_with_pval$P < 0.05|is.na(cor_bor_with_pval$P)), cor_bor_with_pval$r, 0)
cor_matrix_os_resid = cor(as.matrix(os_resid_cor_go))

heatmap_genes <- cor_matrix_os_resid
heatmap_genes_group = merge(as.data.frame(heatmap_genes),one_set_result[c("feature","Type","big_group")],by.y = "feature", by.x = 0)

#tree_clust

# Generate clusters
tree <- hclust(as.dist(1-abs(cor_matrix_os_resid)),method = "ward.D2")
tree_clust = cutree(as.hclust(tree), k = 3)
cluster1 = tree_clust[tree_clust == 1]
cluster2 = tree_clust[tree_clust == 3]
cluster3 = tree_clust[tree_clust == 2]
clusters = c(cluster1, cluster2, cluster3)

#clusters

# save the names of the genes in the cluster
df1 = data.frame(T_cell = names(cluster1))
write.table(df1, "1_figures//revised_plots_Axel/Prolif_cluster_genes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

df2 = data.frame(T_cell = names(cluster2))
write.table(df2, "1_figures//revised_plots_Axel/TGBF_cluster_genes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)       

wdf <- one_set_result[,c("feature", "big_group")] %>% filter(feature %in% top_os_resid)
wdf$w = 3
# TODO match order of the figures with the order of the clusters!!!!!
ord = match(names(clusters), wdf$feature)
wdf = wdf[tree$order,]
wdf$clustord = 1:nrow(wdf)
wdf = wdf %>% arrange(clustord)
breaks <- wdf %>% 
  mutate(cumw = cumsum(w),
         pos = .5 * (cumw + lag(cumw, default = 0))) %>% 
  select(feature, pos)

os_resid_cor_go = os_resid_cor_go[match(wdf$feature, names(clusters)),]
pd <- as.data.frame(cor_matrix_os_resid)
pd$feature <- rownames(pd)
pd.m <- pivot_longer(pd, -feature, names_to = "variable")

pd.m <- pd.m %>% 
  left_join(breaks, by = c("variable" = "feature")) %>% 
  rename(y = pos) %>% 
  left_join(breaks, by = c("feature" = "feature")) %>% 
  rename(x = pos)

pd.m <- pd.m %>% 
  left_join(wdf, by = c("feature" = "feature")) 

pd.m <- pd.m %>% 
  left_join(wdf[,c(1,2)], by = c("variable" = "feature")) 

RNA_features = wdf[wdf$big_group == "RNA",]$feature

pd.m$feature <- factor( pd.m$feature, levels = colnames(os_resid_cor_go)[ord] )
pd.m$variable <- factor(pd.m$variable, levels = colnames(os_resid_cor_go)[ord] )
pd.m$height = ifelse((
     pd.m$variable=="clinical_meta_hasSystemicPreTreatment2"|
     pd.m$variable=="clinical_pre_treated"|
     pd.m$variable=="clinical_systemic_composite"), 
    35,ifelse(pd.m$variable %in% RNA_features, 3,1))

pd.m$width =ifelse(
    (pd.m$feature=="clinical_meta_hasSystemicPreTreatment2"|
     pd.m$feature=="clinical_pre_treated"|
     pd.m$feature=="clinical_systemic_composite"), 
    35,ifelse(pd.m$feature %in% RNA_features,3,1))

pd.m$Type = ifelse(pd.m$big_group.x == "RNA", "Expression", ifelse(pd.m$big_group.x == "Somatic", "SNVs/Indels", "Clinical"))
pd.m$Type2 = ifelse(pd.m$big_group.y == "RNA", "Expression", ifelse(pd.m$big_group.y == "Somatic", "SNVs/Indels", "Clinical"))

pd.m = pd.m %>% 
    mutate(Type = case_when(variable %in% names(cluster1) ~ "cluster1",
                            variable %in% names(cluster2) ~ "cluster2",
                            variable %in% names(cluster3) ~ "cluster3",
                            !clusters[feature] %in% c(1,2) ~ big_group)) %>% 
    mutate(Type = factor(Type, levels = c("cluster1", "cluster2", "cluster3")))

pd.m = pd.m %>% 
    mutate(Type2 = case_when(feature %in% names(cluster1) ~ "cluster1",
                          feature %in% names(cluster2) ~ "cluster2",
                          feature %in% names(cluster3) ~ "cluster3",
                          !clusters[feature] %in% c(1,2) ~ big_group)) %>% 
    mutate(Type2 = factor(Type2, levels = c("cluster1", "cluster2", "cluster3")))

colnames(pd.m)[3]<-"correlation"

fig3b <- ggplot(pd.m, aes(x = x, y = y)) +
  geom_tile(aes(width = width,  height = height, fill = correlation)) +
  scale_x_reverse(breaks = breaks$pos, labels = breaks$feature, expand = c(0, 0.1)) +
  scale_y_continuous(breaks = breaks$pos, labels = breaks$feature, expand = c(0, 0.1))+
  scale_fill_gradientn(limits = c(-1, 1), colours = c("blue", "white", "white",  "red"), values = c(0, 0.46, 0.54, 1)) +
  theme(panel.border=element_rect(fill = NA, colour="black",linewidth=5))+
  theme_bw(base_size = 14) +
  facet_grid2(Type ~ Type2, scales = "free", space = "free", switch = "both",render_empty = TRUE,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8DD3C7","#BEBADA","lightgrey"), size= 1),
                                    background_y = elem_list_rect(fill = c("#8DD3C7","#BEBADA", "lightgrey"), size= 1)))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        panel.spacing.x = unit(0, "mm"),
        panel.spacing.y = unit(0, "mm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "black", linewidth = 10),
        panel.grid = element_blank(),
        legend.position="right") + 
    xlab("") + ylab("")

fig3b <- annotate_figure(fig3b, fig.lab = "B", fig.lab.size = 20, fig.lab.face = "bold")
fig3b

# read more interpretable geneset names: 
geneset_names = data.table::fread(paste0(REF_DIR, "CPI_genesets.txt"))
colnames(geneset_names)[1] = "feature"
geneset_names = sapply(geneset_names, \(x) gsub("gene_set_", "", x)) |> as.data.frame()

read_features = ready %>% select(any_of(unique(ingredients3$feature)))

# get the features of each cluster
features_cl3 = wdf %>% filter(big_group == "RNA") %>% pull(feature)

mean_cl1 = rowMeans(read_features %>% select(all_of(names(cluster1))))
mean_cl2 = rowMeans(read_features %>% select(all_of(names(cluster2))))

cluster_cors = tibble(feature = colnames(read_features),
                   `Expression cluster 1` = cor(mean_cl1, read_features, use = "pairwise.complete.obs") %>% as.numeric(),
                   `Expression cluster 2` = cor(mean_cl2, read_features, use = "pairwise.complete.obs") %>% as.numeric())

cluster_cors = cluster_cors %>% filter(feature %in% unique(ingredients3$feature))
i_cors = left_join(ingredients3, cluster_cors,by = "feature" )

# merge and filter the labels to match value: 
label_cors = left_join(ingredients3, cluster_cors,by = "feature" ) %>% 
    filter(log10_p > threshold)

label_cors = label_cors %>% 
    mutate(big_group = case_when(
        grepl("gene_set", feature) & `Expression cluster 1` > 0.75 & log10_p > threshold ~ "cluster 1 spec. gene set",
        grepl("gene_set", feature) & `Expression cluster 2` > 0.75 & log10_p > threshold ~ "cluster 2 spec. gene set",
        .default = big_group))

label_cors = label_cors %>% 
    select(starts_with("Expr"), big_group, log10_p, feature) %>% 
    pivot_longer(cols = starts_with("Expr"), names_to = "cluster", values_to = "correlation") %>% 
    filter(correlation > 0.85) %>% 
    filter(grepl("gene_set", feature)) %>% 
    filter(!grepl("rand", feature))

# merge the label names with simple names, if they are avaialble
label_cors$label = gsub("gene_set_|isofox_", "", label_cors$feature)
label_cors$label = geneset_names$Short_name[match(label_cors$label, geneset_names$feature)]

i_cors = i_cors %>% 
    filter(big_group == "RNA") %>% 
    mutate(big_group = case_when(
        grepl("gene_set", feature) ~ "gene set",
        .default = big_group)) %>% 
    mutate(big_group = case_when(
        grepl("gene_set", feature) & `Expression cluster 1` > 0.75 & log10_p > threshold ~ "cluster 1 spec. gene set",
        grepl("gene_set", feature) & `Expression cluster 2` > 0.75 & log10_p > threshold ~ "cluster 2 spec. gene set",
        .default = big_group))

input_cors = i_cors %>%  
    filter(model == "os" & dataset == "all" & covariates == "residuals") %>% 
    select(starts_with("Expression"), big_group, log10_p, feature) %>% 
    pivot_longer(cols = starts_with("Expression"), names_to = "cluster", values_to = "correlation")

fig3c = (ggplot(input_cors, aes(x = correlation, y = log10_p,
                           size = correlation > 0.5, alpha = correlation > 0.5)) + 
    geom_point(aes(fill = big_group), shape = 21, color = "white", stroke = 0.3) + 
    facet_rep_grid(. ~ cluster ) 
  +  geom_label_repel(data = label_cors, aes(label = label, fill = big_group),
                      min.segment.length =  0.2, alpha = 0.85, size = 4, show.legend = FALSE)
  + exhaustive_theme_left_legend
  + scale_y_continuous(breaks = c(seq(0,8,2)), limits = c(0,8)) 
  + scale_x_continuous(breaks = c(-.5, 0, 1), limits = c(-.5,1))  
  + scale_size_manual(values = c(1, 2)) 
  + scale_alpha_manual(values = c(0.5, 1))
  + scale_fill_manual(values = c("#8DD3C7", "#BEBADA","#B3DE69",  "#FC913A") )
  + geom_vline(xintercept = 0.5, linetype="dashed", color = "black", linewidth = .1)
  + geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth = .1)
  + guides(alpha = "none", color = "none", size = "none", # remove the 'pval size'
           fill = guide_legend(title = "Group", override.aes = list(size = 4)))
  + theme(strip.background = element_blank())
  + xlab("Correlation to cluster mean") + ylab("-Log10 p-value")
)
fig3c

fig3c <- annotate_figure(fig3c, fig.lab = "C", fig.lab.size = 20, fig.lab.face = "bold")
ggsave("1_figures/revised_plots_Axel/plots/figure3c.png", fig3c, width = 13, height = 5)

cl3_cors_to_meancl2 = input_cors  %>%
filter(feature %in% names(cluster3),
    cluster == "Expression cluster 2")  %>% 
    mutate(`cluster origin` = "cluster 3")

cl2_cors_to_meancl2 = input_cors  %>%
    filter(feature %in% names(cluster2),
    cluster == "Expression cluster 2")  %>% 
    mutate(`cluster origin` = "cluster 2")

cors_cl2_3 = rbind(cl3_cors_to_meancl2, cl2_cors_to_meancl2)

cor_to_cl2 = ggplot(cors_cl2_3, aes(x = `cluster origin`, y = correlation, fill = `cluster origin`)) + 
    geom_boxplot(outlier.shape = NULL) +
    geom_jitter(width = 0.3) + 
    theme_bw() + 
    scale_fill_manual(values = c("#BEBADA", "lightgrey")) + 
    xlab("") + ylab("Correlation genes to mean expression cluster 2")
cor_to_cl2
ggsave("../plots/Supplementary_figureX.pdf", cor_to_cl2, width = 4.5, height = 4.5)
ggsave("../plots/Supplementary_figureX.png", cor_to_cl2, width = 4.5, height = 4.5)

ingredients_selected <- ingredients3 %>%
    filter(abs(cor_pretreat) < 0.3 & abs(cor_tmb) < 0.2 & abs(cor_tcell)<0.3 & abs(cor_prolif)<0.3 & abs(cor_tgfb)<0.3)

color_map_type$RNA = "#FC913A" # change the RNA feature value to a more direct value. 
color_map_type$`CNV/SVs` <- "#B3DE69"
color_map_type$Somatic <- "#80B1D3"
color_map_type$CNV <- color_map_type$SVs <- NULL

# example script to cut out the functions
fig3d = ggplot(ingredients_selected, aes(x = plot_est, y = log10_p, fill = big_group,
                                 size = log10_p > threshold, alpha = log10_p > threshold)) + 
  geom_point(shape = 21, color = "white", stroke = .5) + 
  geom_vline(xintercept = 1, linetype="dashed", color = "black", linewidth = .1) +
  geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth = .2) +
  scale_fill_manual(values = as.vector(color_map_type)) + 
  scale_size_manual(values = c(2, 2.5)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = c(0,5,10), lim = c(0,10))  +
  scale_x_continuous(breaks = c(1,2), limits = c(.4,2.5)) + 

  guides(alpha = "none", color = "none", size = "none", # remove the 'pval size'
           fill = guide_legend(override.aes = list(size=4))) +
  exhaustive_theme_left_legend + 
  theme(legend.position = c(0.85, 0.25), legend.background = element_blank(),
          legend.box.background = element_blank()) + 
 
  labs(fill = "Group") + 
  xlab("1 / OS Hazard Estimate") + 
  ylab("-Log10 p-value") + 
  ggtitle("OS residuals vs remaining features")
fig3d

fig3d <- annotate_figure(fig3d, fig.lab = "D", fig.lab.size = 20, fig.lab.face = "bold")
ggsave("/workspace/users/arosendahl/plots/Fig3E.png", fig3d, width = 7, height = 7, dpi = 300)

library(grImport2)
library(rsvg)
rsvg_svg("/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_5.svg", "/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_5_cairo.svg")

SVGlogo <- readPicture("/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_5_cairo.svg")
svg_pic = pictureGrob(SVGlogo, ext="gridSVG", delayContent=FALSE)
svg_pic = annotate_figure(svg_pic, fig.lab = "E", fig.lab.size = 20, fig.lab.face = "bold")

pl = list(fig3a, fig3b, fig3c, fig3d)
pl = lapply(pl, "+", theme(plot.margin=margin(2,2,2,2, unit = "mm")))

lay <- rbind(c(1,2),
             c(3,3),
             c(4,5))

fig <- arrangeGrob(pl[[1]], pl[[2]],pl[[3]], pl[[4]],svg_pic,  layout_matrix=lay)  

#ggsave("/home/mandrianova/therapy_biomarkers/immune_biomarkers/1_figures/review_new/Final_plots/Fig3.png", fig, width = 17, height = 20, dpi = 300)
ggsave("1_figures//revised_plots_Axel/plots/Fig3.pdf", fig, width = 15, height = 17, dpi = 300)
ggsave("1_figures//revised_plots_Axel/plots/Fig3.png", fig, width = 15, height = 17, dpi = 300)
