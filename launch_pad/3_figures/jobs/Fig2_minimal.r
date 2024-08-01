wd <- dirname(dirname(getwd()))
setwd(wd)

list.files()

source(paste0(wd,"/mission_control/treasure_map.R"))
#source(paste0(wd,"/mission_control/helpers/figures/exhaustive_plots.R"))
#source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
source(paste0(wd,"/mission_control/helpers/figures/themes.R"))

.libPaths("/home/jusset/anaconda3/envs/biomarker_r/lib/R/library")

library(tidyverse)
library(ggrepel)
library(data.table)
library(gridExtra)
library(ggpubr)
library(Hmisc)
library(ggh4x)
library(lemon)

ingredients <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
labels <- readRDS(paste0(TMP_DIR,"exhaustive-plots-labels.Rds"))
threshold <- -log10(ingredients$by_05_fdr)[1]

# filter the systemic clinical composite part from the ingredients file
# pool CNVs and SV's in the same category (they are not very important in this context)
ingredients = ingredients %>%
    filter(feature != "clinical_systemic_composite") %>% 
    filter(!grepl("cibersort", feature)) %>% 
    mutate(big_group = case_when(big_group == "CNV" ~ "CNV/SVs",
                                 big_group =="SVs" ~ "CNV/SVs", .default = big_group))

labels$clean_label <- ifelse(labels$feature %in% c("isofox_gene_set_t_cell_effector"), "T-cell effector gene set", as.character(labels$clean_label))
labels$clean_label <- ifelse(labels$feature %in% c("clinical_meta_hasSystemicPreTreatment2"), "Prior systemic therapy", as.character(labels$clean_label))

groups_to_labels = data.frame(names=c("Clinical","HLA","RNA", "Somatic", "CNV/SVs"),
                                      labels = c("Clinical", "HLA", "Expression", "SNVs/Indels", "CNV/SVs"))

groups_to_labels$labels[order(match(groups_to_labels$names,sort(unique(ingredients$big_group))))]
groups_to_labels$color <- c("#FFFF99","#FDB462", "#FC913A", "#80B1D3", "B3DE69")

color_map_type$RNA = "#FC913A" # change the RNA feature value to a more direct value. 
color_map_type$`CNV/SVs` <- "#B3DE69"
color_map_type$Somatic <- "#80B1D3"
color_map_type$CNV <- color_map_type$SVs <- NULL
set.seed(12345) # set seed to keep random order
ingredients$pval_significance = ifelse(ingredients$p_val <= ingredients$by_05_fdr, "sign", "non-sign")
ingredients =  ingredients[sample(nrow(ingredients)),]

level = ingredients$big_group |> table() |> sort(decreasing = TRUE) |> names()
order = order(factor(ingredients$big_group, levels = level)) # make a factor for ordering the data
ingredients = ingredients[order, ]

ingredients_plot = ingredients %>%
    filter(Group != "Gene Set",
           dataset == "all",
           model == "bor" ,
           covariates == "age_biopsy_purity_tissue")

fig2a = ggplot(ingredients_plot, aes(x = plot_est, y = log10_p,
                           size = log10_p > threshold, 
                           alpha = log10_p > threshold)) + 
    geom_point(aes(fill = big_group), shape = 21, stroke = 0.3, color = "white") + 
    scale_y_continuous(breaks = c(seq(0,8,2)), limits = c(0,8)) + 
    scale_x_continuous(breaks = c(1, 2), limits = c(0.4,2.5))  +
    scale_size_manual(values = c(2, 3))  +
    scale_alpha_manual(values = c(0.4, 1)) + 
    exhaustive_theme_left_legend + 
    theme(legend.position = c(0.85, 0.25),   legend.background = element_blank(),
          legend.box.background = element_blank()) + 
    scale_fill_manual(values = c("#9E7788", "#B3DE69","#FFFF99", "#FC913A" ,"#80B1D3") ) + 
    geom_vline(xintercept = 1, linetype="dashed", color = "black", linewidth=.1) + 
    geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth=.1) + 
    guides(alpha = "none", size = "none", # remove the 'pval size'
           fill = guide_legend(title = "Group", override.aes = list(size=4)),
           color = guide_legend(title = "")) + 
  xlab("Response Odds Ratio Estimate") + ylab("-Log10 p-value") + 
  ggtitle("BOR vs All Features") 
    
fig2a

fig2a <- annotate_figure(fig2a, fig.lab = "A", fig.lab.size = 20, fig.lab.face = "bold")
ggsave("../plots/Fig2A.png", fig2a ,bg='black', width = 15, height = 13, units = "cm")

ready <- readRDS(paste0(I_DIR, "cpi_go.Rds"))
results <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))
cutoff<-results$by_05_fdr[1]

# filter clinical features to have only a single result
results = results %>% filter(feature != "clinical_systemic_composite")

one_set_result<- results %>% filter(model == "os", dataset == "all", covariates == "age_biopsy_purity_tissue") 
one_set_result<-unique(one_set_result) # filter for unique values

bor_set_result<- results %>% filter(model == "bor", dataset == "all", covariates == "age_biopsy_purity_tissue") 
bor_set_result<-unique(bor_set_result)

drop_features <- c("tcell", "tmb", "prolif","tgfb","pretreat","pretreat_comp","tissue",
                   "tissue_full","age","biopsy","purity","pdl1")

top_bor_df <- results %>% 
    filter(model == "bor", dataset == "all", covariates == "age_biopsy_purity_tissue") %>% 
    filter(!feature %in% drop_features) %>%
    filter(!grepl("gene_set_", feature)) %>%
    filter(!grepl("cibersort", feature)) %>%
    filter(p_val < cutoff)
top_bor = top_bor_df  %>%  pull(feature)

bor_cor_go <- ready %>% select(all_of(top_bor)) %>% drop_na() 
cor_bor_with_pval <- rcorr(as.matrix(bor_cor_go))
cor_matrix_bor = cor(as.matrix(bor_cor_go), use = "pairwise.complete.obs")
colnames(cor_matrix_bor) <- colnames(bor_cor_go)
rownames(cor_matrix_bor) <- colnames(bor_cor_go)

for_heatmap_within_genes_bor <- cor_matrix_bor
for_heatmap_within_genes_bor_with_group = merge(as.data.frame(for_heatmap_within_genes_bor),
                                                one_set_result[c("feature","Type","big_group")],
                                                by.y = "feature", by.x = 0)

# cluster the correlation values
set.seed(12)
ord <- hclust(dist(1-cor_matrix_bor))$order

# Cluster RNA values to be able to split RNA features into a T-cell matrix and two outliers
cor_rna = cor_matrix_bor %>% as.data.frame() %>%
    rownames_to_column("feature") %>%
    filter(grepl("isofox", feature)) %>% 
    select(contains("isofox"), feature) %>% 
    column_to_rownames("feature")
cor_tree = hclust(dist(1-abs(cor_rna)), method = "ward.D2")
# set names for the tree-clustering method:
clusters = cutree(cor_tree, 2)

wdf <- one_set_result[,c("feature", "big_group")] %>% filter(feature %in% top_bor)
wdf$w = ifelse(wdf$big_group == "Clinical", 40, ifelse(wdf$big_group == "RNA", 3,1))
wdf = wdf[ord,]
wdf$big_group_fct = factor(wdf$big_group, levels = c("RNA", "Clinical", "Somatic"))
wdf = wdf %>% 
    mutate(cluster = big_group) %>% 
    mutate(cluster = case_when(clusters[feature] == 1 ~ "cluster 1",
                              clusters[feature] == 2 ~ "cluster 2",
                              .default = cluster)) %>% 
    mutate(cluster = factor(cluster, levels = c("cluster 1", "cluster 2", "Clinical", "Somatic")))

wdf$clustord = 1:nrow(wdf)
wdf = wdf %>% arrange(cluster, clustord)
breaks <- wdf %>% 
  mutate(cumw = cumsum(w),
         pos = .5 * (cumw + lag(cumw, default = 0))) %>% 
  select(feature, pos)

bor_cor_go = bor_cor_go %>% select(all_of(wdf$feature))
pd <- as.data.frame(cor_matrix_bor) %>% select(all_of(wdf$feature))
pd.m = pd %>%
    rownames_to_column("feature") %>% 
    pivot_longer(cols = -feature)

# Join the the 'pd.m' dataframe and the 'breaks' dataframe on feature
pd.m <- pd.m %>%  
  left_join(breaks, by = "feature") %>% 
  rename(y = pos) %>% 
  left_join(breaks, by = c("name" = "feature")) %>% 
  rename(x = pos)

pd.m <- pd.m %>% 
  left_join(wdf, by = "feature")
pd.m <- pd.m %>% 
  left_join(wdf[,c(1,2,5)], by = c("name" = "feature"))

RNA_features = wdf[wdf$big_group == "RNA",]$feature
pd.m$feature <- factor( pd.m$feature, levels = colnames(bor_cor_go))
pd.m$name <- factor(pd.m$name, levels = colnames(bor_cor_go) )
pd.m$height = ifelse(pd.m$feature %in% c("clinical_pre_treated", "clinical_meta_hasSystemicPreTreatment2"), 40,
                     ifelse(pd.m$feature %in% RNA_features, 3,1))
pd.m$width = ifelse(pd.m$name == c("clinical_pre_treated", "clinical_meta_hasSystemicPreTreatment2"), 40,
                    ifelse(pd.m$name %in% RNA_features,3,1))

# Add the 'n value' to each column name in the plot
occurrences = table(wdf$cluster)
exp = paste0("Cluster 3:\nExpression \n(n = ", occurrences[1], ")")
snv_indel = paste0("Cluster 1: SNVs/Indels \n(n = ", occurrences[4], ")")

pd.m = pd.m %>% 
    mutate(Type = case_when(cluster.x == "cluster 1" ~ exp,
                            cluster.x == "cluster 2" ~ "cl2",
                            cluster.x == "Somatic" ~ snv_indel, 
                            cluster.x == "Clinical" ~ "Clinical")) %>% 
    mutate(Type2 = case_when(cluster.y == "cluster 1" ~ exp,
                            cluster.y == "cluster 2" ~ "cl2",
                            cluster.y == "Somatic" ~ snv_indel, 
                            cluster.y == "Clinical" ~ "Clinical"))

# change the factor levels to order the plot in the same way
pd.m = pd.m %>%
    mutate(Type2 = factor(Type2, levels = c(snv_indel,"Clinical",exp, "cl2")), 
           Type = factor(Type, levels = c(snv_indel,"Clinical",exp, "cl2")))

# Some aestetic changes to the plot
colnames(pd.m)[3] <- "correlation"
colour_breaks <- c(10, 20, 30)
colours <- c("darkblue", "white", "darkred")

# plot the heatmap
fig2b <- ggplot(pd.m, aes(x = x, y = y)) +
  geom_tile(aes(width = width,  height = height, fill = correlation)) +
  scale_x_reverse(breaks = breaks$pos, labels = breaks$feature, expand = c(0, 0.0)) +
  scale_y_continuous(breaks = breaks$pos, labels = breaks$feature, expand = c(0, 0.0)) +
  scale_fill_gradientn(limits = c(-1, 1), 
                       colours = c("blue", "white", "white", "red"), 
                       values = c(0, 0.43, 0.57, 1.01)) +
  theme_bw(base_size = 15) +
  facet_grid2(Type ~ Type2, scales = "free", space = "free", switch = "both",render_empty = TRUE,
               strip = strip_themed(background_x = elem_list_rect(fill = c("#80B1D3","#9E7788","#FB8072", "lightgrey"), size= 1),
                                    background_y = elem_list_rect(fill = c("#80B1D3","#9E7788", "#FB8072", "lightgrey"), size= 1))) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        panel.spacing.x = unit(0, "mm"),
        panel.spacing.y = unit(0, "mm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "black", linewidth = 0),
        panel.grid  = element_blank(),
        legend.position = "right") + 
        xlab("") + ylab("")
fig2b
fig2b <- annotate_figure(fig2b, fig.lab = "B", fig.lab.size = 20, fig.lab.face = "bold") 

ggsave("1_figures/revised_plots_Axel/plots/Fig2B.png", fig2b ,bg='white', width = 15, height = 13, units = "cm")

read_features = ready %>% select(any_of(unique(ingredients$feature)))

# get the features of each cluster
features_cl1 = wdf %>% filter(big_group == "Somatic") %>% pull(feature)
features_cl2 = wdf %>% filter(big_group == "Clinical") %>% pull(feature)
features_cl3 = wdf %>% filter(big_group == "RNA" & cluster == "cluster 1") %>% pull(feature)

mean_cl1 = rowMeans(ready %>% select(all_of(features_cl1)))
mean_cl2 = rowMeans(ready %>% select(all_of(features_cl2)))
mean_cl3 = rowMeans(ready %>% select(all_of(features_cl3)))

cluster_cors = data.table(feature = colnames(read_features),
                   `Cluster 1: SNVs/indels` = cor(mean_cl1, read_features, use = "pairwise.complete.obs") %>% as.numeric(),
                   `Cluster 2: Clinical` = cor(mean_cl2, read_features, use = "pairwise.complete.obs") %>% as.numeric(),
                   `Cluster 3: Expression` = cor(mean_cl3, read_features, use = "pairwise.complete.obs") %>% as.numeric())

i_cors = ingredients %>% filter(model == "bor" & dataset == "all" & covariates == "age_biopsy_purity_tissue")
cluster_cors = cluster_cors %>% filter(feature %in% unique(i_cors$feature))
i_cors = left_join(i_cors, cluster_cors,by = "feature" )

# merge and filter the labels to match value: 
label_cors = labels %>%
    filter(model == "bor" & dataset == "all" & covariates == "age_biopsy_purity_tissue")
label_cors = left_join(label_cors, cluster_cors,by = "feature" ) %>% 
    filter(log10_p > threshold) %>% 
    filter(!labels %in% c("rna", "data_driven"))

label_cors = label_cors %>% 
    filter(!clean_label %in% c("TGFB", "Proliferation", "Prior Systemic Therapy", "T-cell")) 

label_cors = label_cors %>% 
    select(starts_with("Cluster"), big_group, log10_p, clean_label) %>% 
    pivot_longer(cols = starts_with("Cluster"), names_to = "cluster", values_to = "correlation") %>% 
    filter(correlation > 0.75)

input_cors = i_cors %>%  
    filter(model == "bor" & dataset == "all" & covariates == "age_biopsy_purity_tissue") %>% 
    select(starts_with("Cluster"), big_group, log10_p, feature) %>% 
    pivot_longer(cols = starts_with("Cluster"), names_to = "cluster", values_to = "correlation")

input_cors = input_cors %>% 
    filter(big_group %in% c("Somatic", "CNV/SVs") & cluster == "Cluster 1: SNVs/indels"|
           big_group %in% c("Clinical", "HLA") & cluster == "Cluster 2: Clinical" | 
           big_group == "RNA" & cluster == "Cluster 3: Expression")

input_cors = input_cors  %>% 
    mutate(in_cluster = case_when(
        feature %in% c(features_cl1, features_cl2, features_cl3) ~ "in cluster", 
        .default = 'not in cluster'))

fig2e = ggplot(input_cors, aes(x = correlation, y = log10_p,
                           size = correlation > 0.5, 
                           alpha = correlation > 0.5)) + 
    geom_point(aes(fill = big_group, color = in_cluster), shape = 21, stroke = 0.3) + 
    facet_rep_grid(. ~ cluster )  + 
    geom_label_repel(data = label_cors, aes(label = clean_label),min.segment.length =  0.2, 
                      alpha = 0.85, size = 5,force =  5, box.padding = 2) + 
    exhaustive_theme_left_legend + 
    scale_y_continuous(breaks = c(seq(0,8,2)), limits = c(0,8)) +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1))  +
    scale_size_manual(values = c(1, 2)) +
    scale_alpha_manual(values = c(0.5, 1)) + 
    scale_color_manual(values = c("black", "white")) + 
    scale_fill_manual(values = c("#9E7788", "#B3DE69","#FFFF99", "#FC913A" ,"#80B1D3")) + 
    geom_vline(xintercept = 0.5, linetype="dashed", color = "black", linewidth=.1) + 
    geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth=.1) + 
    guides(alpha = "none", size = "none", # remove the 'pval size'
           fill = guide_legend(title = "Group", override.aes = list(size=4)),
           color = guide_legend(title = "")) +
   theme(strip.background = element_blank()) + 
   xlab("Correlation to cluster mean") + ylab("-Log10 p-value") 
fig2e

fig2e <- annotate_figure(fig2e, fig.lab = "E", fig.lab.size = 20, fig.lab.face = "bold") 
ggsave("../plots/figure2E.png", fig2e, width = 13, height = 5)

# load gene sets
gsea_sets <- readRDS(paste0(REF_DIR,"GSEA_gene_sets.Rds"))
cpi_gene_sets <- readRDS(paste0(REF_DIR,"cpi1000_gene_sets.Rds"))
gsea_sets = gsea_sets[setdiff(names(gsea_sets), names(cpi_gene_sets))]
gene_sets = c(gsea_sets, cpi_gene_sets) # combine the two geneset lists

# Read the interpretable geneset names
geneset_names = data.table::fread(paste0(REF_DIR, "CPI_genesets.txt"))
colnames(geneset_names)[1] = "feature"

# Select T-cell effector geneset as an example to show the correlation
select_sets = gene_sets[c("gene_set_t_cell_effector")]
select_sets = lapply(select_sets, \(x) paste0("isofox_", x)) # add isofox to identify genes
cluster_genes = colnames(pd)[grep("isofox", colnames(pd))] 
total_sets = c(select_sets, list(cluster_mean = cluster_genes))
total_sets = sapply(total_sets, \(x) x[(x %in% colnames(ready))]) # select the number of genes not in the data

# generate cluster means: 
ready_m = ready |> select(starts_with("isofox")) |> as.matrix()# make matrix to allow for looping
patient_means = sapply(total_sets, \(x) rowMeans(ready_m[,x]))
colnames(patient_means)[1] = geneset_names$Short_name[match(colnames(patient_means)[1], geneset_names$feature)]

pm = patient_means |> 
    as.data.frame() |>
    filter(!is.na(cluster_mean))

# For now, name the intermediate figure fig2bc
fig2c = ggplot(pm, aes(x = `T-Cell effector`, y = cluster_mean)) + 
  geom_point(color = "white", fill = "grey20", shape = 21, size = 3) + 
  exhaustive_theme_left+
  ggpubr::stat_cor(aes(label = after_stat(r.label)), size = 6) + 
  theme(strip.background = element_blank(), strip.text = element_text(face="bold")) + 
  ylab("Mean expression\nT-Cell effector geneset by patient") +
  xlab("Mean expression cluster 3 genes by patient") + 
  ggtitle("T-Cell effector")
    
fig2c <- annotate_figure(fig2c, fig.lab = "C", fig.lab.size = 20, fig.lab.face = "bold")
fig2c

ggsave("/workspace/users/arosendahl/plots/Figure2C.png", fig2c, height = 6, width = 4)

# get the significant RNA genes
top_bor_RNA = top_bor_df  %>% 
    filter(big_group == "RNA")  %>% 
    pull(feature)

top_genes <- ready %>% 
    select(all_of(top_bor_RNA)) %>% 
    select(!all_of(c("isofox_BAIAP2", "isofox_LATS2", "isofox_GGT5"))) 

top_genes$cluster_mean = rowMeans(top_genes)
cor_genes<-cor(top_genes, use="pairwise.complete.obs")

cluster_genes <- str_split_i(colnames(top_genes), "_",2)
# combine gsea and self-cluster genesets to a joint set and summarize the information in this gene set
gene_sets_info = data.frame(
    feature = sub("gene_set_", "", names(gene_sets)), 
    Length = lengths(gene_sets), 
    Intersection_cluster = sapply(gene_sets, \(x) sum(x %in% cluster_genes)),
    subset = c(rep("gsea_genesets", length(gsea_sets)),
               rep("cpi_genesets", length(cpi_gene_sets))))

gene_sets_info = gene_sets_info %>%   mutate(cluster_genes_proportion = Intersection_cluster/Length)
gene_sets_info = gene_sets_info %>% filter(!grepl('_rand|impres|vhio', feature))
gene_sets_info = gene_sets_info %>% filter(!grepl("character", feature))

gene_sets_expr <- ready %>% select(contains("gene_set_")) 
colnames(gene_sets_expr) <- sub("isofox_", "", colnames(gene_sets_expr))

top_genes_with_genesets = cbind(top_genes$cluster_mean, gene_sets_expr)
colnames(top_genes_with_genesets)[1]<-"Cluster_mean"

cor_to_cluster_mean <- data.frame("feature" = colnames(top_genes_with_genesets), 
                                  "Cor_with_cluster" = cor(top_genes_with_genesets,top_genes_with_genesets$Cluster_mean,use="pairwise.complete.obs")[,1])
cor_to_cluster_mean = cor_to_cluster_mean |> mutate(feature = sub("gene_set_", "", feature))

result <- merge(cor_to_cluster_mean, gene_sets_info, by="feature")
colnames(result)[1]<-"feature"
result$cluster_genes_proportion = result$Intersection/result$Length

# merge the feature names with the more interpretable names
geneset_names$feature = gsub("gene_set_", "", geneset_names$feature)
result = left_join(result, geneset_names, by = "feature")
colnames(result)

bor_set_result_genesets <- bor_set_result %>% filter(grepl("gene_set", feature))
bor_set_result_genesets$feature = str_split_i(bor_set_result_genesets$feature, "set_",2)

result2<- merge(result, bor_set_result_genesets, by="feature")
result2$selected = ifelse(result2$Cor_with_cluster > 0.8, "selected", "no")

fig2d <- (ggplot(result2, aes(x = Cor_with_cluster, y = log10_p, col = selected)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab("Correlation gene set with expression cluster") +
  ylab("-Log10 (p-value)") +
  #     scale_color_gradient2(low="grey", high="#FB8072")+
  scale_color_manual(values=c("darkgrey","#FB8072")) +
  geom_vline(xintercept = 0.8, linetype="dashed", color = "black", linewidth =.1) +
  geom_hline(yintercept = threshold, linetype="dashed", color = "black", linewidth =.2) +
  geom_text_repel(aes(label = if_else((Cor_with_cluster > 0.8 & log10_p > threshold), Short_name, "")), size=2.7, col="grey20" ,max.overlaps=20)+
  theme(legend.position="none") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size = 16))
  + scale_y_continuous(breaks = c(2,4,6,8), lim = c(0,8)) 
  + exhaustive_theme_left
)
fig2d <- annotate_figure(fig2d, fig.lab = "D", fig.lab.size = 20, fig.lab.face = "bold")
fig2d

ggsave("/workspace/users/arosendahl/plots/Figure2D.png", fig2d, height = 6, width = 6)

# combine the plots with the SVG image:
library(grImport2)
library(rsvg)
rsvg_svg("/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_3.svg", "/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_3_cairo.svg")
SVGlogo <- readPicture("/workspace/users/arosendahl/immune_biomarkers/ref/latent_factor_3_cairo.svg")
svg_pic = pictureGrob(SVGlogo, ext="gridSVG", delayContent=FALSE)
svg_pic = annotate_figure(svg_pic, fig.lab = "F", fig.lab.size = 20, fig.lab.face = "bold")

pl = list(fig2a, fig2b,  fig2c, fig2d, fig2e)
pl = lapply(pl, "+", theme(plot.margin=margin(3,3,3,3, unit = "mm")))
pl[[2]] = pl[[2]] + theme(plot.margin=margin(0,0,0,0, unit = "mm"))

# plot alternative plot (own idea): 
lay = rbind(c(1,1,1,2,2,2),
            c(3,3,3,3,3,3), 
            c(4,4,5,5,6,6))

lay = rbind(c(1,2), 
            c(3,4),
            c(5,5), 
            c(6,NA))

pl[[3]] = pl[[3]] + theme(plot.margin=margin(5,5,5,5, unit = "mm"))
fig <- arrangeGrob(pl[[1]], pl[[2]], pl[[3]], pl[[4]],pl[[5]], svg_pic, layout_matrix=lay)  
ggsave("1_figures/revised_plots_Axel/plots/Fig2.png", fig, width = 14, height = 20, dpi = 200)
ggsave("1_figures/revised_plots_Axel/plots/Fig2.pdf", fig, width = 14, height = 20, dpi = 200)
