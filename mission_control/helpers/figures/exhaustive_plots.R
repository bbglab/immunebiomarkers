### Filter Function 
## Data and Labels 

exhaustive_plots_filter <- function(gold, i, j, x_axis = NULL, subset = NULL, k = 0.7){

    gold <- gold %>% filter(model == i, dataset == j)

    if (j == "all_adj" & i == "bor") {
        gold <- gold %>% filter(plot_est > 0.4)
    } else if (j == "all_adj" & i == "os") {
        gold <- gold %>% filter(plot_est > 0.6)
    } else if (j == "tmb_tcell_clin_adj") {
        gold <- gold %>% filter( cor_tmb < k, cor_tcell < k, cor_pretreat < k, Type != "Somatic")
    } else if (j == "six_adj") {
        gold <- gold %>% filter(cor_tmb < k, cor_tcell < k, cor_pretreat < k, cor_prolif < k, cor_tgfb < k, Type != "Somatic")
    }

    if (subset == "somatic") {
        gold <- gold %>% filter( grepl("Somatic", Type))
    } else if (subset == "rna") {
        gold <- gold %>% filter(grepl("RNA", Type))
    }

    if(x_axis == "cor_prolif") {
        gold <- gold %>% filter(Group != "Cibersort")
    } else if( x_axis == "cor_tcell" ){
        gold$Type_plus_main <- ifelse(gold$feature %in% c("somatic_TMB", "isofox_gene_set_t_cell_effector", "clinical_meta_hasSystemicPreTreatment2"), "Main", as.character(gold$Type))
    }
    gold
}

annotation_filter <- function(annotate, i, j, x_axis, subset = NULL, simple = "True", k = 0.5) {

    annotate <- annotate %>% filter(model == i, dataset == j)

    if (simple == "True") {
        annotate <- annotate %>% filter(size >= 4)
    }

    if (j == "tmb_tcell_clin_adj") {
        annotate <- annotate %>% filter(cor_tmb < k, cor_tcell < k, cor_pretreat < k, Type != "Somatic")
    } else if (j == "six_adj") {
        annotate <- annotate %>% filter(cor_tmb < k, cor_tcell < k, cor_pretreat < k, cor_prolif < k, cor_tgfb < k, Type != "Somatic")
    }

    if (x_axis == "cor_tcell") {
        annotate <- annotate %>% filter(grepl("T-cell", Type) | size > 3.5 | grepl("cibersort", feature))
    } else if (x_axis == "cor_prolif") {
        annotate <- annotate %>% filter(grepl("TGFB", Type) | grepl("Prolif", Type))
    }

    if (subset == "overall") {
        annotate <- annotate %>% filter(grepl("overall", labels) | labels == "data_driven")
    } else if (subset == "rna") {
        annotate <- annotate %>% filter(grepl("rna", labels))
    } else if (subset == "somatic") {
        annotate <- annotate %>% filter(grepl("tmb", labels))
    } else if (subset == "somatic_no_sig") {
        annotate <- (annotate %>% filter(!grepl("_SBS",feature), grepl("tmb",labels) | labels == "data_driven_somatic"))
    }

    annotate <- (annotate %>% select(-labels) %>% distinct())
    annotate
}

### Plots 
## Main and Correlation 

main_plot <- function( gold, i, j, x_axis, y_axis, subset = NULL, k = .7, simple = "True" ){
    
    gold <- exhaustive_plots_filter(gold, i, j, x_axis, subset, k)
    
    aesthetic <- aes( x = .data[[x_axis]], 
                      y = .data[[y_axis]], 
                      color = Main, 
                      fill = Type,
                      alpha = Type_plus_main, 
                      size  = Type_plus_main)
    gg <- (
        ggplot( gold, aesthetic)
            + geom_point(shape = 21)
            + scale_color_manual(values = c("black", "white"))
            + scale_fill_manual(values = unlist(color_map))
            + scale_alpha_manual(values = unlist(alpha_map)) 
            + scale_size_manual(values = unlist(size_map))   

            + labs( x = x_labeller(i, "plot_est"), y = "-Log10 (p-value)" ) 
            #+ geom_vline(xintercept = 1, linetype="dashed", color = "black", size=.2) 
            #+ geom_hline(yintercept = -log10(gold$by_05_fdr)[1], linetype="dashed", color = "black", size=.3) 
        )
    gg
}

cor_plot <- function( gold, i, j, x_axis, y_axis, subset = NULL, k = .7, simple = "True" ){
    
    gold <- exhaustive_plots_filter(gold, i, j, x_axis, subset, k)
    
    aesthetic <- aes( x = .data[[x_axis]], 
                      y = .data[[y_axis]], 
                      color = Type_plus_main,
                      fill = Type, 
                      alpha = Type_plus_main, 
                      size  = Type_plus_main, 
                      shape = Direction)
    gg <- (  ggplot( gold, aesthetic) + geom_point()

            + scale_color_manual(values=unlist(outline_map))
            + scale_fill_manual(values=unlist(color_map))
            + scale_alpha_manual(values = unlist(alpha_map)) 
            + scale_size_manual(values = unlist(size_map)) 
            + scale_shape_manual(values = c(24,25))
            
            + labs( x = x_labeller(i, x_axis), y = "-Log10 (p-value)")
            # + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=.3)
            # + geom_hline(yintercept = -log10(gold$by_05_fdr)[1], linetype="dashed", color = "black", size=.3)
        )
    gg 
    
}

### Add annotation  
add_annotation <- function (annotate, i, j, x_axis, y_axis, subset, simple = "True", k = 0.5, nudge_x = NULL, nudge_y = NULL, repel = FALSE, type = "cor", push = 1, pull = 1, seed = 622) {
    
    annotate <- annotation_filter(annotate, i, j, x_axis, subset, simple, k)
    
    aesthetic <- aes(x = .data[[x_axis]], y = .data[[y_axis]], label = clean_label)
    
    repel_labels <- (
        geom_text_repel(
                data = annotate, 
                aesthetic, 
                size = annotate %>% pull(size), 
                col = "black",
                point.padding = 1e-06, 
                alpha = 1, 
                force = push, 
                force_pull = pull, 
                nudge_y = nudge_y, 
                nudge_x = nudge_x,
                max.overlaps = 100,
                min.segment.length = unit(0, 'lines'), 
                inherit.aes = x_translate(x_axis),
                seed = seed
        )
    )
    normal_labels <- (
        geom_text(
                data = annotate, 
                aes(x = .data[[x_axis]], y = .data[[y_axis]], label = clean_label), 
                size = annotate %>% pull(size), 
                col = "black",
                alpha = 1, 
                nudge_y = nudge_y, 
                nudge_x = nudge_x,
                inherit.aes = x_translate(x_axis)
        )
    )
    
    if(repel){ 
	repel_labels
    } else {
        normal_labels
    }
}

super_plot <- function(ingredients, labels, type = "main", i, j, x_axis, y_axis, subset, simple = "True", k = 0.5, nudge_x = NULL, nudge_y = NULL, repel = FALSE,  push = 1, pull = 1, seed = 622){

    if( type == "main") {
        gg <- main_plot( ingredients, i, j, x_axis, y_axis, subset, k, simple )
    } else if (type == "cor"){
        gg <- cor_plot( ingredients, i, j, x_axis, y_axis, subset, k, simple )
    }
   annotate <- add_annotation( labels, i, j, x_axis, y_axis, subset, simple, k, nudge_x, nudge_y, repel, push = 1, pull = 1, seed = seed ) 
   gg + annotate

}

drop_features <- c("tcell", "tmb", "prolif","tgfb","pretreat","pretreat_comp","tissue",
                   "tissue_full","age","biopsy","purity","pdl1", "clinical_systemic_composite" )
