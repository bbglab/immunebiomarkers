pos <- function(df, i){
    ests <- df %>% filter(feature == i, dataset %in% c("skin","lung","bladder","other")) %>% pull(est)
    data.frame("feature" = i, pos_dir = sum(ests > 0))
}
add_dirs <- function(df){
    dirs <- data.frame(); 
    features <- unique(df %>% pull(feature))
    for (i in features) {
        dirs <- rbind(dirs, pos(df, i))
    }
    df %>% inner_join(dirs, by = "feature")
}
organizer <- function( df ) {
    #df$feature_group <- unlist(lapply(df$feature, function(i) strsplit(i, "_")[[1]][1]))
    (df %>% mutate(log10_p = -log10(p_val))
        %>% relocate(group, dataset, model, model_type,
        col_type, feature, cor_pretreat, cor_tmb, cor_tcell, cor_prolif, cor_tgfb, cor_purity, est, se, p_val, log10_p, 
        gene_set_typ))
}
clean_stats <- function(df, cor, gene_map){
    df <- (add_dirs(df) 
               %>% left_join(cor, by="feature") 
               %>% left_join(gene_map, by="feature") 
               %>% drop_na(p_val) 
           )
    organizer(df)
}
clean_stats2 <- function(df, cor, gene_map){
    out <- tryCatch( clean_stats( df, cor, gene_map), error = function(e) NULL)
    return(out)
}
