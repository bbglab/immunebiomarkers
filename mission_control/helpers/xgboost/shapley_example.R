patient_pred_plot <- function( df, type = "lr" ) {
    
    if( type == "lr" ) {
        base <- (
            ggplot(df, aes(y=patient_id, x=pred, fill = col, alpha = more)) +
            scale_x_continuous( breaks = c(.1,.5,1), limits = c(0,1)) + 
            ggtitle( paste0( "Probability of Response (", round(df$pred,2), ")") ) + 
            scale_fill_manual(values = c('Low' = my_palette[1], 'Medium' = my_palette[2],'High' = my_palette[3])) +
            scale_alpha_manual(values = c(1)) 
            #+geom_text(size = 7, vjust = -2)
        )
    } else {
        base <- (
            ggplot(df, aes(y=patient_id, x=pred, fill = col2, alpha = more)) +
            scale_x_continuous( breaks = c(.5,1.5), limits = c(0,6)) + 
            ggtitle( paste0( "Predicted Hazard (", round(df$pred,1), ")")  ) +
            scale_fill_manual(values = c('Low' = my_palette[3], 'Medium' = my_palette[2],'High' = my_palette[1])) +
            scale_alpha_manual(values = c(.1))
        )
    }
    gg <- ( base + 
            geom_bar(position="dodge", stat="identity", color = "black") + 
            theme_patient
    )
    gg

}

patient_shap_plot <- function( df, type = "lr") {
    gg <- (
        ggplot(df, aes(y=model, fill = feature, group = feature, x=shap_feature)) 
            + geom_bar(position="dodge", stat="identity", color = "black")
            + scale_fill_manual(values=unlist(color_map))
            + theme_shaps
            + scale_x_continuous(limits = c(-1.8,1.8))
    )
    
    if( type == "lr"){
        gg + ggtitle("BOR Shapley Values" ) + labs(x = "Log Odds of Response")
    } else {
        gg + ggtitle("OS Shapley Values" ) + labs(x = "Log Hazard")
    }
    
}

patient_plots <- function( df, id ){
    
    tmp <- df %>% filter(patient_id == id)
    print(tmp %>% select(dataset, age, gender, bor, os) %>% head(1))
    
    a <- patient_pred_plot( tmp %>% filter(model == "lr"), "lr")
    b <- patient_shap_plot( tmp %>% filter(model == "lr"), "lr")
    c <- patient_pred_plot( tmp %>% filter(model == "os"), "os")
    d <- patient_shap_plot( tmp %>% filter(model == "os"), "os")
    
    lay <- rbind(c(1,2),c(3,4),c(3,4),c(3,4))
    arrangeGrob(a,c,b,d, layout_matrix = lay)

}
