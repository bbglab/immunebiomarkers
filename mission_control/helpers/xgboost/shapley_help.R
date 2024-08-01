
get_Y2 <- function( df, model, complete = FALSE){
    if (complete) df <- df %>% drop_na()
    
    if (model == "lr") {
        as.matrix( df %>% select(response))
    }
    else if (model == "pfs") {
        as.matrix( df %>% select(pfs))
    }
    else if (model == "os") {
        as.matrix( df %>% select(os))
    } 
}
get_X2 <- function(df, model_features, complete = FALSE){
    if (complete) df <- df %>% drop_na()
    as.matrix( df %>% select(model_features) )
}

# Below - originally published at https://liuyanguu.github.io/post/2018/10/14/shap-visualization-for-xgboost/

std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}

shap_formatter <- function( shaps, X ){
    require(data.table)
    
    ### Put shaps into long format ###
    shaps <- as.data.table(shaps)[,BIAS:=NULL]
    mn_shaps <- colMeans(abs(shaps))[order(colMeans(abs(shaps)), decreasing = T)]
    shaps <- shaps[, names(mn_shaps), with = F]
    shaps_long <- melt.data.table(shaps, measure.vars = colnames(shaps))
    
    ### features ### 
    X <- as.data.table(X)[,names(mn_shaps), with = F]
    X_long <- melt.data.table(X, measure.vars = colnames(X))
    X_long[, stdfvalue := std1(value), by = "variable"]
    names(X_long) <- c("variable", "rfvalue", "stdfvalue" )    

    ### combine into one long table
    shaps_long2 <- cbind(shaps_long, X_long[,c('rfvalue','stdfvalue')])
    shaps_long2[, mean_value := mean(abs(value)), by = variable]
    setkey(shaps_long2, variable)
    list("nice_shaps" = shaps_long2, "mn_shaps" = mn_shaps)
}

var_importance <- function(mean_shap_score, top_n=12){
  var_importance=tibble(var=names(mean_shap_score), importance=mean_shap_score)
  var_importance=var_importance[1:top_n,]
  
  ggplot(var_importance, aes(x=reorder(var,importance), y=importance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_light() + theme(axis.title.y=element_blank()) 
}

plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value)) + .05
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    # sina plot: 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    # print the mean absolute value: 
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) + 
    scale_color_gradient(low="#FFCC33", high="#6600CC", breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position="bottom") + 
    geom_hline(yintercept = 0) + # the vertical line
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}

######### Combined plotting function
plot_store <- function( attic, model, type, tissue ){
    a <- var_importance(attic[[model]][[type]][[tissue]][['nice_shaps']]$mn_shaps)
    b <- plot.shap.summary(attic[[model]][[type]][[tissue]][['nice_shaps']]$nice_shaps %>% filter(mean_value != 0))
    
    a <- (a + ggtitle(paste(model, type, tissue, sep = "-")) + 
      theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = .5)))
    
    grid.arrange( a, b, nrow = 2)
}
