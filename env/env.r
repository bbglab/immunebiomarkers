#!/usr/bin/env Rscript

cat('Installing R dependencies\n')

install.packages('https://cran.r-project.org/src/contrib/Archive/ggborderline/ggborderline_0.1.0.tar.gz', repos=NULL, type='source')
install.packages('https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.0.1.tar.gz', repos=NULL, type='source')
install.packages('https://bioconductor.org/packages/3.16/data/experiment/src/contrib/geneLenDataBase_1.34.0.tar.gz', repos=NULL, type='source')
install.packages('https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.6.tar.gz', repos=NULL, type='source')