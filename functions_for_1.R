library(cowplot)
library(tidyr)

# tbl <- origData$mouse

plotBoxplotsFromTibble <- function(tbl, which_col = 108:124, title = NULL) {
    myPlot <- select(tbl, c(1, which_col)) %>% #ID in 1st column
        gather(key = "geo_id", value = "expression_level", -ID) %>%
        ggplot(aes(x = geo_id, y = expression_level)) +
        geom_boxplot(outlier.shape = NA) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
        labs(x = NULL, y = "expression level", title = title)
    return(myPlot)
}

replaceInf <- function(x) {
    x[x == -Inf] <- min(x[x != -Inf])
    return(x)
}


