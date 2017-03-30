setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)
library(cowplot)
library(tidyr)
library(future); plan(multiprocess)
library(viridis)


mappingTable <- read_tsv("data/mouse_rat_ampping.txt")

t0 <- Sys.time() # 25 sec
normData <- list(
    mouse = read_tsv("data/quantNormData_mouse_1sd.tsv", progress = FALSE),
    rat = read_tsv("data/quantNormData_rat_1sd.tsv", progress = FALSE)
)
Sys.time() - t0


# variable probe homology -----------------------------------
variableProbes <- map(normData, ~.x$ID)

# few seconds
withHomolog <- list(
    mouse = map_lgl(variableProbes$mouse, function(x) {
        homolog <- filter(mappingTable, mouse_ps == x)$rat_ps
        return(any(homolog %in% variableProbes$rat))
    }),
    rat = map_lgl(variableProbes$rat, function(x) {
        homolog <- filter(mappingTable, rat_ps == x)$mouse_ps
        return(any(homolog %in% variableProbes$mouse))
    })
)

map(withHomolog, summary)

# hclust cluster homology -----------------------------------

load("data/hclust_data.RData")

geneClusters <- function(probes_id, annotations, what = "Ensembl") {
    dplyr::filter(annotations, ID %in% probes_id) %>%
        dplyr::select_(what) %>%
        unlist %>%
        strsplit(" ", fixed = TRUE) %>%
        map_chr(first) %>%
        unname() %>%
        grep("ENS", ., value = TRUE)
}

myGeneClusters <- map2(rowClusters, rowColors, function(x, y) {
    names(y) <- names(x)
    map(unique(y), function(z) {
        names(y[y == z])
    })
})

map(myGeneClusters, ~map_int(.x, length))

getJaccardIndexOfHomolog <- function(mouseGeneList, ratGeneList, annotationTable = mappingTable) {
    mouse_b <- map_lgl(mouseGeneList, function(x) {
        homolog <- filter(annotationTable, mouse_ps == x)$rat_ps
        return(any(homolog %in% ratGeneList))
    })
    rat_b <- map_lgl(ratGeneList, function(x) {
        homolog <- filter(annotationTable, rat_ps == x)$mouse_ps
        return(any(homolog %in% mouseGeneList))
    })
    inner <- mean(sum(mouse_b), sum(rat_b)) # average of the two. They can be different because of one2many homology...
    jaccard <- inner / (length(mouse_b) + length(rat_b) - inner)
    pval <- 1 - phyper(
        inner,
        length(mouse_b),
        (length(unique(annotationTable$mouse_ps)) + length(unique(annotationTable$rat_ps))) / 2 -length(mouse_b), # average of probe number in annotation table
        length(rat_b)
    )
    return(list(jaccard = jaccard, pval = pval))
}

getJaccardIndexOfHomolog(myGeneClusters$mouse[[1]], myGeneClusters$rat[[1]])

jaccard_df <- data_frame(
    cluster_mouse = rep(1:7, each = 7),
    cluster_rat = rep(1:7, times = 7)
)

t0 <- Sys.time() # 2 minutes
jaccard_output <- map_df(
    seq_len(nrow(jaccard_df)),
    function(i) getJaccardIndexOfHomolog(
        myGeneClusters$mouse[[jaccard_df$cluster_mouse[i]]],
        myGeneClusters$rat[[jaccard_df$cluster_rat[i]]]
    )
)
Sys.time() - t0

jaccard_df <- bind_cols(jaccard_df, jaccard_output)

clusterColors <-  colorspace::rainbow_hcl(7, c=90, l=50)
image(matrix(1:7), col = clusterColors)

jaccard_df <- mutate(jaccard_df,
    apval = p.adjust(pval, method = "bonferroni")
) %>% mutate(.,
    mstars = case_when(
        .$apval < 0.0001 ~ "***",
        .$apval < 0.001 ~ "**",
        .$apval < 0.01 ~ "*",
        .$apval < 0.05 ~ "+",
        TRUE ~ ""
    )
)

map_int(myGeneClusters$rat, length)
xlabs <- vector()
xlabs[1] <- expression(atop(bold(1), scriptstyle("(13)"  )))
xlabs[2] <- expression(atop(bold(2), scriptstyle("(533)" )))
xlabs[3] <- expression(atop(bold(3), scriptstyle("(1086)")))
xlabs[4] <- expression(atop(bold(4), scriptstyle("(419)" )))
xlabs[5] <- expression(atop(bold(5), scriptstyle("(31)"  )))
xlabs[6] <- expression(atop(bold(6), scriptstyle("(32)"  )))
xlabs[7] <- expression(atop(bold(7), scriptstyle("(2)"   )))

map_int(myGeneClusters$mouse, length)
ylabs <- vector()
ylabs[1] <- expression(paste(scriptstyle("(16)"  ), bold(1)))
ylabs[2] <- expression(paste(scriptstyle("(1032)"), bold(2)))
ylabs[3] <- expression(paste(scriptstyle("(111)" ), bold(3)))
ylabs[4] <- expression(paste(scriptstyle("(22)"  ), bold(4)))
ylabs[5] <- expression(paste(scriptstyle("(984)" ), bold(5)))
ylabs[6] <- expression(paste(scriptstyle("(1606)"), bold(6)))
ylabs[7] <- expression(paste(scriptstyle("(6)"   ), bold(7)))

myPlot <- ggplot(jaccard_df, aes(x = factor(cluster_rat), y = factor(cluster_mouse), fill = jaccard, label = mstars)) +
    geom_tile() +
    scale_fill_viridis(name = "Jaccard\nindex") +
    geom_text(size = 3, color = "firebrick1") +
    labs(x = "Rat clusters", y = "Mouse clusters", title = "Homology relationships\nbetween clusters") +
    scale_x_discrete(labels = xlabs) +
    scale_y_discrete(labels = ylabs) +
    theme(
        axis.text.x = element_text(colour = clusterColors),
        axis.text.y = element_text(colour = clusterColors)
    )

ggsave("plots/fig2X_hclust_homology.svg", myPlot, device = svglite, width = 6, height = 4, units = "cm", scale = 1.8)

map(rowColors, ~length(which(.x == "#A86B00")))


