setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(future); plan(multiprocess(workers = 1 + 32))

load("data/hclust_data_heart.RData")

anno <- list(
    mouse = read_csv("data/Mouse430_2.na36.annot.csv",  comment = "#", progress = FALSE),
    rat = read_csv("data/Rat230_2.na36.annot.csv", comment = "#", progress = FALSE)
)
anno <- map(anno, function(x) {
    colnames(x)[1] <- "ID"
    return(x)
})

geneClusters <- function(probes_id, annotations, what = "Ensembl") {
    dplyr::filter(annotations, ID %in% probes_id) %>%
        dplyr::select_(what) %>%
        unlist %>%
        strsplit(" ", fixed = TRUE) %>%
        map_chr(first) %>%
        unname() %>%
        grep("ENS", ., value = TRUE)
}

myGeneClusters <- map(rowClusters, function(x) {
    map(unique(x), function(y) {
        names(x[x == y])
    })
})

map(myGeneClusters, ~map_int(.x, length))

ENSgIDs <- map2(myGeneClusters, anno, function(x, y) {
    map(x, function(z) geneClusters(z, y))
})

map(ENSgIDs, ~map_int(.x, length))

walk(seq_along(ENSgIDs), function(i) {
    map(seq_along(ENSgIDs[[i]]), function(j) {
        write.table(
            ENSgIDs[[i]][[j]],
            file = paste0("data/heart_gl_hclust_", names(ENSgIDs)[i], "_cluster", j, "_", length(ENSgIDs[[i]][[j]]), "genes.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE
        )
    })
})


