setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)
library(cowplot)
library(tidyr)
library(annotationTools)
library(future); plan(multiprocess)

t0 <- Sys.time() # 25 sec
normData <- list(
    mouse = read_tsv("data/quantNormData_mouse_1sd.tsv", progress = FALSE),
    rat = read_tsv("data/quantNormData_rat_1sd.tsv", progress = FALSE)
)
Sys.time() - t0

normData_m <- map(
    normData,
    function(x) {
        y <- as.matrix(x[, c(-1, -2, -3)])
        rownames(y) <- x$ID
        return(y)
    }
)

t0 <- Sys.time() # 1 min
corVal <- map(normData_m, ~cor(t(.x), method = "spearman"))
Sys.time() - t0

# distribution of correlation values
corDistrib <- map(seq_along(corVal), function(x) {
    as_data_frame(corVal[[x]]) %>%
        gather %>%
        filter(value != 1) %>%
        ggplot(aes(x = value)) +
        geom_density() +
        geom_vline(xintercept = 0.5, color = "darkred") +
        labs(x = "Correlation (Spearman)", y = "Denisty", title = names(corVal)[x])
})

figS2 <- plot_grid(
    corDistrib[[1]], corDistrib[[2]],
    ncol = 2, labels = LETTERS[1:2], label_size = 22, hjust = -0.7, vjust = 1.2
)
ggsave("plots/figS2_corDens.svg", figS2, device = svglite, width = 10, height = 4, units = "cm", scale = 1.8)

# get pairs of probes ---------------
getPairOfProbes <- function(mat, threshold = 0.5) {
    bmat <- mat >= threshold
    bmat[!lower.tri(bmat)] <- FALSE
    res <- as_data_frame(bmat) %>%
        mutate(probe_id = rownames(mat)) %>%
        gather(probe_id2, th, -probe_id) %>%
        filter(th) %>%
        select(probe_id, probe_id2)
    return(res)
}

links <- map(corVal, getPairOfProbes)

# homology links from homologene---------------------
# wget ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data
# ^ 2016-02-27
# homoloGene <- read_tsv("data/homologene.data", col_names = FALSE) %>% as.data.frame # ps2ps doesn't like tibbles
# mouseID <- 10090
# ratID <- 10116
#
# # affy annotation form pia / affymetrix
# annotMouse <- read_csv("data/Mouse430_2.na36.annot.csv",  comment = "#") %>% as.data.frame
# annotRat <- read_csv("data/Rat230_2.na36.annot.csv", comment = "#") %>% as.data.frame
#
# # long... 15 minutes
# mappingTableMouse %<-% ps2ps(annotMouse, annotRat, homoloGene, ratID)
# mappingTableRat %<-% ps2ps(annotRat, annotMouse, homoloGene, mouseID)
#
# head(mappingTableMouse)
# head(mappingTableRat)
# dim(mappingTableMouse)
# dim(mappingTableRat)
# write.table(mappingTableMouse, file = "data/mappingTableMouse.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(mappingTableRat, file = "data/mappingTableRat.txt", sep = "\t", quote = FALSE, row.names = FALSE)

mappingTableMouse <- read_tsv("data/mappingTableMouse.txt")
mappingTableRat <- read_tsv("data/mappingTableRat.txt")
dim(mappingTableMouse)
dim(mappingTableRat)
dim(annotMouse)
dim(annotRat)

t0 <- Sys.time() # 1 min
mappingTableMouse2 <- map_df(seq_len(nrow(mappingTableMouse)), function(x) {
    x <- mappingTableMouse[x, ]
    return(data_frame(
        ps_1 = x$ps_1, gid_1 = x$gid_1, gid_2 = x$gid_2, ps_2 = unlist(strsplit(x$ps_2, ",", fixed = TRUE))
    ))
})
Sys.time() - t0
t0 <- Sys.time() # 1 min
mappingTableRat2 <- map_df(seq_len(nrow(mappingTableRat)), function(x) {
    x <- mappingTableRat[x, ]
    return(data_frame(
        ps_1 = x$ps_1, gid_1 = x$gid_1, gid_2 = x$gid_2, ps_2 = unlist(strsplit(x$ps_2, ",", fixed = TRUE))
    ))
})
Sys.time() - t0
mappingTableMouse2 <- select(mappingTableMouse2, ps_1, ps_2) %>% rename(mouse_ps = ps_1, rat_ps = ps_2)
mappingTableRat2 <- select(mappingTableRat2, ps_1, ps_2) %>% rename(rat_ps = ps_1, mouse_ps = ps_2)

mappingTable <- bind_rows(mappingTableMouse2, mappingTableRat2) %>% distinct

dim(mappingTableMouse2)
dim(mappingTableRat2)
dim(mappingTable)

write.table(
    mappingTable,
    file = "data/mouse_rat_ampping.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
)


# building file for SCHYPE -------------------

# sorting node-1 node-2 alphabeticaly
sortInteractions <- function(df) {
    map_df(seq_len(nrow(df)), function(x) {
        sorted_pair <- unlist(df[x, ]) %>% sort
        return(data_frame(node_1 = sorted_pair[1], node_2 = sorted_pair[2]))
    }) %>% distinct
}

# slow...
sorted_links <- map(links, sortInteractions)

# converting rat probe_id into mouse
ratHomologousInterations <- left_join(sorted_links$rat, mappingTable, by = c("node_1" = "rat_ps")) %>%
    rename(node_1_mouse = mouse_ps) %>%
    left_join(mappingTable, by = c("node_2" = "rat_ps")) %>%
    rename(node_2_mouse = mouse_ps)

# sorting intreactions
t0 <- Sys.time() # 20 minutes
ratHomologousInterations2 <- map_df(seq_len(nrow(ratHomologousInterations)), function(x) {
    line <- ratHomologousInterations[x, ]
    sorted_pair <- unlist(line)[3:4] %>% sort
    return(data_frame(node_1 = line$node_1, node_2 = line$node_2, node_1_mouse = sorted_pair[1], node_2_mouse = sorted_pair[2]))
}) %>% distinct
Sys.time() - t0

filter(ratHomologousInterations2, node_1 != node_2)

homolog_links <- inner_join(
    mutate(sorted_links$mouse, key = paste(node_1, node_2)),
    mutate(ratHomologousInterations2, key = paste(node_1_mouse, node_2_mouse)),
    by = "key"
) %>% rename(
    node_1_m = node_1.x, node_2_m = node_2.x,
    node_1_r = node_1.y, node_2_r = node_2.y
)

schype_input <- with(homolog_links, data_frame(
    nodes = paste0(
        node_1_m, "_mouse ", node_2_m, "_mouse | ",
        node_1_r, "_rat ", node_2_r, "_rat"
    )
))

write.table(
    schype_input,
    file = "data/fileForSCHypeThreshold0.5.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

# running schype, instant
java -jar ../../programmes/schype-master/SCHype.jar -hgfile data/fileForSCHypeThreshold0.5.txt -dir true -output data/schype_output_0.5th minclustsize 10




