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
    mouse = read_tsv("data/heart_quantNormData_mouse_1sd.tsv", progress = FALSE),
    rat = read_tsv("data/heart_quantNormData_rat_1sd.tsv", progress = FALSE)
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
ggsave("plots/heart_corDens.svg", figS2, device = svglite, width = 10, height = 4, units = "cm", scale = 1.8)

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
links0_75 <- map(corVal, ~getPairOfProbes(.x, threshold = 0.75))

# building file for SCHYPE -------------------
mappingTable <- read_tsv("data/mouse_rat_ampping.txt")

# sorting node-1 node-2 alphabeticaly
sortInteractions <- function(df) {
    map_df(seq_len(nrow(df)), function(x) {
        sorted_pair <- unlist(df[x, ]) %>% sort
        return(data_frame(node_1 = sorted_pair[1], node_2 = sorted_pair[2]))
    }) %>% distinct
}

# 40 minutes
t0 <- Sys.time()
sorted_links <- map(links, sortInteractions)
Sys.time() - t0

# 3 minutes
t0 <- Sys.time()
sorted_links0_75 <- map(links0_75, sortInteractions)
Sys.time() - t0

# converting rat probe_id into mouse
ratHomologousInterations <- left_join(sorted_links$rat, mappingTable, by = c("node_1" = "rat_ps")) %>%
    rename(node_1_mouse = mouse_ps) %>%
    left_join(mappingTable, by = c("node_2" = "rat_ps")) %>%
    rename(node_2_mouse = mouse_ps)

ratHomologousInterations0_75 <- left_join(sorted_links0_75$rat, mappingTable, by = c("node_1" = "rat_ps")) %>%
    rename(node_1_mouse = mouse_ps) %>%
    left_join(mappingTable, by = c("node_2" = "rat_ps")) %>%
    rename(node_2_mouse = mouse_ps)

# sorting intreactions
t0 <- Sys.time() # 4 minutes
ratHomologousInterations2 <- map_df(seq_len(nrow(ratHomologousInterations)), function(x) {
    line <- ratHomologousInterations[x, ]
    sorted_pair <- unlist(line)[3:4] %>% sort
    return(data_frame(node_1 = line$node_1, node_2 = line$node_2, node_1_mouse = sorted_pair[1], node_2_mouse = sorted_pair[2]))
}) %>% distinct
Sys.time() - t0

t0 <- Sys.time() # 1 minute
ratHomologousInterations2_0_75 <- map_df(seq_len(nrow(ratHomologousInterations0_75)), function(x) {
    line <- ratHomologousInterations0_75[x, ]
    sorted_pair <- unlist(line)[3:4] %>% sort
    return(data_frame(node_1 = line$node_1, node_2 = line$node_2, node_1_mouse = sorted_pair[1], node_2_mouse = sorted_pair[2]))
}) %>% distinct
Sys.time() - t0

filter(ratHomologousInterations2, node_1 != node_2)
filter(ratHomologousInterations2_0_75, node_1 != node_2)

homolog_links <- inner_join(
    mutate(sorted_links$mouse, key = paste(node_1, node_2)),
    mutate(ratHomologousInterations2, key = paste(node_1_mouse, node_2_mouse)),
    by = "key"
) %>% rename(
    node_1_m = node_1.x, node_2_m = node_2.x,
    node_1_r = node_1.y, node_2_r = node_2.y
)
homolog_links0_75 <- inner_join(
    mutate(sorted_links0_75$mouse, key = paste(node_1, node_2)),
    mutate(ratHomologousInterations2_0_75, key = paste(node_1_mouse, node_2_mouse)),
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
schype_input0_75 <- with(homolog_links0_75, data_frame(
    nodes = paste0(
        node_1_m, "_mouse ", node_2_m, "_mouse | ",
        node_1_r, "_rat ", node_2_r, "_rat"
    )
))

write.table(
    schype_input,
    file = "data/fileForSCHypeThreshold0.5_heart.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)
write.table(
    schype_input0_75,
    file = "data/fileForSCHypeThreshold0.75_heart.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

# running schype, instant
module load compilers/java/current
java -jar ../../programmes/schype-master/SCHype.jar -hgfile data/fileForSCHypeThreshold0.5_heart.txt -dir true -output data/schype_output_0.5th_heart
java -jar ../../programmes/schype-master/SCHype.jar -hgfile data/fileForSCHypeThreshold0.75_heart.txt -dir true -output data/schype_output_0.75th_heart

