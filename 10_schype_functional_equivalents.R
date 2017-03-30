setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)
library(cowplot)
library(tidyr)
library(future); plan(multiprocess)

mappingTable <- read_tsv("data/mouse_rat_ampping.txt")

anno <- list(
    mouse = read_csv("data/Mouse430_2.na36.annot.csv",  comment = "#"),
    rat = read_csv("data/Rat230_2.na36.annot.csv", comment = "#")
)
anno <- map(anno, function(x) {
    colnames(x)[1] <- "ID"
    return(x)
})

homoloGene <- read_tsv("data/homologene.data", col_names = FALSE)
mouseID <- 10090
ratID <- 10116
colnames(homoloGene) <- c("homology_group", "tax_ID", "gene_ID", "gene_symbol", "protein_ID", "protein_acc")

homoloGene <- filter(homoloGene, tax_ID %in% c(mouseID, ratID))

# identification of homology groups with more than 2 members ----------------------------
table(
    (homoloGene %>% group_by(homology_group) %>% summarise(n = n()))$n
)

# 1     2     3     4     5     6     7     8     9    10    11    12    13
# 1489 17367   704   189    67    51    34    13    11    10     4     7    10
# 14    15    16    17    18    19    20    22    23    27    29    30    36
# 5     2     5     3     2     1     2     3     3     2     1     1     1
# 38    44    46    68    80   110   296
# 1     1     1     1     1     1     1

nGene_per_group <- homoloGene %>%
    group_by(homology_group) %>%
    summarise(n = n())

big_groups <- filter(nGene_per_group, n >= 3)

multiH <- filter(homoloGene, homology_group %in% big_groups$homology_group)
table(
    (multiH %>%  group_by(homology_group) %>% summarise(species = n_distinct(tax_ID)))$species
)
# 1   2
# 250 888

n_species <- multiH %>%
    group_by(homology_group) %>%
    summarise(species = n_distinct(tax_ID)) %>%
    filter(species == 2)

multiH <- filter(multiH, homology_group %in% n_species$homology_group)

# SCHype output filtering --------------------------------

schypeOutput <- read_tsv("data/schype_output_0.5th.nodes.txt", col_names = c("probes", "cluster_id"))
schSpe <- list(
    mouse = filter(schypeOutput, grepl("_mouse", probes, fixed = TRUE)),
    rat = filter(schypeOutput, grepl("_rat", probes, fixed = TRUE))
)
schSpe <- map(schSpe, function(x) {
    x$probes <- gsub("_mouse|_rat", "", x$probes)
    return(x)
})

getGeneID <- function(probes_id, annotations, what = "Entrez Gene") {
    dplyr::filter(annotations, ID %in% probes_id) %>%
        .[, what] %>%
        unlist %>%
        strsplit(" /// ", fixed = TRUE) %>%
        unlist %>%
        unname %>%
        unique
}
getGeneID(schSpe$mouse$probes, anno$mouse) %>% head


schSpeCluster <- map(unique(schypeOutput$cluster_id), function(x) {
    myCluster <- filter(schypeOutput, cluster_id == x)
    mir <- list(
        mouse = filter(myCluster, grepl("_mouse", probes, fixed = TRUE)) %>%
            mutate(probes = gsub("_mouse", "", probes)),
        rat = filter(myCluster, grepl("_rat", probes, fixed = TRUE)) %>%
            mutate(probes = gsub("_rat", "", probes))
    )
    mr <- list(
        mouse = getGeneID(mir$mouse$probes, anno$mouse),
        rat = getGeneID(mir$rat$probes, anno$rat)
    )
})

clusterSize <- map_int(schSpeCluster, function(x) {
    min(length(x$mouse), length(x$rat))
})
# some cluster of 2 probes are for a single entrez gene
schSpeCluster <- schSpeCluster[clusterSize >= 2]

multiHlist <- map(unique(multiH$homology_group), ~filter(multiH, homology_group == .x))
names(multiHlist) <- paste0("h_", unique(multiH$homology_group))

getSchYpeClusterFor <- function(multiHgroup, schypeC = schSpeCluster) {
    mouse_genes <- filter(multiHgroup, tax_ID == mouseID)$gene_ID %>% as.character
    rat_genes <- filter(multiHgroup, tax_ID == ratID)$gene_ID %>% as.character
    return(
        map_lgl(schypeC, function(x) {
            any(mouse_genes %in% x$mouse) & any(rat_genes %in% x$rat)
        }) %>% which
    )
}

myOverlaps <- map(multiHlist, getSchYpeClusterFor)
myOverlaps <- myOverlaps[map_lgl(myOverlaps, ~length(.x) >= 1)]
length(myOverlaps) # 18

multiHlist$h_3938
schSpeCluster[[1]]

getMagicTable <- function(myOverlap, Overlaps = myOverlaps, homologs = multiHlist, clusters = schSpeCluster) {
    ot <- Overlaps[[myOverlap]]
    myCluster <- clusters[ot]
    myHomolog <- homologs[[myOverlap]]
    myHomolog$species <- ifelse(myHomolog$tax_ID == mouseID, "mouse", "rat")
    newColumns <- map(myCluster, function(x) {
        map_lgl(seq_len(nrow(myHomolog)), function(i) {
            myHomolog[i, "gene_ID"] %in% x[[unlist(myHomolog[i, "species"])]]
        })
    })
    names(newColumns) <- paste0("SCHype_cluster_", ot)
    newColumns <- do.call(data_frame, newColumns)
    return(bind_cols(myHomolog, newColumns))
}

names(myOverlaps)
getMagicTable("h_3938") %>% glimpse
getMagicTable("h_10699") %>% glimpse

ortholog_resolution <- map(names(myOverlaps), getMagicTable)

invisible(map(ortholog_resolution, ~write.table(
    .x,
    file = paste0("result/ortholog_resolution_", unlist(.x[1, "gene_symbol"]), ".txt"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, eol = "\r\n"
)))

