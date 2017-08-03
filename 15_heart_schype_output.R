setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)
library(cowplot)
library(tidyr)
library(future); plan(multiprocess)
library(ComplexHeatmap)
library(circlize)
library(cba)
library(dendextend)

# thanks so http://stackoverflow.com/a/8197703/5463724
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

mappingTable <- read_tsv("data/mouse_rat_ampping.txt")
schypeOutput <- read_tsv("data/schype_output_0.5th_heart.nodes.txt", col_names = c("probes", "cluster_id"))

schSpe <- list(
    mouse = filter(schypeOutput, grepl("_mouse", probes, fixed = TRUE)),
    rat = filter(schypeOutput, grepl("_rat", probes, fixed = TRUE))
)

map_int(schSpe, nrow)
# 974   757
map_int(schSpe, ~length(unique(.x$probes)))
#   454   307

seriesInfo <- list(
    mouse = read_tsv("data/gse_gsm_mouse.tsv", progress = FALSE),
    rat = read_tsv("data/gse_gsm_rat.tsv", progress = FALSE)
)

# cluster histogram ----------------
ngenePerCluster <- map_df(
    unique(schypeOutput$cluster_id),
    function(x) {
        data_frame(
            cluster_id = x,
            n_mouse = nrow(filter(schSpe$mouse, cluster_id == x)),
            n_rat = nrow(filter(schSpe$rat, cluster_id == x))
        )
    }
) %>% arrange(desc(n_mouse)) %>%
    filter(n_mouse >= 10 | n_rat >= 10) %>%
    mutate(new_cluster_id = seq_len(length(cluster_id))) %>%
    gather("species", "n_gene", -cluster_id, -new_cluster_id)

mp <- ggplot(ngenePerCluster, aes(x = factor(new_cluster_id, levels = 30:1), y = n_gene, fill = species)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(x = "cluster", y = "Number of probes", title = "SCHype clustering") +
    scale_fill_manual(
        values = c("grey25", "grey75"),
        name = "Species",
        breaks = c("n_mouse", "n_rat"),
        labels = c("Mouse", "Rat")
    ) +
    coord_flip() +
    theme(legend.position = "top")

ggsave("plots/heart_SCHypeClusterSize.svg", mp, device = svglite, width = 4, height = 8, units = "cm", scale = 1.8)


# schype heatmaps ------------------
clusterToKeep <- c("c1" = 1, "c2" = 0, "c3" = 12, "c4" = 4)

listOfGene <- map(clusterToKeep, function(x) {
    map(schSpe, function(y) {
        dplyr::filter(y, cluster_id == x)$probes %>%
            gsub("_mouse|_rat", "", .)
    })
})

load("data/hclust_data_heart.RData")

listOfGene2 <- map(c(mouse = "mouse", rat = "rat"), function(x) {
    list(
        c1 = listOfGene$c1[[x]],
        c2 = listOfGene$c2[[x]],
        c3 = listOfGene$c3[[x]],
        c4 = listOfGene$c4[[x]]
    )
})

schye_cluster <- map2(log2FC_m, listOfGene2, function(x, y) {
    map(y, function(z) {
        x[z, ]
    })
})

schye_cluster_fused <- map(schye_cluster, ~do.call(rbind, .x))
schye_cluster_indices <- map(schye_cluster, function(x) {
    c(
        rep("c1", nrow(x$c1)),
        rep("c2", nrow(x$c2)),
        rep("c3", nrow(x$c3)),
        rep("c4", nrow(x$c4))
    )
})

hclust_original_clustering <- map2(schye_cluster_fused, rowClusters, function(x, y) {
    y[rownames(x)]
})

ggpalette <-  colorspace::rainbow_hcl(7, c=90, l=50)

colDend <- map(schye_cluster_fused, function(x) {
    mdis <- dist(t(x))
    mhcl <- hclust(mdis)
    moop <- order.optimal(mdis, mhcl$merge)
    mhcl$merge  <- moop$merge
    mhcl$order   <- moop$order
    return(mhcl)
})
colDend2 <- map(colDend, ~set(as.dendrogram(.x), "branches_lwd", 0.5))

# mouse heatmap -------------
table(hclust_original_clustering[[1]])
table(rowClusters[[1]])
image(matrix(1:7), col = colorspace::rainbow_hcl(7, c=90, l=50))

hclust_anno_mouse <- rowAnnotation(
    df = data.frame("hclust" = paste0("c", hclust_original_clustering[[1]])),
    name = "hclust",
    col = list(hclust = c("c1" =  ggpalette[1], "c2" = ggpalette[6], "c3" = ggpalette[3], "c4" = ggpalette[2], "c5" = ggpalette[4])),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_name_gp = list("cex" = 0.8)
)

# b& w anno ---------
accNames <- seriesInfo$mouse[match(labels(colDend2$mouse), seriesInfo$mouse$gsm),]$gse %>%
    unique
names(accNames) <- accNames
bg <- "grey95"
fg <- "grey15"
cols <- map(accNames, function(x) {
    mcols <- rep(bg, length(accNames))
    names(mcols) <- accNames
    mcols[x] <- fg

    if(bg == "grey95") bg <<- "grey85" else bg <<- "grey95"
    if(fg == "grey15") fg <<- "black" else fg <<- "grey15"

    return(mcols)

})
df <- as.data.frame(matrix(
    seriesInfo$mouse$gse,
    ncol = length(cols),
    nrow = length(seriesInfo$mouse$gse)
))
names(df) <- names(cols)

decGeoAcc <- columnAnnotation(
    df = df,
    name = "GEO series",
    col = cols,
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_height = unit(rep(0.05, length(cols)), "cm")
    #gap =  unit(rep(0.02, length(cols)), "cm")
)
decGeoAcc_mouse <- decGeoAcc

# heatmap --------------------
t0 <- Sys.time() # 20 s
mh <- Heatmap(
    matrix = schye_cluster_fused[[1]],
    col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
    name = "Expression\nlevel",
    split = schye_cluster_indices[[1]],
    show_row_dend = FALSE,
    cluster_columns = colDend2$mouse,
    column_dend_height = unit(2, "cm"),
    row_title = "Probes",
    column_title = "Experiments",
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = decGeoAcc,
    bottom_annotation_height = unit(3, "cm"),
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 4
)
png(file = "plots/Heatmap_mouse_dend_SChype3_heart.png", width = 15, heigh = 12, units = "cm", res = 600)
hclust_anno_mouse + mh
dev.off()
Sys.time() - t0

# rat heatmap -------------
table(hclust_original_clustering[[2]])
table(rowClusters[[2]])
image(matrix(1:7), col = colorspace::rainbow_hcl(7, c=90, l=50))


hclust_anno_rat <- rowAnnotation(
    df = data.frame("hclust" = paste0("c", hclust_original_clustering[[2]])),
    name = "hclust",
    col = list(hclust = c("c2" =  ggpalette[3], "c3" = ggpalette[6], "c5" = ggpalette[7])),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_name_gp = list("cex" = 0.8)
    # annotation_height = unit(rep(1.4/3, 3), "cm")
)

accNames <- seriesInfo$rat[match(labels(colDend2$rat), seriesInfo$rat$gsm),]$gse %>%
    unique
names(accNames) <- accNames
bg <- "grey95"
fg <- "grey15"
cols <- map(accNames, function(x) {
    mcols <- rep(bg, length(accNames))
    names(mcols) <- accNames
    mcols[x] <- fg

    if(bg == "grey95") bg <<- "grey85" else bg <<- "grey95"
    if(fg == "grey15") fg <<- "black" else fg <<- "grey15"

    return(mcols)
})
df <- as.data.frame(matrix(
    seriesInfo$rat$gse,
    ncol = length(cols),
    nrow = length(seriesInfo$rat$gse)
))
names(df) <- names(cols)

decGeoAcc <- columnAnnotation(
    df = df,
    name = "GEO series",
    col = cols,
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_height = unit(rep(0.05, length(cols)), "cm")
    #gap =  unit(rep(0.02, length(cols)), "cm")
)
decGeoAcc_rat <- decGeoAcc

# heatmap rat ---------------
t0 <- Sys.time() # 1 minutes
mh <- Heatmap(
    matrix = schye_cluster_fused[[2]],
    col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
    name = "Expression\nlevel",
    split = schye_cluster_indices[[2]],
    show_row_dend = FALSE,
    cluster_columns = colDend2$rat,
    column_dend_height = unit(2, "cm"),
    row_title = "Probes",
    column_title = "Experiments",
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = decGeoAcc,
    bottom_annotation_height = unit(3, "cm"),
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 4
)
png(file = "plots/Heatmap_rat_dend_SChype3_heart.png", width = 15, heigh = 12, units = "cm", res = 600)
hclust_anno_rat + mh
dev.off()
Sys.time() - t0

# gene lists for go analysis ---------------
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

map(listOfGene2, ~map_int(.x, length))
# $mouse
# c1 c2 c3 c4
# 50 48 29 28
#
# $rat
# c1 c2 c3 c4
# 30 22 16 21

ENSgIDs <- map2(listOfGene2, anno, function(x, y) {
    map(x, function(z) geneClusters(z, y))
})

map(ENSgIDs, ~map_int(.x, length))
# $mouse
# c1 c2 c3 c4
# 49 46 25 26
#
# $rat
# c1 c2 c3 c4
# 26 14 13 17

map(seq_along(ENSgIDs), function(i) {
    map(seq_along(ENSgIDs[[i]]), function(j) {
        write.table(
            ENSgIDs[[i]][[j]],
            file = paste0("data/gl_SCHype_clust_heart_", names(ENSgIDs)[i], "_cluster", j, "_", length(ENSgIDs[[i]][[j]]), "genes.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE
        )
    })
}) %>% invisible

# cluster by cluster analysis ----------------------

hclust_original_clustering2 <- map2(schye_cluster, rowClusters, function(x, y) {
    map(x, function(z) {
        y[rownames(z)]
    })
})

col_clusters_2 <- map(schye_cluster, function(species) {
    map(species, function(myClusters) {
        mdis <- dist(t(myClusters))
        mhcl <- hclust(mdis)
        moop <- order.optimal(mdis, mhcl$merge)
        mhcl$merge  <- moop$merge
        mhcl$order   <- moop$order
        mhcl <- set(as.dendrogram(mhcl), "branches_lwd", 0.5)
        mhcl <- color_branches(mhcl, k = 2)
        return(mhcl)
    })
})


col_clusters_3 <- map(schye_cluster, function(species) {
    map(species, function(myCluster) {
        mdis <- dist(t(myCluster))
        mhcl <- hclust(mdis)
        moop <- order.optimal(mdis, mhcl$merge)
        mhcl$merge  <- moop$merge
        mhcl$order   <- moop$order
        mhcl <- set(as.dendrogram(mhcl), "branches_lwd", 0.5)
        mhcl <- color_branches(mhcl, k = 3)
        return(mhcl)
    })
})

geo_bottom_dec <- list(mouse = decGeoAcc_mouse, rat = decGeoAcc_rat)

hclust_colors <- list(
    mouse = list(hclust = c("c1" =  ggpalette[1], "c2" = ggpalette[6], "c3" = ggpalette[3], "c4" = ggpalette[2], "c5" = ggpalette[4])),
    rat = list(hclust = c("c2" =  ggpalette[3], "c3" = ggpalette[6], "c5" = ggpalette[7]))
)

myHms2 <- pmap(
    list(schye_cluster, col_clusters_2, geo_bottom_dec, hclust_original_clustering2, hclust_colors),
    function(species, dendro, decoration, hclust, colors) {
        map(1:4, function(i) {
            hclust_anno <- rowAnnotation(
                df = data.frame("hclust" = paste0("c", hclust[[i]])),
                name = "hclust",
                col = colors,
                show_legend = FALSE,
                show_annotation_name = FALSE,
                annotation_name_gp = list("cex" = 0.8)
            )
            mh <- Heatmap(
                matrix = species[[i]],
                col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
                name = "log(FC)",
                show_row_dend = FALSE,
                cluster_columns = dendro[[i]],
                column_dend_height = unit(1, "cm"),
                row_title = "Probes",
                column_title = "Experiments",
                show_row_names = FALSE,
                show_column_names = FALSE,
                bottom_annotation = decoration,
                bottom_annotation_height = unit(3, "cm"),
                use_raster = TRUE,
                raster_device = "CairoPNG",
                raster_quality = 4
            )
            return(hclust_anno + mh)
        })
    }
)

myHms3 <- pmap(
    list(schye_cluster, col_clusters_3, geo_bottom_dec, hclust_original_clustering2, hclust_colors),
    function(species, dendro, decoration, hclust, colors) {
        map(1:4, function(i) {
            hclust_anno <- rowAnnotation(
                df = data.frame("hclust" = paste0("c", hclust[[i]])),
                name = "hclust",
                col = colors,
                show_legend = FALSE,
                show_annotation_name = FALSE,
                annotation_name_gp = list("cex" = 0.8)
            )
            mh <- Heatmap(
                matrix = species[[i]],
                col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
                name = "log(FC)",
                show_row_dend = FALSE,
                cluster_columns = dendro[[i]],
                column_dend_height = unit(1, "cm"),
                row_title = "Probes",
                column_title = "Experiments",
                show_row_names = FALSE,
                show_column_names = FALSE,
                bottom_annotation = decoration,
                bottom_annotation_height = unit(3, "cm"),
                use_raster = TRUE,
                raster_device = "CairoPNG",
                raster_quality = 4
            )
            return(hclust_anno + mh)
        })
    }
)

png(file = "plots/schypeHeatmaps/hm2_heart_%03d.png", width = 10, heigh = 7, units = "cm", res = 600)
myHms2
dev.off()
png(file = "plots/schypeHeatmaps/hm3_heart_%03d.png", width = 10, heigh = 7, units = "cm", res = 600)
myHms3
dev.off()
