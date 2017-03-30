setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library(future); plan(multiprocess)
library(viridis)
library(dendextend)
library(cba)

# source("scripts/functions_for_2.R")

# mean variance plot and 1SD selection ----------------

seriesInfo <- list(
    mouse = read_tsv("data/GEOMouse.txt", progress = FALSE),
    rat = read_tsv("data/GEORat.txt", progress = FALSE)
)

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

log2FC_m <- map(
    normData_m,
    function(x) {
        t(apply(x, 1, function(l) l - mean(l)))
    }
)

t0 <- Sys.time() # 3 minutes
# rowDend <- map(log2FC_m, ~hclust(dist(.x)))
# colDend <- map(log2FC_m, ~hclust(dist(t(.x))))
rowDend <- map(log2FC_m, function(x) {
    mdis <- dist(x)
    mhcl <- hclust(mdis)
    moop <- order.optimal(mdis, mhcl$merge)
    mhcl$merge  <- moop$merge
    mhcl$order   <- moop$order
    return(mhcl)
})
colDend <- map(log2FC_m, function(x) {
    mdis <- dist(t(x))
    mhcl <- hclust(mdis)
    moop <- order.optimal(mdis, mhcl$merge)
    mhcl$merge  <- moop$merge
    mhcl$order   <- moop$order
    return(mhcl)
})
rowDend2 <- map(rowDend, ~set(as.dendrogram(.x), "branches_lwd", 0.5))
colDend2 <- map(colDend, ~set(as.dendrogram(.x), "branches_lwd", 0.5))
rowDend2 <- map2(rowDend2, c(7, 7), ~color_branches(.x, k = .y))
colDend2 <- map2(colDend2, c(4, 4), ~color_branches(.x, k = .y))
Sys.time() - t0

rowClusters <- map2(rowDend, c(7, 7), ~cutree(.x, k = .y))
colClusters <- map2(colDend, c(4, 4), ~cutree(.x, k = .y))

rowColors <- map(rowDend2, get_leaves_branches_col)
colColors <- map(colDend2, get_leaves_branches_col)

save(rowDend, colDend, rowClusters, colClusters, rowColors, colColors, log2FC_m, file = "data/hclust_data.RData")


# mouse GEO accession color generation, a wee bit of work --------------

# thanks so http://stackoverflow.com/a/8197703/5463724
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

accNames <- data_frame(
    geo_id = seriesInfo$mouse$GEO_Accession[colDend$mouse$order],
    rank = seq_len(length(geo_id))
) %>% group_by(geo_id) %>% summarise(rank = median(rank)) %>%
    arrange(rank)
accNames <- accNames$geo_id

cols <- list()
template <- rep("#EEEEEE", length(unique(seriesInfo$mouse$GEO_Accession)))
i <- 1
j <- 1
maxi <- length(unique(seriesInfo$mouse$GEO_Accession))

while(i <= maxi) {
    cols[[j]] <- template
    cols[[j]][i:min(c(i + 9, maxi))] <- gg_color_hue(length(cols[[j]][i:min(c(i + 9, maxi))]))
    i <- i + 10
    j <- j +1
}
cols <- map(cols, function(x) {
    names(x) <- accNames
    return(x)
})
names(cols) <- paste0("GEO", c("1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-62"))

df <- as.data.frame(matrix(
    seriesInfo$mouse$GEO_Accession,
    ncol = length(cols),
    nrow = length(seriesInfo$mouse$GEO_Accession)
))
names(df) <- names(cols)

decGeoAcc <- columnAnnotation(
    df = df,
    col = cols,
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = list("cex" = 0.5),
    annotation_height = unit(rep(0.2, 7), "cm")
)
# b & w decoration--------
accNames <- seriesInfo$mouse[match(labels(colDend2$mouse), seriesInfo$mouse$Experiment),]$GEO_Accession %>%
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
    seriesInfo$mouse$GEO_Accession,
    ncol = length(cols),
    nrow = length(seriesInfo$mouse$GEO_Accession)
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

t0 <- Sys.time() # 2 minutes
mh <- Heatmap(
    matrix = log2FC_m[[1]],
    col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
    name = "log(FC)",
    cluster_rows = rowDend2$mouse,
    cluster_columns = colDend2$mouse,
    bottom_annotation = decGeoAcc,
    bottom_annotation_height = unit(3, "cm"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_title = "Probes",
    column_title = "Experiments",
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 4
)
png(file = "plots/testHeatmap_mouse_dend_geoAcc3.png", width = 15, heigh = 15, units = "cm", res = 600)
draw(mh)
dev.off()
Sys.time() - t0

#rat GEO accession color generation, a wee bit of work --------------

accNames <- data_frame(
    geo_id = seriesInfo$rat$GEO_Accession[colDend$rat$order],
    rank = seq_len(length(geo_id))
) %>% group_by(geo_id) %>% summarise(rank = median(rank)) %>%
    arrange(rank)
accNames <- accNames$geo_id

cols <- list()
template <- rep("#EEEEEE", length(unique(seriesInfo$rat$GEO_Accession)))
i <- 1
j <- 1
maxi <- length(unique(seriesInfo$rat$GEO_Accession))

while(i <= maxi) {
    cols[[j]] <- template
    cols[[j]][i:min(c(i + 9, maxi))] <- gg_color_hue(length(cols[[j]][i:min(c(i + 9, maxi))]))
    i <- i + 10
    j <- j +1
}
cols <- map(cols, function(x) {
    names(x) <- accNames
    return(x)
})
names(cols) <- paste0("GEO", c("1-10", "11-20", "21-28"))

df <- as.data.frame(matrix(
    seriesInfo$rat$GEO_Accession,
    ncol = length(cols),
    nrow = length(seriesInfo$rat$GEO_Accession)
))
names(df) <- names(cols)

decGeoAcc <- columnAnnotation(
    df = df,
    col = cols,
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = list("cex" = 0.5),
    annotation_height = unit(rep(1.4/3, 3), "cm")
)

t0 <- Sys.time() # 1 minutes
mh <- Heatmap(
    matrix = log2FC_m[[2]],
    col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
    name = "Expression\nlevel",
    cluster_rows = rowDend2$rat,
    cluster_columns = colDend2$rat,
    bottom_annotation = decGeoAcc,
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_title = "Probes",
    column_title = "Experiments",
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 4
)
png(file = "plots/testHeatmap_rat_dend_geoAcc.png", width = 15, heigh = 15, units = "cm", res = 600)
draw(mh)
dev.off()
Sys.time() - t0


# b & w decoration--------
accNames <- seriesInfo$rat[match(labels(colDend2$rat), seriesInfo$rat$Experiment),]$GEO_Accession %>%
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
    seriesInfo$rat$GEO_Accession,
    ncol = length(cols),
    nrow = length(seriesInfo$rat$GEO_Accession)
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

t0 <- Sys.time() # 2 minutes
mh <- Heatmap(
    matrix = log2FC_m[[2]],
    col = colorRamp2(c(-3, 0, 3), colors = c("green", "black", "magenta")),
    name = "log(FC)",
    cluster_rows = rowDend2$rat,
    cluster_columns = colDend2$rat,
    bottom_annotation = decGeoAcc,
    bottom_annotation_height = unit(3, "cm"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_title = "Probes",
    column_title = "Experiments",
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 4
)
png(file = "plots/testHeatmap_rat_dend_geoAcc3.png", width = 15, heigh = 15, units = "cm", res = 600)
draw(mh)
dev.off()
Sys.time() - t0




