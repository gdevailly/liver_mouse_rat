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

seriesInfo <- list(
    mouse = read_tsv("data/gse_gsm_mouse.tsv", progress = FALSE),
    rat = read_tsv("data/gse_gsm_rat.tsv", progress = FALSE)
)

map_df(
    seriesInfo,
    ~c("gse" = n_distinct(.x$gse), "gsm" = n_distinct(.x$gsm))
)
# mouse   rat
# <int> <int>
# 20    19
# 248  1202

seriesTable <- map(seriesInfo, ~table(.x$gse) %>%
                       as_data_frame %>%
                       arrange(desc(n)))

write.table(
    seriesTable$mouse,
    file = "supp_table_2_mouse_heart.txt",
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
)
write.table(
    seriesTable$rat,
    file = "supp_table_2_rat_hear.txt",
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
)

t0 <- Sys.time() # 3 sec
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

log2FC_m <- map(
    normData_m,
    function(x) {
        t(apply(x, 1, function(l) l - mean(l)))
    }
)


t0 <- Sys.time() # 23 minutes
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
colDend2 <- map2(colDend2, c(4, 5), ~color_branches(.x, k = .y))
Sys.time() - t0

rowClusters <- map2(rowDend, c(7, 7), ~cutree(.x, k = .y))
colClusters <- map2(colDend, c(4, 5), ~cutree(.x, k = .y))

rowColors <- map(rowDend2, get_leaves_branches_col)
colColors <- map(colDend2, get_leaves_branches_col)

save(rowDend, colDend, rowClusters, colClusters, rowColors, colColors, log2FC_m, file = "data/hclust_data_heart.RData")

# mouse GEO accession color generation, a wee bit of work --------------
# b & w decoration--------
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
png(file = "plots/testHeatmap_mouse_dend_geoAcc3_heart.png", width = 15, heigh = 15, units = "cm", res = 600)
draw(mh)
dev.off()
Sys.time() - t0


#rat GEO accession color generation, a wee bit of work --------------
# b & w decoration--------
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
png(file = "plots/testHeatmap_rat_dend_geoAcc3_heart.png", width = 15, heigh = 15, units = "cm", res = 600)
draw(mh)
dev.off()
Sys.time() - t0



