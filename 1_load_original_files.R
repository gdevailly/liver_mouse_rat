setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(preprocessCore)
library(purrr)
library(svglite)

source("scripts/functions_for_1.R")

# loading and cleaning data -----------------

t0 <- Sys.time() # 30s
origData <- list(
    mouse = read_tsv("data/matrixLiverMouse.txt", progress = FALSE),
    rat = read_tsv("data/rat_full_data_matrix.txt", progress = FALSE)
)
Sys.time() - t0

# plot of data before normalisation
bp_before_norm <- lapply(
    seq_along(origData),
    function(x) plotBoxplotsFromTibble(origData[[x]], title = paste0("Before normalisation (", names(origData)[x], ")"))
)

# replacing infinite values
cleanData <- map(origData, ~mutate_if(.x, is.numeric, replaceInf))

# lower quartile trimming
t0 <- Sys.time() # 20 sec
quartNormData <- map(cleanData, ~mutate_if(.x, is.numeric,
                                             function(x) {
                                                 myThreshold <- quantile(x, 0.25)
                                                 return(
                                                     ifelse(x < myThreshold, myThreshold, x)
                                                 )
                                             }))
Sys.time() - t0

bp_quart_norm <- lapply(
    seq_along(cleanData),
    function(x) plotBoxplotsFromTibble(quartNormData[[x]], title = paste0("Lower quartile trimming (", names(quartNormData)[x], ")"))
)

# quantile normalisation with mouse as a target
quartNormData_m <- map(quartNormData, function(x) {
    y <- as.matrix(x[, -1])
    rownames(y) <- x$ID
    return(y)
})

t0 <- Sys.time() # 41 sec
quantNormData_m <- list(mouse = normalize.quantiles(quartNormData_m$mouse))
target <- normalize.quantiles.determine.target(quantNormData_m$mouse, target.length = nrow(quartNormData_m$rat))
quantNormData_m$rat <- normalize.quantiles.use.target(quartNormData_m$rat, target)
quantNormData_m <- map2(quantNormData_m, quartNormData_m, function(x, y) {
    colnames(x) <- colnames(y)
    rownames(x) <- rownames(y)
    return(x)
})
Sys.time() - t0

quantNormData <- map(quantNormData_m, function(x) {
    y <- data.frame(ID = rownames(x), x) %>% as_data_frame
})

bp_quant_norm <- lapply(
    seq_along(quantNormData),
    function(x) plotBoxplotsFromTibble(quantNormData[[x]], title = paste0("Quantile normalisation (", names(quantNormData)[x], ")"))
)

# multipanel plot
figS1 <- plot_grid(
    bp_before_norm[[1]], bp_before_norm[[2]],
    bp_quart_norm[[1]], bp_quart_norm[[2]],
    bp_quant_norm[[1]], bp_quant_norm[[2]],
    ncol = 2, labels = LETTERS[1:6], label_size = 22, hjust = -0.7, vjust = 1.2
)
ggsave("plots/figS1_normalisations.svg", figS1, device = svglite, width = 10, height = 12, units = "cm", scale = 1.8)
# system("firefox plots/figS1_normalisations.svg &")

# saving data
write.table(
    quantNormData$mouse, file = "data/quantNormData_mouse.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
write.table(
    quantNormData$rat, file = "data/quantNormData_rat.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# experiment / series ----------------
seriesInfo <- list(
    mouse = read_tsv("data/GEOMouse.txt", progress = FALSE),
    rat = read_tsv("data/GEORat.txt", progress = FALSE)
)

# seriesInfoFiltered <- list(
#     mouse = filter(seriesInfo$mouse, Experiment %in% colnames(origData$mouse)),
#     rat = filter(seriesInfo$rat, Experiment %in% colnames(origData$rat))
# )

map(seriesInfo, ~length(unique(.x$GEO_Accession)))



