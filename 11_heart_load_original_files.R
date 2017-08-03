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
    mouse = read_tsv("data/mouse_heart_data", progress = FALSE, skip = 1), # missing values in ctrl probes, to delete last 2 rows
    # lots of duplicate column, are values the same ?
    rat = read_tsv("data/rat_heart_expression", progress = FALSE, skip = 1)
)
Sys.time() - t0

map_df(origData, dim)
origData$mouse <- dplyr::slice(origData$mouse, 1:45099)
map_df(origData, dim)

dplyr::select(origData$mouse, GSM147451, GSM147451_1)
dplyr::select(origData$mouse, GSM147452, GSM147452_1, GSM147452_2)
grep("_", colnames(origData$mouse), value = TRUE, fixed = TRUE) # ID_REF !
grep("_", colnames(origData$rat), value = TRUE, fixed = TRUE) # no dup in rat

map_df(origData, dim)
origData$mouse <- dplyr::select(origData$mouse, -grep("_", colnames(origData$mouse), value = FALSE, fixed = TRUE)[-1])
map_df(origData, dim)
origData <- map(origData, function(x) {
    colnames(x)[1] <- "ID"
    return(x)
})

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

# quantile normalisation with liver mouse as a target

liver_mouse <- read_tsv("data/quantNormData_mouse.tsv", progress = FALSE)
liver_mouse_m <- as.matrix(liver_mouse[, -1])
rownames(liver_mouse_m) <- liver_mouse$ID

quartNormData_m <- map(quartNormData, function(x) {
    y <- as.matrix(x[, -1])
    rownames(y) <- x$ID
    return(y)
})

t0 <- Sys.time() #
targetMouse <- normalize.quantiles.determine.target(liver_mouse_m, target.length = nrow(quartNormData_m$mouse))
targetRat <- normalize.quantiles.determine.target(liver_mouse_m, target.length = nrow(quartNormData_m$rat))
quantNormData_m <- list(
    mouse = normalize.quantiles.use.target(quartNormData_m$mouse, targetMouse),
    rat = normalize.quantiles.use.target(quartNormData_m$rat, targetRat)
)
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
    bp_before_norm[[1]] + coord_cartesian(ylim = c(0, 500)), bp_before_norm[[2]] + coord_cartesian(ylim = c(0, 500)),
    bp_quart_norm[[1]] + coord_cartesian(ylim = c(0, 500)), bp_quart_norm[[2]] + coord_cartesian(ylim = c(0, 500)),
    bp_quant_norm[[1]], bp_quant_norm[[2]],
    ncol = 2, labels = LETTERS[1:6], label_size = 22, hjust = -0.7, vjust = 1.2
)
ggsave("plots/heart_normalisations.svg", figS1, device = svglite, width = 10, height = 12, units = "cm", scale = 1.8)

# saving data
write.table(
    quantNormData$mouse, file = "data/heart_quantNormData_mouse.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
write.table(
    quantNormData$rat, file = "data/heart_quantNormData_rat.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
