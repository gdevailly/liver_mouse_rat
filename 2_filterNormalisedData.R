setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(svglite)

source("scripts/functions_for_2.R")

# mean variance plot and 1SD selection ----------------
t0 <- Sys.time() # 25 sec
normData <- list(
    mouse = read_tsv("data/quantNormData_mouse.tsv", progress = FALSE),
    rat = read_tsv("data/quantNormData_rat.tsv", progress = FALSE)
)
Sys.time() - t0

t0 <- Sys.time() # 20 sec
mean_variance <- map(
    normData,
    function(x) {
        y <- as.matrix(x[, -1])
        mutate(x, mean = apply(y, 1, mean), sd = apply(y, 1, sd)) %>%
            select(ID, mean, sd, everything())
    }
)
Sys.time() - t0

mv_plots <- map(seq_along(mean_variance), function(x) plotMeanVariance(mean_variance[[x]], title = names(mean_variance)[x]))
mv_plots <- plot_grid(mv_plots[[1]], mv_plots[[2]], ncol = 2, labels = c("A", "B"))
ggsave("plots/fig1_meanvariance.png", width = 18, height = 10, units = "cm", scale = 0.8)

variableProbes <- map(
    mean_variance,
    ~filter(.x, sd >= 1)
)

# saving data
write.table(
    variableProbes$mouse, file = "data/quantNormData_mouse_1sd.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
write.table(
    variableProbes$rat, file = "data/quantNormData_rat_1sd.tsv",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)


# files for GO analysis -------------------

anno <- list(
    mouse = read_csv("data/Mouse430_2.na36.annot.csv",  comment = "#"),
    rat = read_csv("data/Rat230_2.na36.annot.csv", comment = "#")
)
anno <- map(anno, function(x) {
    colnames(x)[1] <- "ID"
    return(x)
})

getGeneList <- function(probes_id, annotations, what = "Ensembl") {
    dplyr::filter(annotations, ID %in% probes_id) %>%
        dplyr::select_(what) %>%
        unlist %>%
        strsplit(" ", fixed = TRUE) %>%
        map_chr(first) %>%
        unname() %>%
        grep("ENS", ., value = TRUE)
}

t0 <- Sys.time()
all_gene <- list(
    mouse = getGeneList(mean_variance$mouse$ID,  anno$mouse) ,
    rat = getGeneList(mean_variance$rat$ID,  anno$rat)
)
Sys.time() - t0

t0 <- Sys.time()
variable_gene <- list(
    mouse = getGeneList(variableProbes$mouse$ID,  anno$mouse) ,
    rat = getGeneList(variableProbes$rat$ID,  anno$rat)
)
Sys.time() - t0

write.table(all_gene$mouse, file = "data/gl_all_gene_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(all_gene$rat, file = "data/gl_all_gene_rat.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(variable_gene$mouse, file = "data/gl_variable_gene_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(variable_gene$rat, file = "data/gl_variable_gene_rat.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

