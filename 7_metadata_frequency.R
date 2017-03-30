setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(future); plan(multiprocess(workers = 1 + 32))
library(cowplot)
library(svglite)
library(tm)

load("data/hclust_data.RData")

geoMeta <- list(
    mouse = read_tsv("data/mouseCharacteristics.txt"),
    rat = read_tsv("data/ratCharacteristics.txt")
)

getWordCount <- function(myvect) {
    mr <- paste(myvect, collapse = " ") %>%
        VectorSource %>%
        Corpus %>%
        tm_map(content_transformer(tolower)) %>%
        tm_map(removePunctuation) %>%
        tm_map(stripWhitespace) %>%
        TermDocumentMatrix %>%
        as.matrix
    return(data_frame(term = rownames(mr), word_count = mr[, 1]))
}

expPerhclustCluster <- map(colClusters, function(x) {
    map(unique(x), function(y) {
        names(x[which(x == y)])
    })
})

myCharcthclust <- map2(expPerhclustCluster, geoMeta, function(x, y) {
    map(x, function(z) {
        filter(y, gsm_id %in% z)$characteristics
    })
})


# mouse frequency plot ------------------
myWordCounthclustMouse <- map(myCharcthclust$mouse, getWordCount)

myWordCounthclustMouseMerged <- full_join(myWordCounthclustMouse[[1]], myWordCounthclustMouse[[2]], by = "term") %>%
    full_join(myWordCounthclustMouse[[3]], by = "term") %>%
    full_join(myWordCounthclustMouse[[4]], by = "term")
colnames(myWordCounthclustMouseMerged) <- c("term", "c1", "c2", "c3", "c4")

myWordCounthclustMouseMerged$c1[is.na(myWordCounthclustMouseMerged$c1)] <- 0
myWordCounthclustMouseMerged$c2[is.na(myWordCounthclustMouseMerged$c2)] <- 0
myWordCounthclustMouseMerged$c3[is.na(myWordCounthclustMouseMerged$c3)] <- 0
myWordCounthclustMouseMerged$c4[is.na(myWordCounthclustMouseMerged$c4)] <- 0

map(expPerhclustCluster$mouse, length)
myWordCounthclustMouseMerged <- mutate(
    myWordCounthclustMouseMerged,
    c1f = c1/length(expPerhclustCluster$mouse[[1]]),
    c2f = c2/length(expPerhclustCluster$mouse[[2]]),
    c3f = c3/length(expPerhclustCluster$mouse[[3]]),
    c4f = c4/length(expPerhclustCluster$mouse[[4]])
)

myWordCounthclustMouseMerged <- mutate(
    myWordCounthclustMouseMerged,
    diff_c1_c3 = c1f - c3f
) %>% arrange(diff_c1_c3)

myWordCounthclustMouseMergedClean <- filter(
    myWordCounthclustMouseMerged,
    !(term %in% c("diet", "strain", "tissue", "gender", "liver", "for",
                  "dose", "mgkg", "route", "genotype", "months", "age",
                  "mixed", "high", "wks", "weeks", "hrs", "compound",
                  "sex", "8wk", "days", "treatment", "background", "strainbackground",
                  "status", "disease", "the", "genotypevariation", "genetic", "experiment",
                  "start", "mice", "mus", "musculus", "with"))
)
plot(myWordCounthclustMouseMergedClean$diff_c1_c3)

ggpalette <-  colorspace::rainbow_hcl(4, c=90, l=50)

p_mouse <- bind_rows(list(
    filter(
        myWordCounthclustMouseMergedClean,
        term %in% c("male", "female")
    ),
    head(myWordCounthclustMouseMergedClean, 5),
    tail(myWordCounthclustMouseMergedClean, 5)
)) %>%
    distinct %>%
    select(term, c1f, c3f) %>%
    gather("cluster", "frequency", -term) %>%
    mutate(term = factor(term, levels = rev(unique(c("female", "male", unique(term)))))) %>%
    ggplot(aes(x = term, y = frequency, fill = cluster)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.6) +
    labs(x = "", y = "Frequency", title = "Mouse") +
    scale_fill_manual(
        values = c(ggpalette[1], ggpalette[3]),
        name = "",
        breaks = c("c1f", "c3f"),
        labels = c("Cluster 1", "Cluster 3")
    ) +
    coord_flip(ylim = c(0, 1)) +
    theme(legend.position = "top")

# rat frequency plot ------------------
myWordCounthclustRat <- map(myCharcthclust$rat, getWordCount)

myWordCounthclustRatMerged <- full_join(myWordCounthclustRat[[1]], myWordCounthclustRat[[2]], by = "term") %>%
    full_join(myWordCounthclustRat[[3]], by = "term") %>%
    full_join(myWordCounthclustRat[[4]], by = "term")
colnames(myWordCounthclustRatMerged) <- c("term", "c1", "c2", "c3", "c4")

myWordCounthclustRatMerged$c1[is.na(myWordCounthclustRatMerged$c1)] <- 0
myWordCounthclustRatMerged$c2[is.na(myWordCounthclustRatMerged$c2)] <- 0
myWordCounthclustRatMerged$c3[is.na(myWordCounthclustRatMerged$c3)] <- 0
myWordCounthclustRatMerged$c4[is.na(myWordCounthclustRatMerged$c4)] <- 0

map(expPerhclustCluster$rat, length)
myWordCounthclustRatMerged <- mutate(
    myWordCounthclustRatMerged,
    c1f = c1/length(expPerhclustCluster$rat[[1]]),
    c2f = c2/length(expPerhclustCluster$rat[[2]]),
    c3f = c3/length(expPerhclustCluster$rat[[3]]),
    c4f = c4/length(expPerhclustCluster$rat[[4]])
)

myWordCounthclustRatMerged <- mutate(
    myWordCounthclustRatMerged,
    diff_c1_c3 = c1f - c3f
) %>% arrange(diff_c1_c3)

myWordCounthclustRatMergedClean <- filter(
    myWordCounthclustRatMerged,
    !(term %in% c("diet", "strain", "tissue", "gender", "liver", "for",
                  "dose", "mgkg", "route", "genotype", "months", "age",
                  "mixed", "high", "wks", "weeks", "hrs", "compound",
                  "sex", "8wk", "days", "treatment", "background", "strainbackground",
                  "status", "disease", "the", "genotypevariation", "genetic", "experiment",
                  "start", "rat", "rats", "norvegicus", "with", "time", "hour", "µgkg", "group",
                  "and", "sacrificed", "from", "feeded", "wky", "fed", "source", "mlkg"))
)
plot(myWordCounthclustRatMergedClean$diff_c1_c3)

ggpalette <-  colorspace::rainbow_hcl(4, c=90, l=50)

p_rat <- bind_rows(list(
    filter(
        myWordCounthclustRatMergedClean,
        term %in% c("male", "female")
    ),
    head(myWordCounthclustRatMergedClean, 5),
    tail(myWordCounthclustRatMergedClean, 5)
)) %>%
    distinct %>%
    select(term, c1f, c3f) %>%
    gather("cluster", "frequency", -term) %>%
    mutate(term = factor(term, levels = rev(unique(c("female", "male", unique(term)))))) %>%
    ggplot(aes(x = term, y = frequency, fill = cluster)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.6) +
    labs(x = "", y = "Frequency", title = "Rat") +
    coord_cartesian(xlim = c(0, 1)) +
    scale_fill_manual(
        values = c(ggpalette[2], ggpalette[3]),
        name = "",
        breaks = c("c1f", "c3f"),
        labels = c("Cluster 1", "Cluster 3")
    ) +
    coord_flip(ylim = c(0, 1)) +
    theme(legend.position = "top")

# plot export -------------------------
mp <- plot_grid(p_mouse, p_rat, ncol = 2, align = "hv")
ggsave("plots/fig2_frequency_hclust.svg", mp, device = svglite, width = 10, height = 5, units = "cm", scale = 1.8)
