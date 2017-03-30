setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

library(dplyr)
library(readr)
library(purrr)
library(future); plan(multiprocess(workers = 1 + 32))
library(GEOquery)

seriesInfo <- list(
    mouse = read_tsv("data/GEOMouse.txt", progress = FALSE),
    rat = read_tsv("data/GEORat.txt", progress = FALSE)
)

getCharacterics <- function(gse) {
    gdt <- getGEO(gse, GSEMatrix = FALSE, destdir = "geo/")
    mch <- map(GSMList(gdt), function(x) {
        mr <- Meta(x)
        mr <- mr$characteristics_ch1 %>% paste(collapse = " ")
        return(mr)
    })
    message(gse)
    return(data_frame(
        gsm_id = names(mch),
        characteristics = unlist(mch)
    ))
}

# a few minutes
# getGEO can be quite obnoxious, try a few hours later if frozen
# delete temp files first
mouseCharacteristics <- future_lapply(
    unique(seriesInfo$mouse$GEO_Accession),
    getCharacterics
)
mouseCharacteristics <- bind_rows(mouseCharacteristics)
ratCharacteristics <- lapply(
    unique(seriesInfo$rat$GEO_Accession),
    getCharacterics
)
ratCharacteristics <- bind_rows(ratCharacteristics)

write.table(
    mouseCharacteristics,
    file = "data/mouseCharacteristics.txt",
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
)
write.table(
    ratCharacteristics,
    file = "data/ratCharacteristics.txt",
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
)


