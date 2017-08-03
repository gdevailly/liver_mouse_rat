library(GEOmetadb)
library(dplyr)
library(readr)
library(purrr)

setwd("/groups2/joshi_grp/guillaume/otherProject/pia_liver/")

t0 <- Sys.time()
getSQLiteFile(destdir = "data/")
Sys.time() - t0

t0 <- Sys.time() # 25 sec
normData <- list(
    mouse = read_tsv("data/heart_quantNormData_mouse.tsv", progress = FALSE),
    rat = read_tsv("data/heart_quantNormData_rat.tsv", progress = FALSE)
)
Sys.time() - t0

gseTable <- list(
    mouse = data_frame(gsm = colnames(normData$mouse[-1])),
    rat = data_frame(gsm = colnames(normData$rat[-1]))
)

con <- dbConnect(SQLite(), "data/GEOmetadb.sqlite")
dbListTables(con)
dbListFields(con, "gse_gsm")

mouseGseTable <- dbGetQuery(
    con,
    paste0(
        "SELECT * FROM gse_gsm WHERE gsm IN (",
        paste(
            paste0("'", gseTable$mouse$gsm, "'"),
            collapse = ", "
        ),
        ")"
    )
)

ratGseTable <- dbGetQuery(
    con,
    paste0(
        "SELECT * FROM gse_gsm WHERE gsm IN (",
        paste(
            paste0("'", gseTable$rat$gsm, "'"),
            collapse = ", "
        ),
        ")"
    )
)

# some GSM are in several different gse...
map_int(gseTable, nrow)
nrow(mouseGseTable)
nrow(ratGseTable)
nrow(distinct(mouseGseTable))
nrow(distinct(ratGseTable))

gseTables <- list(
    mouse = group_by(mouseGseTable, gsm) %>%
        sample_n(1) %>%
        ungroup(),
    rat = group_by(ratGseTable, gsm) %>%
        sample_n(1) %>%
        ungroup()
)

map_int(gseTable, nrow)
map_int(gseTables, nrow) # yes !

write_tsv(gseTables$mouse, path = "data/gse_gsm_mouse.tsv")
write_tsv(gseTables$rat, path = "data/gse_gsm_rat.tsv")

dbDisconnect(con)
