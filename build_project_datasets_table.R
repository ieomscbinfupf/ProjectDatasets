tab <- readRDS("selGEOstudies.rds")

## discard datasets for which either phenotypic data from GEO could not be retrieved,
geogseids <- gsub(".rds", "", list.files("SEs", pattern="*.rds"))
tab <- tab[tab$GSE %in% geogseids, ]

duppubs <- tab$PubmedID[duplicated(tab$PubmedID)]
stopifnot(identical(duplicated(tab$PubmedID), duplicated(tab$Title))) ## QC
tabdups <- tab[tab$PubmedID %in% duppubs, ]
tabuniq <- tab[!tab$PubmedID %in% duppubs, ]

dtf <- data.frame("Publication DOI & PMID"=paste(sprintf("<a href=\"https://doi.org/%s\" target=\"_blank\">%s</a>",
                                                         tabuniq$DOI, tabuniq$DOI),
                                                 sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term=%s\" target=\"_blank\">%s</a>",
                                                         tabuniq$PubmedID, tabuniq$PubmedID), sep=" - "),
                  Title=tabuniq$Title,
                  Date=format(tabuniq$Date, "%Y/%m/%d"),
                  Citations=as.integer(tabuniq$Citations),
                  "#Samples"=as.character(tabuniq$NSamples),
                  GEO=sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\" target=\"_blank\">%s</a>", tabuniq$GSE, tabuniq$GSE),
                  SEobj=sprintf("<a href=\"https://functionalgenomics.upf.edu/courses/IEO/projects/datasets/%s.rds\">link</a>", tabuniq$GSE),
                  DEE2=sprintf("<a href=\"https://functionalgenomics.upf.edu/courses/IEO/projects/datasets/%s.zip\">link</a>", tabuniq$GSE),
                  stringsAsFactors=FALSE, check.names=FALSE)

gsebypub <- split(tabdups$GSE, tabdups$PubmedID)
nsmbypub <- split(tabdups$NSamples, tabdups$PubmedID)
tabdups <- tabdups[!duplicated(tabdups$PubmedID), ]

dtfdups <- data.frame("Publication DOI & PMID"=paste(sprintf("<a href=\"https://doi.org/%s\" target=\"_blank\">%s</a>",
                                                         tabdups$DOI, tabdups$DOI),
                                                 sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term=%s\" target=\"_blank\">%s</a>",
                                                         tabdups$PubmedID, tabdups$PubmedID), sep=" - "),
                  Title=tabdups$Title,
                  Date=format(tabdups$Date, "%Y/%m/%d"),
                  Citations=as.integer(tabdups$Citations),
                  "#Samples"=sapply(nsmbypub[as.character(tabdups$PubmedID)], paste, collapse=", "),
                  GEO=sapply(gsebypub[as.character(tabdups$PubmedID)],
                             function(gseids) {
                               paste(sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\" target=\"_blank\">%s</a>", gseids, gseids), collapse=", ")
                             }),
                  SEobj=sapply(gsebypub[as.character(tabdups$PubmedID)],
                               function(gseids) {
                                 paste(sprintf("<a href=\"https://functionalgenomics.upf.edu/courses/IEO/projects/datasets/%s.rds\">link</a>", gseids), collapse=", ")
                               }),
                  DEE2=sapply(gsebypub[as.character(tabdups$PubmedID)],
                              function(gseids) {
                                paste(sprintf("<a href=\"https://functionalgenomics.upf.edu/courses/IEO/projects/datasets/%s.zip\">link</a>", gseids), collapse=", ")
                              }),
                  stringsAsFactors=FALSE, check.names=FALSE)

dtf <- rbind(dtf, dtfdups)

## order by date (most recent first)
dtfdate <- as.Date(dtf$Date, "%Y/%m/%d")
dtf <- dtf[order(dtfdate, decreasing=TRUE), ]
rownames(dtf) <- 1:nrow(dtf)
dtf[["#"]] <- 1:nrow(dtf)
dtf <- dtf[, c("#", colnames(dtf)[-ncol(dtf)])]

library(reactable)
library(htmltools)
library(htmlwidgets)

rtab <- reactable(dtf, bordered=TRUE, striped=TRUE, highlight=TRUE, searchable=TRUE,
                  defaultPageSize=25,
                  columns=list("#"=colDef(minWidth=60),
                               "Publication DOI & PMID"=colDef(html=TRUE, minWidth=110),
                               Title=colDef(minWidth=300),
                               Date=colDef(minWidth=120),
                               Citations=colDef(minWidth=100),
                               "#Samples"=colDef(minWidth=110),
                               GEO=colDef(html=TRUE, minWidth=150),
                               SEobj=colDef(html=TRUE, minWidth=100),
                               DEE2=colDef(html=TRUE, minWidth=100)),
                  theme=reactableTheme(borderColor = "#dfe2e5",
                                       stripedColor = "#f6f8fa",
                                       highlightColor = "#90ee90",
                                       cellPadding = "8px 12px",
                                       style = list(fontFamily = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"),
                                       searchInputStyle = list(width = "100%")
                  ))
rtab <- prependContent(rtab,
                       h1(class="title",
                          "Selected human polyA+ bulk RNA-seq datasets publicly available at NCBI GEO with at least 10 samples, described in one publication and meeting some minimum data quality standards",
                          style="font-family: -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"))
saveWidget(rtab, "SEs/index.html", selfcontained=TRUE)
