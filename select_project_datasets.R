dat <- read.csv(gzfile("bulk_rnaseq_samples.csv.gz"))

gsmbygse <- split(dat$sample_id, dat$accession)
length(gsmbygse)
## [1] 2752
mt <- match(names(gsmbygse), dat$accession)
tab <- data.frame(GSE=names(gsmbygse),
                  PubmedID=dat$pubmed_id[mt],
                  Title=dat$title[mt],
                  NSamples=lengths(gsmbygse))

dee2 <- read.csv(gzfile("dee2_hsapiens_metadata.csv.gz"),
                 stringsAsFactors=FALSE, row.names=1)
dim(dee2)
## [1] 768802    23
dee2 <- dee2[dee2$GEO_series %in% tab$GSE, ]
dim(dee2)
## [1] 115979    23
dee2qcsummary <- split(dee2$QC_summary, dee2$GEO_series)
length(dee2qcsummary)
## [1] 2752
mask <- sapply(dee2qcsummary, function(x) all(substr(x, 1, 4) == "PASS" | substr(x, 1, 4) == "WARN"))
sum(mask)
## [1] 1936

## attempt to count samples at the level of biological replicates
nactualsamples <- lengths(lapply(split(dee2$sample_alias, dee2$GEO_series), unique))
stopifnot(identical(names(mask), names(nactualsamples))) ## QC
sum(mask & nactualsamples >= 10)
## [1] 631

selectedGSE <- names(dee2qcsummary)[mask & nactualsamples >= 10]
length(selectedGSE)
## [1] 631

stopifnot(all(selectedGSE %in% tab$GSE)) ## QC
dim(tab)
## 2752    4
tab <- tab[tab$GSE %in% selectedGSE, ]
dim(tab)
## [1] 631   4
tab$NActualSamples <- nactualsamples[tab$GSE]

## add publication information
source("pubmed.R")
cits <- fetchCitations(tab$PubmedID[1:300])
cits2 <- fetchCitations(tab$PubmedID[301:nrow(tab)])
cits <- rbind(cits, cits2)
stopifnot(identical(cits$PMID, tab$PubmedID)) ## QC
tab$Title <- cits$TITLE
tab$DOI <- cits$DOI
tab$Date <- cits$DATE
tab$Citations <- cits$CIT
tab$NSamples <- tab$NActualSamples
tab$NActualSamples <- NULL

saveRDS(tab, file="selGEOstudies.rds")
tab <- readRDS("selGEOstudies.rds")
