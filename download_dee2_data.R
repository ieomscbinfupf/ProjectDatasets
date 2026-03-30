suppressPackageStartupMessages({
library(SummarizedExperiment)
library(getDEE2)
})

options(timeout=600)

tab <- readRDS("selGEOstudies.rds")

dee2 <- read.csv(gzfile("dee2_hsapiens_metadata.csv.gz"),
                 stringsAsFactors=FALSE, row.names=1)
dee2 <- dee2[dee2$GEO_series %in% tab$GSE, ]
maskQC <- substr(dee2$QC_summary, 1, 4) %in% c("PASS", "WARN")

dir.create("GSE")
for (i in seq_along(tab$GSE)) {
  gseid <- tab$GSE[i]
  cat(sprintf("GSE %s (%d/%d)\n", gseid, i, nrow(tab)))
  if (!file.exists(sprintf("GSE/%s.rds", gseid))) {
    mask <- (dee2$GEO_series %in% gseid) & maskQC
    srrlist <- dee2[mask, "SRR_accession"]
    tryCatch({
      gseobj <- getDEE2("hsapiens", srrlist, metadata=dee2, outfile=sprintf("GSE/%s", gseid))
      saveRDS(gseobj, file=sprintf("GSE/%s.rds", gseid))
    }, error=function(e) {
      cat(sprintf("%s FAILED!!\n", gseid))
      unlink(sprintf("GSE/%s.rds", gseid))
      unlink(sprintf("GSE/%s.zip", gseid))
    })
  }
}
