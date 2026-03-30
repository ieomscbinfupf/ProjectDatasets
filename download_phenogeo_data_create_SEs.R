suppressPackageStartupMessages({
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(GEOquery)
})

options(timeout=600)

srr2gsm <- read.csv(gzfile("SRR2GSM_2026.csv.gz"))
dim(srr2gsm)
## [1] 7112962      2
srr2gsm[[2]] <- gsub("[\\._].+$", "", srr2gsm[[2]])

rr <- readRDS("ensgenesv90_entrez.rds")

## set rownames to the Entrez gene identifiers and remove
## unnecessary columns
ensids <- names(rr)
names(rr) <- rr$entrezid
rr$seq_coord_system <- rr$gene_name <- rr$entrezid <- NULL


dir.create("SEs")
fnames <- list.files(path="GSE", pattern="GSE.*.rds$")
for (i in seq_along(fnames)) {
  gseid <- sub(".rds", "", fnames[i])
  cat(sprintf("GSE %s (%d/%d)\n", gseid, i, length(fnames)))
  if (!file.exists(file.path("SEs", paste0(gseid, ".rds")))) {
    cat(sprintf("%s\n", gseid))
    gse <- readRDS(file.path("GSE", fnames[i]))
    eset <- tryCatch(getGEO(gseid),
                     error=function(cond) {
                       cat(sprintf("couldn't download GEO data for %s\n", gseid))
                       return(NULL)
                     })

    if (!is.null(eset)) {
      cnt <- as.matrix(assay(gse))
      cnt <- cnt[ensids, ]
      stopifnot(all(rr$entrezid == ensids)) ## QC
      stopifnot(all(rr$entrezid == names(rr))) ## QC
      rownames(cnt) <- names(rr)
      mt <- match(colnames(cnt), srr2gsm[[1]])
      if (any(is.na(mt)))
        cat(sprintf("%d SRRs could not be mapped to GSMs (%s)\n", sum(is.na(mt)), paste(colnames(cnt)[is.na(mt)], collapse=", ")))

      if (all(is.na(mt)))
        cat(sprintf("No SRRs could be mapped to GSMs, skipping %s\n", gseid))
      else {
        gsmids <- srr2gsm[[2]][mt[!is.na(mt)]]
        if ("list" %in% class(eset)) {
          ## if GEO returned more than one ExpressionSet, select the one
          ## matching the samples to the phenotypic data
          mask <- sapply(eset, function(es, ids) all(ids %in% sampleNames(es)), gsmids)
          if (any(mask))
            eset <- eset[[which(mask)[1]]]
          else {
            cat(sprintf("%s doesn't have GEO data matching the samples\n", gseid))
            eset <- NULL
          }
        } else if ("ExpressionSet" %in% class(eset)) {
          ## check whether samples match to the GEO phenotypic data
          if (!all(gsmids %in% sampleNames(eset))) {
            cat(sprintf("%s doesn't have GEO data matching the samples\n", gseid))
            eset <- NULL
          }
        } else
          stop(sprintf("GEO returned an object class different from ExpressionSet: %s\n", class(eset)))

        if (!is.null(eset)) {
          cdata <- pData(eset)[gsmids, ]
          rownames(cdata) <- colnames(cnt)
          cdata <- cdata[, -grep("contact_", colnames(cdata))]

          mdata <- list(experimentData=experimentData(eset),
                        annotation="org.Hs.eg.db",
                        ensemblVersion="v90",
                        urlProcessedData=sprintf("https://dee2.io/cgi-bin/search.sh?org=hsapiens&accessionsearch=%s",
                                                 gse$SRP_accession[1]))
          stopifnot(all(colnames(cnt) == rownames(cdata))) ## QC
          stopifnot(all(rownames(cnt) == names(rr))) ## QC

          seobj <- SummarizedExperiment(assays=list(counts=cnt),
                                        rowRanges=rr,
                                        colData=cdata,
                                        metadata=mdata)
          saveRDS(seobj, file=file.path("SEs", paste0(gseid, ".rds")))
        }
      }
    }
  } else cat(sprintf("skipping %s\n", gseid))
}
