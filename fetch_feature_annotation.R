suppressPackageStartupMessages({
  library(AnnotationHub)
  library(org.Hs.eg.db)
  library(ensembldb)
  library(GenomeInfoDb)
  library(SummarizedExperiment)
})

## get feature information

## retrieve Ensembl annotations v90
## matching the version reported in
## https://github.com/markziemann/dee2/blob/master/misc/gene_info/curate_gene_tx_info.sh
## fetchTablesFromEnsembl(90, species="human")

ah <- AnnotationHub()
ahEnsDb <- query(ah, c("Homo Sapiens", "EnsDb", 90))
ahEnsDb <- ahEnsDb[[1]]
ensgenes <- genes(ahEnsDb)
length(ensgenes)
## [1] 64661
stopifnot(all(!duplicated(names(ensgenes)))) ## QC
stopifnot(all(!is.na(ensgenes$symbol))) ## QC

## discard alternate chromosome sequences and patches
ensgenes <- keepStandardChromosomes(ensgenes, pruning.mode="coarse")
length(ensgenes)
## [1] 58930

## we translate ENSEMBL gene identifiers into Entrez gene identifiers
## using gene SYMBOL from org.Hs.eg.db to select the most appropriate Entrez
## when a Ensembl gene is annotated to multiple Entrez genes with the
## same symbol assignment, then the entrez with identical symbol will be selected
entrezbyens <- ensgenes$entrezid
allentrez <- as.character(unlist(entrezbyens, use.names=FALSE))
symbolbyens <- rep(NA_character_, length(allentrez))
symbolbyens[!is.na(allentrez)] <- select(org.Hs.eg.db, keys=allentrez[!is.na(allentrez)],
                                         columns="SYMBOL", keytype="ENTREZID")$SYMBOL
symbolbyens <- relist(symbolbyens, entrezbyens)
mtsymbolbyens <- mapply(function(sbye, sym) sbye == sym,
                        symbolbyens, relist(rep(ensgenes$symbol, lengths(symbolbyens)), symbolbyens))
entrezbyens <- mapply(function(ebye, mtsbye) {
                        res <- ebye[1]
                        mtsbye[is.na(mtsbye)] <- FALSE
                        if (any(mtsbye))
                          res <- ebye[mtsbye]
                        res[1]
                      }, entrezbyens, mtsymbolbyens)
stopifnot(identical(names(ensgenes), names(entrezbyens))) ## QC
ensgenes$entrezid <- as.character(entrezbyens)
sum(duplicated(ensgenes$entrezid))
## [1] 33802

symbyentrez <- split(ensgenes$symbol, ensgenes$entrezid)
mask <- sapply(symbyentrez, function(x) all(x==x[1]))
length(mask)
## [1] 25127
sum(mask)
## [1] 25101

## so there are 26 (in 2026) different Ensembl genes with identical
## Entrez IDs but with different symbols. we will use
## org.Hs.eg.db again to decide which Entrez assignment
## do we select
syms <- select(org.Hs.eg.db, keys=names(symbyentrez)[!mask],
               columns="SYMBOL")$SYMBOL
names(syms) <- names(symbyentrez)[!mask]
wh <- which(ensgenes$entrezid %in% names(syms))
wh <- names(ensgenes)[wh][syms[ensgenes$entrezid[wh]] != ensgenes$symbol[wh]]

## discard Ensembl genes with the same Entrez ID but where the
## symbol assignment with org.Hs.eg.db does not match
ensgenes <- ensgenes[-match(wh, names(ensgenes))]
length(ensgenes)
## [1] 58892

## check that all Ensembl genes with identical Entrez IDs also
## have identical symbols
symbyentrez <- split(ensgenes$symbol, ensgenes$entrezid)
mask <- sapply(symbyentrez, function(x) all(x==x[1]))
length(mask)
## [1] 25118
stopifnot(all(mask)) ## QC

## finally resolve duplicated Entrez IDs by selecting the
## Ensembl gene with higher depth across all experiments
fnames <- list.files(path="GSE", pattern="GSE.*.rds$")
gse <- readRDS(file.path("GSE", fnames[1]))
cnt <- matrix(integer(0), nrow=nrow(gse), ncol=0,
              dimnames=list(rownames(gse), NULL))
for (i in seq_along(fnames)) {
  gse <- readRDS(file.path("GSE", fnames[i]))
  stopifnot(identical(rownames(gse), rownames(cnt))) ## QC
  cnt <- cbind(cnt, assay(gse))
}
dim(cnt)
## [1] 58302 15998
cnt <- cnt[rownames(cnt) %in% names(ensgenes), ]
dim(cnt)
## [1] 58205 15998
geneseqdepth <- rowSums(cnt)
ensbyeg <- split(names(ensgenes), ensgenes$entrezid)
lens <- lengths(ensbyeg)
ensbyeg <- sapply(ensbyeg, function(e, d) e[which.max(d[e])],
                  geneseqdepth)
length(ensbyeg)
## [1] 25118

stopifnot(all(!duplicated(ensbyeg))) ## QC
stopifnot(all(ensbyeg %in% names(ensgenes))) ## QC
stopifnot(all(ensbyeg %in% rownames(cnt))) ## QC

ensgenes <- ensgenes[ensbyeg]
length(ensgenes)
## [1] 25118

stopifnot(all(!duplicated(ensgenes$entrezid))) ## QC
sum(duplicated(ensgenes$symbol))
## [1] 2

mcols(ensgenes[ensgenes$symbol %in% ensgenes$symbol[duplicated(ensgenes$symbol)]])[, c("symbol", "entrezid")]
## DataFrame with 4 rows and 2 columns
##                      symbol    entrezid
##                 <character> <character>
## ENSG00000187600     TMEM247   101805491 ("LINC02583" manually checked in NCBI)
## ENSG00000253051     SNORA31   109616966 ("SNORA31B" manually checked in NCBI)
## ENSG00000284701     TMEM247      388946 ("TMEM247" manually checked in NCBI)
## ENSG00000199477     SNORA31      677814 ("SNORA31" manually checked in NCBI)
ensgenes["ENSG00000187600"]$symbol <- "LINC02583"
ensgenes["ENSG00000253051"]$symbol <- "SNORA31B"
stopifnot(all(!duplicated(ensgenes$symbol))) ## QC
stopifnot(all(names(ensgenes) %in% rownames(cnt))) ## QC
saveRDS(ensgenes, file="ensgenesv90_entrez.rds")
