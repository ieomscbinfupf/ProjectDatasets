# Project Datasets

This is the project datasets repository for the IEO subject taught in the
MSc program on Bioinformatics for the Health Sciences at the
Universitat Pompeu Fabra. The goal of this repo is to provide a pipeline to
select publicly available human bulk RNA-seq datasets and put them available to
the students through the URL https://ieomscbinfupf.github.io/ProjectDatasets, so
that the students may choose what dataset want to analyze according to their
preferences of topic (disease, physiological system, molecular mechanism, etc.).

The metadata and phenotype data on these publicly available data sets are
retrieved from the
[Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo), through the
[omicIDX](https://github.com/omicidx) database. The actual raw RNA-seq counts
are retrieved from the
[Digital Expresssion Explorer 2 (DEE2)](https://dee2.io). These human bulk
RNA-seq dataset are selected on the basis of the following criteria:

* They are described in a published article.
* The library preparation selected polyA+ RNA molecules.
* They have a minimum of 10 samples sequenced in Illumina instruments and whose
  sequencing output also meets a minimum standard of quality in terms of depth
  and mapability to the human genome.
* The raw sequencing reads have been processed by the [DEE2](https://dee2.io)
  pipeline, from which the pipeline uses the reported
  [quality metrics](https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md)
  to exclude potentially problematic samples, and downloaded the processed table
  of raw counts.

## Pipeline

1. Download DEE2 metadata that describes the available processed bulk RNA-seq
   samples for human.

```sh
$ curl --output dee2_hsapiens_metadata.tsv https://dee2.io/metadata/hsapiens_metadata.tsv
$ R --vanilla -e 'dat <- read.csv("dee2_hsapiens_metadata.tsv", sep="\t") ; write.csv(dat, "dee2_hsapiens_metadata.csv")'
$ gzip dee2_hsapiens_metadata.csv
$ rm dee2_hsapiens_metadata.tsv
```

2. Process the omicIDX databases with [duckdb](https://duckdb.org). This step will
create two output files: `SRR2GSM.csv` and `bulk_rnaseq_samples.csv`.

```sh
$ duckdb < process_omicidx.sql
$ gzip SRR2GSM.csv
$ gzip bulk_rnaseq_samples.csv
```

2. Select samples meeting some minimum quality standards and build a table
including bibliographic information on the associated publications (this step
requires the R packages, XML, rcrossref and rentrez). This step will create an
output file called `selGEOstudies.rds`.

```sh
$ Rscript select_project_datasets.R
```

3. Download DEE2 data.

```sh
$ Rscript download_dee2_data.R
```

4. Fetch feature annotation for Ensembl v90, which is the version used in DEE2,
and map them to Entrez genes. This step will create an output file called
`ensgenesv90_entrez.rds`.

```sh
$ Rscript fetch_feature_annotation.R
```

5. Download corresponding phenotypic data from GEO and create `SummarizedExperiment`
objects. This step will create a directory called `SEs` with all the
`SummarizedExperiment` objects.

```sh
$ Rscript download_phenogeo_data_create_SEs.R
``` 

6. Build table with the project datasets. This step will create a file called
`SEs/index.html` with the project datasets table.

```sh
$ Rscript build_project_datasets_table.R
```
