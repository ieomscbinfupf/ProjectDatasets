INSTALL httpfs;

-- SET VARIABLE url_geo_series = 'https://data-omicidx.cancerdatasci.org/geo/parquet/geo_series.parquet';
-- SET VARIABLE url_geo_samples = 'https://data-omicidx.cancerdatasci.org/geo/parquet/geo_samples.parquet';
-- SET VARIABLE url_geo_platforms = 'https://data-omicidx.cancerdatasci.org/geo/parquet/geo_platforms.parquet';
-- SET VARIABLE url_sra_accessions = 'https://data-omicidx.cancerdatasci.org/sra/parquet/sra_accessions.parquet';

-- GEO SERIES
CREATE VIEW geo_series AS SELECT accession, title, type,
            pubmed_id, sample_id, sample_taxid, platform_id,
            platform_taxid FROM read_parquet('https://data-omicidx.cancerdatasci.org/geo/parquet/geo_series.parquet');
--            platform_taxid FROM read_parquet(getvariable('url_geo_series'));
SELECT COUNT(*) FROM geo_series; -- 280008

-- GEO SAMPLES
CREATE VIEW geo_samples AS SELECT accession, type, platform_id,
            channels[1].molecule AS molecule, library_source
            FROM read_parquet('https://data-omicidx.cancerdatasci.org/geo/parquet/geo_samples.parquet');
--            FROM read_parquet(getvariable('url_geo_samples'));
SELECT COUNT(*) FROM geo_samples; -- 8384414

-- GEO PLATFORMS
CREATE VIEW geo_platforms AS SELECT *
            FROM read_parquet('https://data-omicidx.cancerdatasci.org/geo/parquet/geo_platforms.parquet');
--            FROM read_parquet(getvariable('url_geo_platforms'));
SELECT COUNT(*) FROM geo_platforms; -- 28186

-- select Illumina platforms for human
CREATE VIEW ilmna_platforms AS SELECT accession
            FROM geo_platforms WHERE title SIMILAR TO '.*Illumina .+ \(Homo sapiens\).*';

-- select samples sequenced in human Illumina platforms whose library
-- prep was with polyA RNA for bulk RNA-seq
CREATE VIEW bulk_samples AS SELECT geo_samples.accession AS accession
            FROM geo_samples JOIN ilmna_platforms ON
            geo_samples.platform_id=ilmna_platforms.accession
            WHERE geo_samples.molecule ILIKE '%polyA RNA%' AND
            geo_samples.library_source ILIKE '"transcriptomic"';
SELECT COUNT(*) FROM bulk_samples; -- 328480

-- SRA GEO ACCESSION
CREATE VIEW sra_geo_accessions AS SELECT accession, alias
            FROM read_parquet('https://data-omicidx.cancerdatasci.org/sra/parquet/sra_accessions.parquet');
--            FROM read_parquet(getvariable('url_sra_accessions'));
SELECT COUNT(*) FROM sra_geo_accessions; -- 145009027

-- dump to a CSV file
COPY (SELECT * FROM sra_geo_accessions WHERE accession LIKE 'SRR%' AND alias LIKE 'GSM%') TO 'SRR2GSM.csv' (header, delimiter ',');

-- select distinct GEO series identifiers occurring in human DEE2 data only
CREATE VIEW dee2gse AS SELECT DISTINCT GEO_series FROM read_csv('dee2_hsapiens_metadata.csv', delim=',') WHERE GEO_series LIKE 'GSE%';
SELECT COUNT(*) FROM dee2gse; -- 22794

-- subset GEO series to those accession occurring in dee2 data only
CREATE VIEW geo_series1 AS SELECT * FROM geo_series JOIN dee2gse ON geo_series.accession=dee2gse.GEO_series;
SELECT COUNT(*) FROM geo_series1; -- 22788

-- unnest platform taxonomic id
CREATE VIEW geo_series2 AS SELECT accession, title, type, pubmed_id,
       unnest(platform_id) AS platform_id, sample_id, platform_taxid
       FROM geo_series1;
SELECT COUNT(*) FROM geo_series2; -- 27010

-- select GEO series from human Illumina platforms
CREATE VIEW geo_series3 AS SELECT geo_series2.accession AS accession,
       title, type, pubmed_id, sample_id, platform_taxid
       FROM geo_series2 JOIN ilmna_platforms
       ON geo_series2.platform_id=ilmna_platforms.accession;
SELECT COUNT(*) FROM geo_series3; -- 22239

-- unnest platform_taxid
CREATE VIEW geo_series4 AS SELECT accession, title, type, pubmed_id,
       sample_id, unnest(platform_taxid) AS platform_taxid
       FROM geo_series3;
SELECT COUNT(*) FROM geo_series4; -- 23982

-- unnest type selecting GEO series with platform_taxid=9606 (human)
CREATE VIEW geo_series5 AS SELECT accession, title,
       unnest(type) as type, pubmed_id, sample_id
       FROM geo_series4 WHERE platform_taxid=9606;
SELECT COUNT(*) FROM geo_series5; -- 26279

CREATE VIEW geo_series6 AS SELECT accession, title,
       pubmed_id[1] AS pubmed_id, sample_id
       FROM geo_series5 WHERE
       type ILIKE 'Expression profiling by high throughput sequencing'
       AND pubmed_id[1] IS NOT NULL;
SELECT COUNT(*) FROM geo_series6; -- 16527

-- unnest sample_id
CREATE VIEW geo_series7 AS SELECT accession, title, pubmed_id,
       unnest(sample_id) AS sample_id FROM geo_series6;
SELECT COUNT(*) FROM geo_series7; -- 1042549

-- select GSE series and samples from polyA bulk RNA-seq 
CREATE VIEW geo_series8 AS SELECT geo_series7.accession AS accession,
       title, pubmed_id, sample_id FROM geo_series7 JOIN bulk_samples
       ON geo_series7.sample_id=bulk_samples.accession;
SELECT COUNT(*) FROM geo_series8; -- 223102

-- dump to a CSV file
-- TODO add the corresponding SRA ID column and propagate it through the pipeline
COPY (SELECT * FROM geo_series8) TO 'bulk_rnaseq_samples.csv' (header, delimiter ',');
