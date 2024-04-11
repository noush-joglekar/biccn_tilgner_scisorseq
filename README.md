# biccn_tilgner_scisorseq

Code for the BRAIN Initiative Cell Census Network (BICCN) - [RF1 Tilgner](https://biccn.org/teams/rf1-tilgner)

Until recently, technological limitations prevented a genome-wide appraisal of the influence of isoform expression
on the establishment and persistence of cell identity in various parts of the brain. Using an enhanced method of
long-read single-cell isoform sequencing (ScISOr-Seq2), we performed a comprehensive analysis of RNA isoforms in
mouse brain across multiple brain regions, cell subtypes, and a range of developmental timepoints from postnatal
day 14 (P14) to adult (P56).

A detailed summary of the samples processed and links to raw data can be found on the [BICCN portal](https://knowledge.brain-map.org/data/Z0GBA7V12N4J4NNSUHA/summary).
Nonetheless we are providing alist of samples and the data types available for each below:

| Age | Sex | Region | Replicates | 10x scRNA-Seq | PacBio | Oxford Nanopore |
| --- | --- | ------ | ---------- | ---------- | ------- | -------- |
| P14 | M | Hippocampus | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P21 | M | Hippocampus | 2 | :heavy_check_mark: | x | :heavy_check_mark: |
| P28 | M | Hippocampus | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P56 | M | Hippocampus | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P14 | M | VisCortex | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P21 | M | VisCortex | 2 | :heavy_check_mark: | x | :heavy_check_mark: |
| P28 | M | VisCortex | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P56 | M | VisCortex | 2 | :heavy_check_mark: | x | :heavy_check_mark: |
| P56 | M | Striatum | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P56 | M | Thalamus | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| P56 | M | Cerebellum | 2 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |


An article describing this data and our findings can be found [here](https://www.nature.com/articles/s41593-024-01616-4)

Custom code code used for processing and summarizing this data, files needed for this processing,
and code used to generate plots for the manuscript is provided in this repo.

Additionally, the following raw and processed data is /will be made available on [NeMO](https://nemoarchive.org/)
Links to mouse data:
- Short read Illumina reads [fastqs](https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/10x_v3/mouse/raw/)
- Cellranger output files [mtx, tsv](https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/10x_v3/mouse/processed/counts/),
- Raw long-read sequencing data sequenced on (PacBio)[https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/] and (ONT)[https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/oxford_nanopore/mouse/raw/]
- AllInfo files per processed (ONT)[https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/oxford_nanopore/mouse/processed/other/] long read containing
  - cell barcode
  - assigned cell type
  - UMI
  - intron chain
  - TSS assignment
  - PolyA site assignment
  - exon chain
- BED12 files per processed (ONT)[https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/oxford_nanopore/mouse/processed/other/] long read

We additionally have data from 6 human hippocampal samples available through the (NeMO)[https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/sncell/]

If you have questions, please do not hesitate to open an issue or reach out over email to [me](mailto:ajoglekar@nygenome.org)
or the corresponding author.
