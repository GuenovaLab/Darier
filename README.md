# Darier
[![DOI](https://zenodo.org/badge/654606763.svg)](https://zenodo.org/badge/latestdoi/654606763)

A single-cell RNA analysis of in-house Darier disease integrated with Psoriasis and Healthy samples (GEO portal number GSE162183). 




## Usage
In order to run the analysis, first put the Darier matrices as well as the Psoriasis and Healthy matrices from GSE162183 in an "input" folder. Then, run the R script.

## Methods

### Single-cell bioinformatic analysis
Raw sequencing reads were subjected to demultiplexing and aligned to the human reference genome (refdata-gex-GRCh38-2020-A) using the cellranger count tool with default parameters. The analysis included a combination of publicly available dataset (GEO portal number GSE162183) (1) comprising psoriasis and healthy control samples, as well as an in-house dataset consisting of four Darrier patients. The datasets were processed together by merging the expression matrices based on common genes.
To ensure data quality, cells with a read mapping rate of over 25% to mitochondrial genes, which is indicative of dying cells, were filtered out. The remaining cells were then subjected to further analysis using the Seurat package (version 4.3.0) (2). Integration of the two datasets was performed using Harmony (version 0.1.1) (3), with the retention of the top 50 principal components.
Subsequently, cell clusters were identified using the Louvain algorithm implemented in Seurat with a resolution parameter of 0.1. To assign biological annotations, known marker genes were utilized. In addition, the T cell cluster was subjected to additional clustering to distinguish T helper and cytotoxic T cell subsets.

(1) Gao Y, Yao X, Zhai Y, Li L et al. Single cell transcriptional zonation of human psoriasis skin identifies an alternative immunoregulatory axis conducted by skin resident cells. Cell Death Dis 2021 May 6;12(5):450. PMID: 33958582   

(2) Integrated analysis of multimodal single-cell data, Hao et al., Cell, 2020  

(3) Fast, sensitive and accurate integration of single-cell data with Harmony, Korsunsky et al., Nature Methods, 2019

