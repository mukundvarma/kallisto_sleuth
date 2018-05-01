# RNA Seq pipeline for pseudoalignment, quantification of reads and differential expression of genes

## Dependencies
0. python 3, R 3.4.4
1. [Snakemake](https://snakemake.readthedocs.io/en/stable/) (Can be installed using conda)
2. Kallisto (Quantification via pseudoalignment)
Complete installation instructions here - https://pachterlab.github.io/kallisto/download, can also be installed using conda if the bioconda channel is active
3. [Sleuth (R Package)](https://github.com/pachterlab/sleuth)
4. Other R packages - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), 
ggplot2, dplyr, pheatmap, reshape2, [org.EcK12.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)

## Usage

`snakemake --snakefile snakefile.py`

To recreate output, just install dependencies, remove the outs directory and run this command.

## Input files

Fastqs go in`ecoli/cdna`

`metadata.tsv` goes in the home directory and should contain two columns - samples, groups. Example included

## Outputs

`outs/counts/` contains abundance estimations calculated using Kallisto

`outs/objects/` contains a sleuth object that can be visualized in a shiny app to see QC metrics,
differentially expressed genes, correlation heatmaps etc.

~~~~
library(sleuth)
obj = readRDS('outs/objects/sleuth.object.RDS')
sleuth_live(obj)
~~~~

This directory also contains a normalized TPM matrix (gene expression level), and a differential expression
output table from *sleuth*'s implementation of the Wald test.

`outs/plots/` contains additional visualizations for correlation, volcano plots and pathway enrichment

## Acknowledgements

Snakemake code heavily borrowed from github.com/slowkow and github.com/saketkc
