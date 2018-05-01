

## Dependencies
0. python 3, R 3.4.4

1. Snakemake (Install using conda)

2. Kallisto (Quantification via pseudoalignment)
Complete installation instructions here - https://pachterlab.github.io/kallisto/download

3. [Sleuth (R Package)](https://github.com/pachterlab/sleuth)

4. Other R packages - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), 
ggplot2, dplyr, pheatmap, reshape2, [org.EcK12.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)

## Usage

`snakemake --snakefile snakefile.py`

## Input files

Fastqs go in`ecoli/cdna`
`metadata.tsv` in the home directory should contain two columns - samples, groups. Example included

## Outputs

`outs/counts/` contains abundance estimations calculated using Kallisto
`outs/plots/` contains visualizations for correlation, volcano plots and pathway enrichment
`outs/objects/` contains sleuth objects that can be used in a shiny app using `sleuth::sleuth_live(obj)`, 
differential expression tables, and normalized TPM matrices
