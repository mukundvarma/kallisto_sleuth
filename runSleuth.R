library(sleuth)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
library(org.EcK12.eg.db)

### Sleuth Differential Expression Analysis

## Create metadata object
sample_id = dir(file.path("outs", "counts"))
outdirs = file.path("outs", "counts", sample_id)
s2c = read.table("metadata.tsv", header=T, stringsAsFactors=F)
s2c = dplyr::mutate(s2c, path = outdirs)

## Read in abundances file and create a sleuth object
so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

## Add Gene names
embl2gene = read.delim("reference/embl_to_genesymbol.tsv", header=F, stringsAsFactors=F)
colnames(embl2gene) = c("target_id", "gene")
so = sleuth_prep(s2c, target_mapping = embl2gene)

## Fit a model
so = sleuth_fit(so, ~group, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
so = sleuth_wt(so, 'groupstate2')
so = sleuth_wt(so, "groupstate2", which_model="full")

## Create a table of differentially expressed transcripts
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
waldresults = sleuth_results(so, 'groupstate2')
print(head(waldresults))

tpm_norm = sleuth_to_matrix(so, "obs_norm", "tpm")

## Calculate correlations and plot them
cor_kendall = cor(tpm_norm$data %>% as.matrix, method="kendall")
cor_pearson = cor(log(1+tpm_norm$data) %>% as.matrix, method="pearson")

pheatmap::pheatmap(cor_kendall, breaks = seq(-1,1,by=0.1),
                   color = colorRampPalette(c("firebrick3",
                                              "lightgrey",
                                              "navy"))(20) %>% rev,
                   main="Kendall's Tau", file = "outs/plots/correlation.kendall.pdf")

pheatmap::pheatmap(cor_pearson, breaks = seq(-1,1,by=0.1),
                   color = colorRampPalette(c("firebrick3",
                                              "lightgrey",
                                              "navy"))(20) %>% rev,
                   main="Pearson's coeff.", file = "outs/plots/correlation.pearson.pdf")

## Create counts QC plot
counts.melted = melt(tpm_norm$data)
colnames(counts.melted) = c("transcript", "sample", "tpm_norm")
ggplot(data = counts.melted, aes(x = tpm_norm)) + 
    geom_histogram(bins = 10) + facet_wrap(~sample, scales = 'free_x') +
    scale_x_log10() + theme_bw()

## Create differential expression volcano plot


## Save objects
saveRDS(so, "outs/objects/sleuth.object.RDS")
write.table(sleuth_table, "outs/objects/sleuth.differential_expression.tsv", quote=F, sep="\t")
write.table(tpm_norm, "outs/objects/abundances.tpm_norm.tsv", quote=F, sep="\t")


## More plotting
pdf("outs/plots/volcanoplot.diffexp.pdf")
ggplot(waldresults, aes(x=b, y=-log10(qval))) +
    geom_point(color="grey") +
    geom_point(data=subset(waldresults, qval<0.01),
               aes(x=b, y=-log10(qval), color="firebrick3")) +
    geom_text(data=subset(waldresults,qval < 0.01),
              aes(label=gene, x=b, y=-log10(qval))) +
    labs(title="Volcano plot of DE genes", subtitle="qval < 0.01 highlighted")
dev.off()

Up = waldresults %>% filter(b > 0 & qval < 0.01) %>% .$gene %>%
    bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb='org.EcK12.eg.db') %>% .$ENTREZID

Down = waldresults %>% filter(b < 0 & qval < 0.01) %>% .$gene %>%
    bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb='org.EcK12.eg.db') %>% .$ENTREZID

genelists = list(Up, Down)
names(genelists) = c("State2", "State1")

clustprof = compareCluster(genelists, fun="enrichGO", OrgDb='org.EcK12.eg.db', ont="BP")

pdf("outs/plots/GO_Enrichment.pdf", height=6, width=10)
dotplot(clustprof) + ggtitle("GO Pathway enrichment for DE genes")
dev.off()


