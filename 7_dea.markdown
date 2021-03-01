---
layout: page
title: 7. Phân tích khác biệt biểu hiện gene (Differential expression analysis)
permalink: /dea/
nav_order: 7
---

Source: 

[https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#dechapter](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#dechapter)

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/05-diffexp.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/05-diffexp.html)

## Differential expression

Now that we've assigned cells into clusters, we'd like to understand what makes each cluster different from other cells in the dataset, or to annotate clusters according to their cell types (as has been previously done for this dataset).

There are several approaches to this task:

- Look for upregulation of marker genes for cell types of interest (compared to the rest of the dataset)
- Compare the complete gene expression profiles between groups
- Use automated methods to compare cells of interest to databases of cell type expression profiles to combine clustering and annotation

Automated methods are a promising advance, but are not yet able to replace careful human curation.

For well-defined cell types, we expect marker genes to show large differences in expression between the cell type of interest and the rest of the dataset, allowing us to use simple methods. We'll focus on this approach for this workshop, while building intuition that is broadly applicable to other approaches.

## Comparing distributions 

Differential expression algorithms represent various approaches to comparing the distribution of gene expression in one group versus another group. Unlike bulk RNA-seq, we generally have a large number of samples (i.e. cells) for each group we are comparing in single-cell experiments. Thus, we can take advantage of the whole distribution of expression values in each group to identify differences between groups rather than only comparing estimates of mean-expression as is standard for bulk RNASeq.

Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

The min.pct argument requires a gene to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a gene to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of genes that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed genes will likely still rise to the top.


```R
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
```

```console
            p_val avg_logFC pct.1 pct.2    p_val_adj
IL7R 1.963908e-63 1.2187234 0.893 0.292 6.586554e-59
TRAC 1.591017e-53 0.9225490 0.932 0.380 5.335953e-49
IL32 7.534694e-52 0.9455882 0.927 0.328 2.526986e-47
LDHB 1.350431e-51 0.8313688 0.971 0.714 4.529077e-47
CD2  1.015116e-50 0.7672408 0.791 0.238 3.404496e-46
```

```R
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))
```
```console
             p_val avg_logFC pct.1 pct.2     p_val_adj
CD3E 6.885642e-121  1.889710 0.988 0.016 2.309307e-116
CD3D 8.068824e-121  1.863903 0.975 0.009 2.706122e-116
TRAC 6.004612e-105  1.930913 0.963 0.067 2.013827e-100
IL7R 3.483779e-104  1.961180 0.894 0.018  1.168390e-99
TCF7 3.630007e-103  1.670272 0.925 0.046  1.217432e-98
```

```R
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
```

```console
Calculating cluster 0
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=05s  
Calculating cluster 1
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Calculating cluster 2
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Calculating cluster 3
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Calculating cluster 4
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Calculating cluster 5
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Calculating cluster 6
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
Calculating cluster 7
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s
```

Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

```R
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

We include several tools for visualizing marker expression. • VlnPlot (shows expression probability distributions across clusters), • and FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring: • RidgePlot, • CellPlot, and • DotPlot as additional methods to view your dataset.

```R
VlnPlot(object = pbmc, features =c("LYZ", "VCAN"))

FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")
```

DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.


```R
library(dplyr)

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, features = top5$gene, label = TRUE)
```

## Assigning cell type identity to clusters

Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types.

```R
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5)
```