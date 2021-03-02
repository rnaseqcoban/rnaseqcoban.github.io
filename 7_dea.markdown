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
                p_val avg_logFC pct.1 pct.2     p_val_adj
S100A12 4.698018e-208  3.102009 0.967 0.015 1.575621e-203
VCAN    7.589032e-207  2.811492 0.980 0.029 2.545210e-202
MNDA    2.392414e-200  2.504461 0.993 0.058 8.023677e-196
MS4A6A  6.548670e-198  2.015251 0.967 0.030 2.196293e-193
FCN1    2.352186e-197  2.547167 0.990 0.049 7.888762e-193
```

```R
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster2.markers <- FindMarkers(object = pbmc, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster2.markers, n = 5))
```
```console
                 p_val avg_logFC pct.1 pct.2     p_val_adj
CD79A    1.048544e-131  2.850362 0.963 0.020 3.516606e-127
MS4A1    4.615700e-125  2.399839 0.936 0.026 1.548013e-120
BANK1    4.149180e-119  1.805717 0.888 0.016 1.391552e-114
HLA-DQA1 6.781299e-119  1.980109 0.973 0.053 2.274312e-114
IGHM     2.290497e-116  3.537477 0.882 0.022 7.681869e-112
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
```

Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

```R
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

We include several tools for visualizing marker expression. • VlnPlot (shows expression probability distributions across clusters), • and FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring: • RidgePlot, • CellPlot, and • DotPlot as additional methods to view your dataset.

```R
VlnPlot(object = pbmc, features =c("LYZ", "CD14"))

```
![](../assets/images/Part7/plot_7_1.png)

```R
FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")
```
![](../assets/images/Part7/plot_7_2.png)
DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.


```R
library(dplyr)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, features = top10$gene, label = TRUE)
```
![](../assets/images/Part7/plot_7_4.png)
## Assigning cell type identity to clusters

Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types.

```R
current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("T cells","Macrophage/Monocyte", "B cells", "GZMK+ T cells", "NK cells","Neutrophil")
pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
```
![](../assets/images/Part7/plot_7_4.png)