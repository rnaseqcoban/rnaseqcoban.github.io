---
layout: page
title: 6. Phân nhóm tế bào (Clustering)
permalink: /clustering/
nav_order: 6
---

Source:

[https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html)

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/04-clustering.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/04-clustering.html)

## Cluster the cells

Seurat now includes an graph-based clustering approach compared to (Macosko et al.). Importantly, the distance metric which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data SNN-Cliq, Xu and Su, Bioinformatics, 2015 and CyTOF data PhenoGraph, Levine et al., Cell, 2015. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar gene expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). To cluster the cells, we apply modularity optimization techniques such as the Louvain algorithm (default) or SLM SLM, Blondel et al., Journal of Statistical Mechanics, to iteratively group cells together, with the goal of optimizing the standard modularity function.

The FindClusters function implements the procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. Latest clustering results will be stored in object metadata under seurat_clusters.

First calculate k-nearest neighbors and construct the SNN graph (FindNeighbors), then run FindClusters.

```R
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
```

```console
Computing nearest neighbor graph
Computing SNN
```

```R
pbmc <- FindClusters(pbmc, resolution = 0.4, algorithm = 1)
```
```console
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 1062
Number of edges: 41527

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8916
Number of communities: 6
Elapsed time: 0 seconds
```

## Visualize by tSNE

```R
DimPlot(object = pbmc, reduction = "tsne")
```

## Visualize by UMAP

```R
DimPlot(pbmc, reduction = "umap")
```
