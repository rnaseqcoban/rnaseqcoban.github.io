---
layout: page
title: 5. Giảm chiều dữ liệu (Dimensionality reduction)
permalink: /reducedim/
nav_order: 5
---

Source:

7.3.2

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/03-dimensionality-reduction.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/03-dimensionality-reduction.html)

## Principle components analysis

Dimensionality reduction methods seek to take a large set of variables and return a smaller set of components that still contain most of the information in the original dataset.

One of the simplest forms of dimensionality reduction is PCA. Principal component analysis (PCA) is a mathematical procedure that transforms a number of possibly correlated (e.g., expression of genes in a network) variables into a (smaller) number of uncorrelated variables called principal components ("PCs").

Mathematically, the PCs correspond to the eigenvectors of the covariance matrix. The eigenvectors are sorted by eigenvalue so that the first principal component accounts for as much of the variability in the data as possible, and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to the preceding components (the figure below is taken from here).

Biologically, this type of dimensionality reduction is useful and appropriate because cells respond to their environment by turning on regulatory programs that result in expression of modules of genes. As a result, gene expression displays structured co-expression, and dimnesionality reduction by principle component analysis groups those co-varying genes into principle components, ordered by how much variation they explain.


Image[]

Now that we have a clean expression matrix, we can use PCA to visualize an overview of the data and assess confounding factors. SCANPY provides several very useful functions to simplify visualisation.

## Perform linear dimensional reduction

refered to Seurat v3 (latest): high variable features are accessed through the function HVFInfo(object). Despite RunPCA has a features argument where to specify the features to compute PCA on, I’ve been modifying its values and the output PCA graph has always the same dimensions, indicating that the provided genes in the features argument are not exactly the ones used to compute PCA. Wether the function gets the HVG directly or does not take them into account, I don’t know.

```R
pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)
```

Seurat v3 provides functions for visualizing: - PCA - PCA plot coloured by a quantitative feature - Scatter plot across single cells - Scatter plot across individual features - Variable Feature Plot - Violin and Ridge plots - Heatmaps

```R
# Examine and visualize PCA results a few different ways
DimPlot(object = pbmc, reduction = "pca")
```

```R
# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")
```

In particular DimHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.

## Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

In Macosko et al, we implemented a resampling test inspired by the jackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of gene scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value genes.

```R
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)

pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20, reduction = "pca")

JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca")
```


A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with ElbowPlot. In this example, it looks like the elbow would fall around PC 5.

```R
ElbowPlot(object = pbmc)
```

## Run Non-linear dimensional reduction (tSNE)

An alternative to PCA for visualizing scRNASeq data is a tSNE plot. tSNE (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) with random walks on the nearest-neighbour network to map high dimensional data (i.e. our 18,585 dimensional expression matrix) to a 2-dimensional space. In contrast with PCA, tSNE can capture nonlinear structure in the data, and tries to preserve the local distances between cells. Due to the non-linear and stochastic nature of the algorithm, tSNE is more difficult to intuitively interpret: while tSNE faithfully represents local relationships, it doesn't always capture the relatioships between more distant cells correctly.

tSNE is a stochastic algorithm which means running the method multiple times on the same dataset will result in different plots. To ensure reproducibility, we fix the "seed" of the random-number generator in the code below so that we always get the same plot.

```R
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = pbmc, reduction = "tsne")
```

## Run UMAP

UMAP (Uniform Approximation and Projection) is another nonlinear dimensionality reduction method. Like tSNE, UMAP is nondeterministic and requires that we fix the random seed to ensure reproducibility. While tSNE optimizes for local structure, UMAP tries to balance the preservation of local and global structure. For this reason, we prefer UMAP over tSNE for exploratory analysis and general visualization.

```R
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
DimPlot(pbmc, reduction = "umap")
```