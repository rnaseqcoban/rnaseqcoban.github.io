---
layout: page
title: 4. Tiền xử lý dữ liệu (Preprocessing)
permalink: /preprocessing/
nav_order: 4
---

Source:

7.7-7.10

[https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html](https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html)

[https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/02-normalization.html](https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/02-normalization.html)

https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html


## Normalizing cell library size 

One factor that contributes variation to single-cell RNA-sequencing experiments is called "Library size variation". Library sizes vary for many reasons, including natural differences in cell size, variation of RNA capture, variation in the efficiency of PCR amplification used to generate enough RNA to create the sequencing library. In addition, because scRNA-seq data is often sequenced on highly multiplexed platforms, and the total reads which are derived from each cell may differ substantially.

As a result, while the volume of a cell is informative of a cell's phenotype, there is much more variation in size due to technical factors, and so cells are commonly normalized to have comparable RNA content, becuase this is known to exclude much more technical than biological variation. However, it is important to note that all reasoning about differences between cells after this normalization occurs is restricted to asking question about the relative, not absolute, abundance of RNA in one cell vs another.

Some quantification methods (eg. Cufflinks, RSEM) incorporate library size when determining gene expression estimates and thus do not require this normalization. However, if another quantification method was used then library size must be corrected for.

There are two main approaches to this correction. Many methods use a simple linear scaling to adjust counts such that each cell (row) has about the same total library size. Examples include converting to counts per million (CPM) and closely related methods such as scran. While simple, these approaches do a reasonable job of correcting for differences in library size.

Other methods are more complex, and generally involve parametric modeling of count data to perform nonlinear normalization. These methods are useful when there are more complex sources of unwanted variation (e.g., for highly heterogeneous populations of cells with different sizes).


After removing unwanted genes cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. There have been many methods to normalize the data, but this is the simplest and the most intuitive. The division by total expression is done to change all expression counts to a relative measure, since experience has suggested that technical factors (e.g. capture rate, efficiency of reverse transcription) are largely responsible for the variation in the number of molecules per cell, although genuine biological factors (e.g. cell cycle stage, cell size) also play a smaller, but non-negligible role. The log-transformation is a commonly used transformation that has many desirable properties, such as variance stabilization (can you think of others?).

```R
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
```

```console
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```

A potential drawback of CPM is if your sample contains genes that are both very highly expressed and differentially expressed across the cells. In this case, the total molecules in the cell may depend of whether such genes are on/off in the cell and normalizing by total molecules may hide the differential expression of those genes and/or falsely create differential expression for the remaining genes. One way to mitigate this is to exclude highly expressed genes from the size factor estimation.

## Detection of variable genes across the single cells

Seurat calculates highly variable genes and focuses on these for downstream analysis. FindVariableFeatures calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko et al.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~3,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.


```R
### seurat <- FindVariableFeatures(object = seurat, selection.method = ?, nfeatures = ?) 
### Task: ?FindVariableFeatures into the console to read about different selection methods
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
```
```console
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```


```R
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
```

```console
 [1] "GNLY"   "IGLC2"  "IGLC3"  "S100A9" "FCGR3A" "S100A8" "CDKN1C" "GZMB"   "ITM2C"  "LYZ"
```

```R
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
```

## Scaling the data and removing unwanted sources of variation

Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and clustering.

We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (see example here) and regress this out as well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content.

Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData.


```R
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
```

```console
Regressing out nCount_RNA, percent.mito
  |================================================================================================================| 100%
Centering and scaling data matrix
  |================================================================================================================| 100%
```