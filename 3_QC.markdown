---
layout: page
title: 3. Kiểm soát chất lượng dữ liệu (Quality control)
permalink: /qc/
nav_order: 3
---

Source:

7.1 - 7.6

[https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html](https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html)

## 1. Giới thiệu

Once we have our expression matrix, it should be examined to remove poor quality cells which were not detected in the initial processing of the raw reads. Failure to remove low quality cells at this stage may add technical noise which has the potential to obscure the biological signals of interest in the downstream analysis.

Since there is currently no standard method for performing scRNAseq, the expected values for the various QC measures that will be presented here can vary substantially from experiment to experiment. Thus, to perform QC we will be looking for cells which are outliers with respect to the rest of the dataset rather than comparing to independent quality standards. Consequently, care should be taken when comparing quality metrics across datasets collected using different protocols.

## 2. Đọc dữ liệu

```R
pbmc <- readRDS(file = "PBMC_1k.RDS")
```

## 3. Filter out low-quality cells

The Seurat object initialization step above only considered cells that expressed at least 350 genes. Additionally, we would like to exclude cells that are damaged. A common metric to judge this (although by no means the only one) is the relative expression of mitochondrially derived genes. When the cells apoptose due to stress, their mitochondria becomes leaky and there is widespread RNA degradation. Thus a relative enrichment of mitochondrially derived genes can be a tell-tale sign of cell stress. Here, we compute the proportion of transcripts that are of mitochondrial origin for every cell (percent.mt), and visualize its distribution as a violin plot. We also use the FeatureScatter function to observe how percent.mt correlates with other metrics.

```R
head(brain@meta.data)
```

```console
                      orig.ident nCount_RNA nFeature_RNA
AAACCCAAGGAGAGTA-1 SeuratProject       8288         2620
AAACGCTTCAGCCCAG-1 SeuratProject       5512         1808
AAAGAACAGACGACTG-1 SeuratProject       4283         1562
AAAGAACCAATGGCAG-1 SeuratProject       2754         1225
AAAGAACGTCTGCAAT-1 SeuratProject       6592         1831
AAAGGATAGTAGACAT-1 SeuratProject       8845         2048
```

Here we calculated the percent mitochondrial reads and added it to the Seurat object in the slot named meta.data. This allowed us to plot using the violin plot function provided by Seurat.

```R
# Tìm mitochondiral genes bắt đầu bằng MT-
mt.genes <- rownames(pbmc)[grep("^MT-",rownames(pbmc))]
# Lấy count matrix từ Seurat
C<-GetAssayData(object = pbmc, slot = "counts")

# Tính phần trăm mt genes
percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
# Đưa kết quả % mt genes vào bảng metadata
pbmc <- AddMetaData(pbmc, percent.mito, col.name = "percent.mito")
```


```R
rb.genes <- rownames(pbmc)[grep("^RP[SL]",rownames(pbmc))]
percent.ribo <-  Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
pbmc <- AddMetaData(pbmc, percent.ribo, col.name = "percent.ribo")
```

```R
VlnPlot(pbmc, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```
![](../assets/images/Part3/plot_3_1.png)
```R
VlnPlot(pbmc, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
```
![](../assets/images/Part3/plot_3_2.png)
```R
VlnPlot(pbmc, features = "percent.mito", pt.size = 0.1) + NoLegend()
```
![](../assets/images/Part3/plot_3_3.png)
```R
VlnPlot(pbmc, features = "percent.ribo", pt.size = 0.1) + NoLegend()
```
![](../assets/images/Part3/plot_3_4.png)
```R
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```
![](../assets/images/Part3/plot_3_5.png)
```R
FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mito")
```
![](../assets/images/Part3/plot_3_6.png)
```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mito < 25)
```

```R
VlnPlot(pbmc, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```
![](../assets/images/Part3/plot_3_7.png)


