---
layout: page
title: 3. Kiểm soát chất lượng dữ liệu (Quality control)
permalink: /qc/
nav_order: 3
---

Source:

7.1 - 7.6

[https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html](https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html)

# 1. Giới thiệu

# 2. Dữ liệu

```R
pbmc <- readRDS(file = "PBMC_1k.RDS")
```

# 3. Tính tỷ lệ mitochondria

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
![](images/part3/plot_3_1.png)
```R
VlnPlot(pbmc, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
```
![](images/part3/plot_3_2.png)
```R
VlnPlot(pbmc, features = "percent.mito", pt.size = 0.1) + NoLegend()
```
![](images/part3/plot_3_3.png)
```R
VlnPlot(pbmc, features = "percent.ribo", pt.size = 0.1) + NoLegend()
```
![](images/part3/plot_3_4.png)
```R
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```
![](images/part3/plot_3_5.png)
```R
FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mito")
```
![](images/part3/plot_3_6.png)
```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mito < 25)
```

```R
VlnPlot(pbmc, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```
![](images/part3/plot_3_7.png)


