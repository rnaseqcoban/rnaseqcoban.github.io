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
brain <- readRDS(file = "brain_seurat.RDS")
```

# 3. Tính tỷ lệ mitochondria

```R
head(brain@meta.data)
```

```console
                           orig.ident nCount_RNA nFeature_RNA  mouse well        type
A1.B003290.3_38_F.1.1   SeuratProject     390075         3359 3_38_F   A1    Striatum
A1.B003728.3_56_F.1.1   SeuratProject     776439         1718 3_56_F   A1    Striatum
A1.MAA000560.3_10_M.1.1 SeuratProject    1616087         3910 3_10_M   A1      Cortex
A1.MAA000564.3_10_M.1.1 SeuratProject     360004         4352 3_10_M   A1    Striatum
A1.MAA000923.3_9_M.1.1  SeuratProject     290282         2248  3_9_M   A1 Hippocampus
A1.MAA000930.3_8_M.1.1  SeuratProject     574628          948  3_8_M   A1      Cortex
```

```R
mt.genes <- rownames(brain)[grep("^Mt-",rownames(brain))]
C<-GetAssayData(object = brain, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
brain <- AddMetaData(brain, percent.mito, col.name = "percent.mito")
```


```R
rb.genes <- rownames(brain)[grep("^Rp[sl]",rownames(brain))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
brain <- AddMetaData(brain, percent.ribo, col.name = "percent.ribo")
```

```R
VlnPlot(brain, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```

```R
VlnPlot(brain, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
```

```R
VlnPlot(brain, features = "percent.ribo", pt.size = 0.1) + NoLegend()
```

```R
FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```