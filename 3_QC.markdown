---
layout: page
title: 3. Kiểm soát chất lượng dữ liệu (Quality control)
permalink: /qc/
nav_order: 3
---

Source:

7.1 - 7.6

[https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html](https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html)

# Kiểm soát chất lượng dữ liệu (Quality control)

## 1. Giới thiệu

Một khi đã có dữ liệu, chúng ta cần kiểm tra và loại bỏ các tế bào có dữ liệu kém. Sai sót trong việc loại bỏ các tế bào có dữ liệu kém có thể làm tăng nhiễu do kĩ thuât và ảnh hưởng tới việc tìm ra các thông tin mang ý nghĩa sinh học.
 
Hiện tại chưa có quy trình chuẩn cho thí nghiệm scRNAseq, vậy nên các chỉ số kiểm tra chất lượng có thể biến đổi giữa các thí nghiệm. Khi kiểm tra chất lượng dữ liệu, chúng ta sẽ tìm những tế bào nằm ngoài vùng phân bố của phần lớn các tế bào còn lại. Lưu ý, khi so sánh chất lượng dữ liệu từ các quy trình khác nhau, chúng ta cần phải cân nhắc cẩn trọng.

## 2. Đọc dữ liệu

```R
pbmc <- readRDS(file = "PBMC_1k.RDS")
```

## 3. Filter out low-quality cells

Bước khởi đầu trong kiểm tra chất lượng Seurat Object là giữ lại các tế bào có tối thiểu 350 genes. Sau đó, dựa vào sự biểu hiện tương đối của các gene từ ti thể (mitochodrial gene) để loại bỏ các tế bào không còn nguyên vẹn. Khi tế bào bước vào chu trình chết tự nhiên gây ra bởi “stress”, ti thể bên trong tế bào này trở nên rò rỉ và RNA bị phá hỏng. Ở trong ví dụ này, chúng ta sẽ tính toán phần trăm các mitochrondial gene trên mỗi tế bào, thêm vào bảng dữ liệu, và biểu diễn chỉ số này bằng biểu đồ violin.

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


