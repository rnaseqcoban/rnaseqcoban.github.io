---
layout: page
title: 5. Giảm chiều dữ liệu (Dimensionality reduction)
permalink: /reducedim/
nav_order: 5
---

Source:

7.3.2

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/03-dimensionality-reduction.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/03-dimensionality-reduction.html)

# Giảm chiều dữ liệu (Dimensionality reduction)

## Phân tích thành phần chính (Principal component analysis - PCA)

Các phương pháp giảm kích thước tìm cách lấy một tập hợp lớn dữ liệu và trả về một tập hợp nhỏ hơn với các thành phần vẫn chứa hầu hết thông tin trong tập dữ liệu gốc.

Một trong những hình thức giảm kích thước đơn giản nhất là PCA. Phân tích thành phần chính (PCA) là một phương pháp toán học giúp biến đổi một số biến số tương quan (ví dụ: gene expression) thành một số (nhỏ hơn) các biến không tương quan được gọi là thành phần chính ("PC").

Về mặt toán học, các PC tương ứng với các vector riêng (eigenvector) của ma trận hiệp phương sai. Các eigenvector được sắp xếp theo trị riêng (eigenvalue) để thành phần chính đầu tiên chiếm càng nhiều sự thay đổi (variability) trong dữ liệu càng tốt và mỗi thành phần tiếp theo lần lượt có phương sai cao nhất có thể với điều kiện là nó trực giao với các thành phần trước đó.

![](../assets/images/Part5/plot_5_1.png)

Về mặt sinh học, phương pháp giảm kích thước này hữu ích và thích hợp với dữ liệu scRNAseq vì tế bào phản ứng với môi trường của chúng bằng cách kích hoạt các chương trình điều hòa dẫn đến sự biểu hiện của mô-đun/gen. Kết quả là, sự biểu hiện gen hiển thị sự đồng biểu hiện có cấu trúc và giảm số chiều bằng cách  phân tích các gene đồng biến đổi đó thành các thành phần chính, được sắp xếp theo mức độ phương sai mà chúng giải thích được.

Khuyến khích các bạn đọc thêm bài này: [Link - Machine learning cơ bản](https://machinelearningcoban.com/2017/06/15/pca/)

## Thực hiện giảm số chiều tuyến tính

Tham khảo phiên bản Seurat v3: các feature (cách gọi thay thế cho gene khi đưa vào mô hình toán học) có tính biến số cao (highly variable - có nghĩa là các gene này có độ biến đổi cao, biểu hiện thấp ở một số tế bào này, nhưng biểu hiện cao ở các tế bào khác) được phát hiện thông qua hàm `FindVariableFeatures` đã được thực hiện ở phần tiền xử lý dữ liệu. Chúng ta sẽ sử dụng các feature này khi đã được scale để thực hiện PCA.

```R
pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)
```

Seurat v3 cung cấp các chức năng để trực quan hóa PCA: 

- Biểu đồ PCA
- Biểu đồ PCA được tô màu bởi một đối tượng định lượng (i.e gene) 
- Biểu đồ phân tán (scatter)  
- Biểu đồ Violin và Ridge 
- Biểu đồ nhiệt

```R
# Examine and visualize PCA results a few different ways
DimPlot(object = pbmc, reduction = "pca")
```
![](../assets/images/Part5/plot_5_2.png)
```R
# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")
```
![](../assets/images/Part5/plot_5_3.png)

Đặc biệt, `DimHeatmap` cho phép dễ dàng khám phá các nguyên nhân chính của sự không đồng nhất trong tập dữ liệu và nó rất hữu ích khi quyết định PC nào nên đưa vào cho các phân tích tiếp theo. Cả tế bào và gene đều được sắp xếp theo điểm số PCA của chúng. PCA là một phân tích rất có giá trị để khám phá các nhóm gen tương quan với nhau.

## Xác định các thành phần chính có ý nghĩa thống kê

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
![](../assets/images/Part5/plot_5_5.png)

A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with ElbowPlot. In this example, it looks like the elbow would fall around PC 5.

```R
ElbowPlot(object = pbmc)
```
![](../assets/images/Part5/plot_5_6.png)
## Run Non-linear dimensional reduction (tSNE)

An alternative to PCA for visualizing scRNASeq data is a tSNE plot. tSNE (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) with random walks on the nearest-neighbour network to map high dimensional data (i.e. our 18,585 dimensional expression matrix) to a 2-dimensional space. In contrast with PCA, tSNE can capture nonlinear structure in the data, and tries to preserve the local distances between cells. Due to the non-linear and stochastic nature of the algorithm, tSNE is more difficult to intuitively interpret: while tSNE faithfully represents local relationships, it doesn't always capture the relatioships between more distant cells correctly.

tSNE is a stochastic algorithm which means running the method multiple times on the same dataset will result in different plots. To ensure reproducibility, we fix the "seed" of the random-number generator in the code below so that we always get the same plot.

```R
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = pbmc, reduction = "tsne")
```
![](../assets/images/Part5/plot_5_7.png)
## Run UMAP

UMAP (Uniform Approximation and Projection) is another nonlinear dimensionality reduction method. Like tSNE, UMAP is nondeterministic and requires that we fix the random seed to ensure reproducibility. While tSNE optimizes for local structure, UMAP tries to balance the preservation of local and global structure. For this reason, we prefer UMAP over tSNE for exploratory analysis and general visualization.

```R
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
DimPlot(pbmc, reduction = "umap")
```
![](../assets/images/Part5/plot_5_8.png)