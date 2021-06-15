---
layout: page
title: 6. Phân nhóm tế bào (Clustering)
permalink: /R_tutorial/clustering/
nav_order: 9
parent: Hướng dẫn phiên bản R
---

# Phân nhóm tế bào (Clustering)

## 1. Giới thiệu

Seurat hiện cung cấp các phương pháp phân nhóm tế bào dựa trên đồ thị. Nói một cách ngắn gọn, các phương pháp này nhúng các tế bào vào đồ thị (graph) - ví dụ đồ thị [K-nearest-neighbor](https://machinelearningcoban.com/2017/01/08/knn/) (KNN), với các dạnh được kết nối với nhau giữa các tế bào có <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> tương tự nhau. Sau đó, thuật toán cố gắng phân vùng graph này thành các `communities` (cộng đồng/nhóm) có tính liên kết với nhau cao. Đầu tiên, Seurat sẽ xây dựng một đồ thị KNN dựa vào [khoảng cách Euclid](https://vi.wikipedia.org/wiki/Kho%E1%BA%A3ng_c%C3%A1ch_Euclid) trong không gian PCA. Sau đó tinh chỉnh số cạnh giữa các vùng lân cận cục bộ của chúng. Để phân nhóm các tế bào, Seurat áp dụng các kỹ thuật tối ưu hóa mô-đun như thuật toán [Louvain](https://python-louvain.readthedocs.io/en/latest/) (mặc định) hoặc [SLM](http://www.ludowaltman.nl/slm/), với mục tiêu tối ưu hóa chức năng mô-đun chuẩn.

Trong Seurat, hàm `FindClusters` có thể thực hiện phân nhóm, và tham số về độ phân giải (resolution) graph nhằm đặt ra mức độ chi tiết các cụm được phân nhóm, với giá trị càng lớn thì số lượng cụm sẽ càng lớn theo. Theo kinh nghiệm phân tích, việc đặt tham số này trong khoảng 0.6 - 1.2 thường trả về kết quả tốt cho các tập dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNAseq</a> khoảng 3000 tế bào. Độ phân giải tối ưu thường tăng với các bộ dữ liệu lớn hơn. Kết quả phân nhóm mới nhất sẽ được lưu trữ trong phần metadata với tên `seurat_clusters`.

Đầu tiên, chúng ta sẽ tính toán K-nearest neighbors và dựng đồ thị SNN (Shared nearest neightbor, một thuật toán khác phát triển từ KNN), sau đó chạy hàm `FindClusters`.

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

## 2. Trực quan hóa sử dụng tSNE

Chúng ta sẽ sử dụng không gian tSNE được tính từ phần trước để trực quan hóa các cụm tế bào được phân nhóm.

```R
DimPlot(object = pbmc, reduction = "tsne")
```
![](../assets/images/Part6/plot_6_1.png)
## 3. Trực quan hóa sử dụng UMAP

Tương tự như vậy với UMAP.

```R
DimPlot(pbmc, reduction = "umap")
```
![](../assets/images/Part6/plot_6_2.png)

--------------------------------------------------------

Nguồn:

[https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html)

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/04-clustering.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/04-clustering.html)
