---
layout: page
title: 7. Phân tích khác biệt biểu hiện gene (Differential expression analysis)
permalink: /R_tutorial/dea/
nav_order: 10
parent: Hướng dẫn phiên bản R
---

# Phân tích khác biệt biểu hiện gen

## 1. Giới thiệu

Sau khi phân tế bào thành các nhóm khác nhau, câu hỏi đặt ra là điều gì đã làm cho mỗi nhóm khác biệt với các nhóm còn lại trong tập dữ liệu, và làm sao để chú thích chính xác các nhóm theo phân loại tế bào (tế bào máu, T-cells, B-cells,...) của chúng.

Có nhiều cách để thực hiện bước này:

- Tìm kiếm các marker gene (gene có tính đại diện, nổi bật) làm tăng/giảm điều hòa (upregulation/downregulation) ở các loại tế bào đang quan tâm (so với phần còn lại của bộ dữ liệu)
- So sách toàn bộ biểu hiện gene giữa các nhóm tế bào với nhau
- Sử dụng các phương pháp tự động so sánh các tế bào trên cơ sở dữ liệu về biểu hiện gene của các loại tế bào đã được công bố giúp dễ dàng hơn trong việc chú giải loại tế bào.

Các phương pháp tự động là một tiến bộ đầy hứa hẹn, nhưng vẫn chưa thể thay đổi được sự cẩn thận của con người.

Đối với các loại tế bào được xác định rõ ràng, chúng ta hy vọng rằng các marker genes sẽ cho thấy rõ sự khác biệt trong biểu hiện của gene giữa các loại tế bào với toàn bộ phần còn lại của dữ liệu. Nó cho phép chúng ta có thể sử dụng những phương pháp đơn giản và hiệu quả. Chúng ta sẽ nhắm vào các phương pháp đó trong tutorial này.

## 2. So sánh phân bố

Các thuật toán phân tích khác biệt biểu hiện gene cung cấp các cách tiếp cận khác nhau để so sánh sự phân bố của biểu hiện gene giữa một nhóm tế bào này với một nhóm tế bào khác. Không giống như [Bulk RNA-seq](https://rpubs.com/ewilkinson_KRISP/479615), chúng ta thường có một số lượng lớn các mẫu (tế bào) cho mỗi nhóm mà chúng ta đang so sánh trong các thí nghiệm <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNAseq</a>. Do đó, chúng ta có thể tận dụng toàn bộ phân phối của các giá trị trong mỗi nhóm để xác định sự khác biệt giữa các nhóm thay vì chỉ so sánh các ước tính của biểu hiện gene trung bình như tiêu chuẩn thường gặp của Bulk RNASeq.

Seurat có thể giúp bạn tìm các gene markers mà nó giúp xác định được phân loại của các nhóm tế bào thông qua sự khác biệt trong biểu hiện gene. Mặc định, đầu tiên thuật toán giúp xác định sự tăng và giảm trong biểu hiện gene của một nhóm tế bào so với tất cả các tế bào khác. Hàm `FindAllMarkers` giúp tự động hóa quy trình này cho tất cả các nhóm, nhưng bạn cũng có thể kiểm tra các nhóm so với nhau hoặc so với toàn bộ tế bào khi sử dụng hàm `FindMarkers`. 

Tham số `min.pct` yêu cầu một gene được phát hiện ở một tỷ lệ phần trăm tối thiểu ở một trong hai nhóm tế bào đang so sánh và tham số `thresh.test` yêu cầu một gene phải được biểu hiện khác biệt (trung bình) theo một số lượng nhất định giữa hai nhóm. Bạn có thể đặt cả hai tham số này là 0, nhưng thời gian chạy hàm sẽ tăng lên đáng kể bởi vì điều này sẽ làm Seurat kiểm tra thêm một số lượng lớn các gen không có khả năng phân biệt cao. Ngoài ra, là một tùy chọn khác để tăng tốc thuật toán này, bạn có thể sử dụng tham số `max.cells.per.ident`. Tham số này giúp giảm số lượng mẫu (tế bào) cho mỗi nhóm. Mỗi nhóm này sẽ không được có số lượng tế bào nhiều hơn giới hạn này. Mặc dù nói chung sẽ có sự mất mát về sức mạnh của thuật toán, nhưng tốc độ tăng vẫn có thể sẽ giúp được đưa các gene có biểu hiện khác biệt nhiều nhất lên trên cùng và hoàn toàn có ý nghĩa thống kê.


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

Seurat có một số tùy chọn cho việc chọn công cụ để phân tích khác biệt biểu hiện gene, nó có thể được xác định bằng tham số `test.use`. Ví dụ: kiểm tra [ROC](https://rstudio-pubs-static.s3.amazonaws.com/267441_5459af9d83ae44f18a13aea4a479f31f.html) trả về 'classification power' (độ mạnh của kết quả phân loại) cho bất kỳ gene marker riêng lẻ nào (từ 0 - ngẫu nhiên, đến 1 - hoàn hảo).

```R
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

Seurat đưa vào khá nhiều các công cụ khác nhau để trực quan hóa kết quả.

- `Vlplot` là hàm dùng để vẽ [Violin plot](https://en.wikipedia.org/wiki/Violin_plot) giúp hiển thị phân bố xác suất biểu hiện gene giữa các nhóm.

- `FeaturePlot` là hàm dùng để hiển thị biểu hiện gene trên biểu đồ tSNE hoặc PCA. Hàm này thường được dùng khá phổ biến vì tính trực quan của nó.

- `CellPlot`, `DotPlot` hay `RidgePlot` cũng là những phương pháp bổ sung để quan sát dữ liệu của bạn.

```R
VlnPlot(object = pbmc, features =c("LYZ", "CD14"))

```
![](../assets/images/Part7/plot_7_1.png)

```R
FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")
```
![](../assets/images/Part7/plot_7_2.png)

Hàm `DoHeatmap` tạo ra một biểu đồ nhiệt cho các tế bào và genes nhất định. Trong trường hợp này, chúng ta sẽ vẽ biểu đồ nhiệt với 20 gene markers hàng đầu (hoặc tất cả các markers nếu số lượng ít hơn 20) cho mỗi nhóm.

```R
library(dplyr)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, features = top10$gene, label = TRUE)
```
![](../assets/images/Part7/plot_7_3.png)

## 3. Gán tên cho các loại tế bào

Thông thường đây là bước khá là gặp khó khăn khi thực hiện, nhưng may mắn đối với loại dữ liệu này, các gene marker biểu hiện rất rõ ràng và đại diện cho các loại tế bào có trong máu.

```R
current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("T cells","Macrophage/Monocyte", "B cells", "GZMK+ T cells", "NK cells","Neutrophil")
pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
```
![](../assets/images/Part7/plot_7_4.png)

Tới đây, bạn đã hoàn thành việc phân tích cơ bản một tập dữ liệu scRNA-seq. Bạn có thể áp dụng cấu trúc phân tích này với các dữ liệu scRNA-seq nhưng cần chỉnh sửa các tham số sao cho phù hợp với bản chất dữ liệu và câu hỏi sinh học bạn đang nhắm tới.

-------------------------------------------------

Nguồn: 

[https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#dechapter](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#dechapter)

[https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/05-diffexp.html](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/05-diffexp.html)
