---
layout: page
title: 2. Dữ liệu scRNAseq
permalink: /tabula/
nav_order: 2
---

Source:

[https://scrnaseq-course.cog.sanger.ac.uk/website/tabula-muris.html](https://scrnaseq-course.cog.sanger.ac.uk/website/tabula-muris.html)


# 1. Dữ liệu Tabula Muris

Trong toàn bộ phần hướng dẫn này, chúng ta sẽ sử dụng dữ liệu scRNAseq từ dự án Tabula Muris. Đây là dự án hợp tác quốc tế nhằm tìm ra thông tin sinh học phân tử của tất của các loại tế bào trên chuột ở cấp độ phiên mã (transcriptomics level). Họ kết hợp giữa giải trình tự thông lượng cao, có độ phủ thấp 10X (tên platform) với giải trình tự thông lượng thấp, độ phủ cao FACS-sorted cells Smartseq2 (tên của platform). Ở blog này, chúng ta sẽ chỉ sử dụng dữ liệu Smartseq2.

Dữ liệu này lần đầu tiên công bố vào tháng 12 năm 2017, bao gồm khoảng 100,000 tế bào thu thập được từ các mô/cơ quan khác nhau của chuột. Các bạn có thể dựa vào code trong phần hướng dẫn để thay đổi loại mô/cơ quan mà các bạn muốn phân tích.

# 2. Tải dữ liệu scRNAseq

Dữ liệu Tabula Muris từ não chuột có thể download được từ đây: [Link Download](https://github.com/chanzuckerberg/scRNA-python-workshop/blob/master/content/data.zip).

Sau khi tải về, bạn giải nén và giữ lại 2 files: `brain_counts.csv` và `brain_metadata.csv`.

# 3. Đọc dữ liệu

Chúng ta có thể đọc dữ liệu ma trận count (count có nghĩa là số lượng reads chồng lên nhau ở mỗi gene khi giải trình tự) từ file .csv. Sau đó quan sát thử dataframe (bảng):
```R
dat <- read.delim("brain_counts.csv", sep=",", header=TRUE, row.names = 1)
dat[1:5,1:5]
```

```console
                        X0610005C13Rik X0610007C21Rik X0610007L01Rik X0610007N19Rik X0610007P08Rik
A1.B003290.3_38_F.1.1                0            125             16              0              0
A1.B003728.3_56_F.1.1                0              0              0              0              0
A1.MAA000560.3_10_M.1.1              0            348              0              0              0
A1.MAA000564.3_10_M.1.1              0             41             36              0              0
A1.MAA000923.3_9_M.1.1               0             53              0              0              0
```

Chúng ta có thể kiểm tra số dòng và số cột của dữ liệu:

```R
dim(dat)
```

```console
[1]  3401 23434
```
Bây giờ, chúng ta có thể tách một số metadata từ tên cột:

```R
# Lấy tất cả các cell id ở tên dòng
cellIDs <- rownames(dat)
# Tách chuỗi bằng dấu chấm
cell_info <- strsplit(cellIDs, "\\.")
# Lấy dữ liệu về batch (phần tử đầu tiên của mỗi mảng)
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
# Lấy dữ liệu về plate
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
# Lấy dữ liệu về chuột
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))
```
Chúng ta có thể kiểm tra phân bố chuột:

```R
summary(factor(Mouse))
```
```console
3_10_M 3_11_M 3_38_F 3_39_F 3_56_F  3_8_M  3_9_M 
   980    253    355    241    111    590    87
```

Cuối cùng chúng ta đọc file metadata của dữ liệu.

```R
# Đọc file metadata
ann <- read.table("brain_metadata.csv", sep=",", header=TRUE)
# Lấy ra những dòng có cột đầu tiên match với cell id
ann <- ann[match(cellIDs, ann[,1]),]
# Tách cell type ở cột thứ 3
celltype <- ann[,3]
```

# 4. Tạo object
Để tính toán hiệu quả hơn thì chúng ta cần đưa các dữ liệu vào các loại object khác nhau và tính toán toàn bộ xung quanh object đó. Ở đây chúng ta sẽ thử 2 loại object khác nhau là SingleCellExperiment và Seurat object (xem thêm ở phần setup).

## 4.1. SingleCellExperiment
Để tạo ra SingleCellExperiment object, chúng ta cần đưa tất cả thông tin chú giải tế bào vào 1 dataframe duy nhất.

```R
library("SingleCellExperiment")
# Tạo dataframe cho metadata
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
# Thế tên sample vào metadata
rownames(cell_anns) <- rownames(dat)
# Đưa tất cả thông tin có được bao gồm count matrix và metadata vào SingleCellExperiment object
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(t(dat))), colData=cell_anns)
```

## 4.2. Seurat
Seurat là một trong những package phổ biến nhất hiện tại trong xử lý dữ liệu scRNAseq.

```R
library("Seurat")
# Tạo dataframe cho metadata
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
# Thế tên sample vào metadata
rownames(cell_anns) <- rownames(dat)
# Tương tự như SCE object
brain <- CreateSeuratObject(counts = as.matrix(t(dat)), meta.data = cell_anns)
brain
```

```console
An object of class Seurat 
23433 features across 3401 samples within 1 assay 
Active assay: RNA (23433 features, 0 variable features)
```

Chúng ta có thể lưu lại Seurat object bằng cách:

```R
saveRDS(brain, file = "brain_seurat.RDS")
```

Sau đó, bạn có thể đọc lại Seurat object bất cứ lúc nào:

```R
brain <- readRDS(file = "brain_seurat.RDS")
```

Từ các tutorial sau, chúng ta sẽ chỉ sử dụng Seurat object để phân tích dữ liệu.