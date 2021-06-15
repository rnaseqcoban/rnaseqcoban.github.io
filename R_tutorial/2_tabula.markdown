---
layout: page
title: 2. Dữ liệu scRNA-seq
permalink: /R_tutorial/tabula/
nav_order: 5
parent: Hướng dẫn phiên bản R
---

# Dữ liệu scRNA-seq

Để hiểu chi tiết về cách đọc dữ liệu, phần đầu tiên bạn sẽ sử dụng dữ liệu Tabula Muris từ file count matrix. Rồi sau đó sẽ tạo 2 loại object là SingleCellExperiment và Seurat. Ở phần thứ 2, bạn sẽ sử dụng cách đọc ngắn gọn hơn với dữ liệu PMBC từ 10X Genomics theo định dạng .H5. Và từ những tutorial sau cũng sẽ chỉ làm việc với dữ liệu này.

# Phần 1: Tabula Muris

## 1. Dữ liệu Tabula Muris

Trong phần hướng dẫn này, chúng ta sẽ sử dụng dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> từ dự án Tabula Muris. Đây là dự án hợp tác quốc tế nhằm tìm ra thông tin sinh học phân tử của tất của các loại tế bào trên chuột ở cấp độ phiên mã (transcriptomics level). Họ kết hợp giữa giải trình tự <a target="_blank" href="https://rnaseqcoban.github.io/def/#thông-lượng" data-tooltip="{{site.data.dict.Thong_luong}}"  data-tooltip-location="top">thông lượng</a> cao, có <a target="_blank" href="https://rnaseqcoban.github.io/def/#coverageđộ-phủ" data-tooltip="{{site.data.dict.Coverage}}"  data-tooltip-location="top">độ phủ</a> thấp 10X (tên platform) với giải trình tự <a target="_blank" href="https://rnaseqcoban.github.io/def/#thông-lượng" data-tooltip="{{site.data.dict.Thong_luong}}"  data-tooltip-location="top">thông lượng</a> thấp, <a target="_blank" href="https://rnaseqcoban.github.io/def/#coverageđộ-phủ" data-tooltip="{{site.data.dict.Coverage}}"  data-tooltip-location="top">độ phủ</a> cao FACS-sorted cells Smartseq2 (tên của protocol). Ở blog này, chúng ta sẽ chỉ sử dụng dữ liệu Smartseq2.

Dữ liệu này lần đầu tiên công bố vào tháng 12 năm 2017, bao gồm khoảng 100,000 tế bào thu thập được từ các mô/cơ quan khác nhau của chuột. Các bạn có thể dựa vào code trong phần hướng dẫn để thay đổi loại mô/cơ quan mà các bạn muốn phân tích.

## 2. Tải dữ liệu scRNA-seq

Dữ liệu Tabula Muris từ não chuột có thể download được từ đây: [Link Download](https://github.com/chanzuckerberg/scRNA-python-workshop/blob/master/content/data.zip).

Sau khi tải về, bạn giải nén và giữ lại 2 files: `brain_counts.csv` và `brain_metadata.csv`.

## 3. Đọc dữ liệu

Chúng ta có thể đọc dữ liệu từ count matrix từ file .csv. Sau đó quan sát thử <a target="_blank" href="https://rnaseqcoban.github.io/def/#dataframe" data-tooltip="{{site.data.dict.Dataframe}}"  data-tooltip-location="top">dataframe</a> (bảng):
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
Bây giờ, chúng ta có thể tách một số <a target="_blank" href="https://rnaseqcoban.github.io/def/#metadata" data-tooltip="{{site.data.dict.Metadata}}"  data-tooltip-location="top">metadata</a> từ tên cột:

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

## 4. Tạo object
Để tính toán hiệu quả hơn thì chúng ta cần đưa các dữ liệu vào các loại object khác nhau và tính toán toàn bộ xung quanh object đó. Ở đây chúng ta sẽ thử 2 loại object khác nhau là SingleCellExperiment và <a target="_blank" href="https://rnaseqcoban.github.io/def/#seurat-object" data-tooltip="{{site.data.dict.Seurat_Object}}"  data-tooltip-location="top">Seurat object</a> (xem thêm ở phần setup).

### 4.1. SingleCellExperiment
Để tạo ra SingleCellExperiment object, chúng ta cần đưa tất cả thông tin chú giải tế bào vào 1 <a target="_blank" href="https://rnaseqcoban.github.io/def/#dataframe" data-tooltip="{{site.data.dict.Dataframe}}"  data-tooltip-location="top">dataframe</a> duy nhất.

```R
library("SingleCellExperiment")
# Tạo dataframe cho metadata
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
# Thế tên sample vào metadata
rownames(cell_anns) <- rownames(dat)
# Đưa tất cả thông tin có được bao gồm count matrix và metadata vào SingleCellExperiment object
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(t(dat))), colData=cell_anns)
```

### 4.2. Seurat
Seurat là một trong những package phổ biến nhất hiện tại trong xử lý dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNAseq</a>.

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

Chúng ta có thể lưu lại <a target="_blank" href="https://rnaseqcoban.github.io/def/#seurat-object" data-tooltip="{{site.data.dict.Seurat_Object}}"  data-tooltip-location="top">Seurat object</a> bằng cách:

```R
saveRDS(brain, file = "brain_seurat.RDS")
```

Sau đó, bạn có thể đọc lại <a target="_blank" href="https://rnaseqcoban.github.io/def/#seurat-object" data-tooltip="{{site.data.dict.Seurat_Object}}"  data-tooltip-location="top">Seurat object</a> bất cứ lúc nào:

```R
brain <- readRDS(file = "brain_seurat.RDS")
```

# Phần 2: PBMC

## 1. Dữ liệu Peripheral blood mononuclear cells (PBMC)

Dữ liệu PBMCs là dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> miễn phí từ công ty 10X Genomics. PBMSs là các tế bào đơn nhân trong máu ngoại vi. Trong bài này chúng ta sẽ sử dụng mẫu từ những người hiến tặng khỏe mạnh.

## 2. Tải dữ liệu scRNA-seq

Sẽ có 2 cách để tải như sau:


- **Cách 1**: Bạn có thể tải miễn phí tại [Link](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3). Để tải thì cần cung cấp một số thông tin cho họ. Sau khi điền xong, bạn sẽ có thể tải file `Feature / cell matrix HDF5 (filtered)`.

- **Cách 2**:

Sử dụng Command promt `cmd` nếu bạn sử dụng hệ điều hành Windows hoặc terminal nếu dùng MacOS:

```console
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
```

- Sử dụng terminal trong Linux

```console
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
```

## 3. Đọc dữ liệu

Seurat cung cấp 1 hàm để đọc dữ liệu từ file `.h5` là `Read10X_h5`.

```R
v3 <- Read10X_h5("./pbmc_1k_v3_filtered_feature_bc_matrix.h5")
```

`v3` lúc này sẽ là <a target="_blank" href="https://rnaseqcoban.github.io/def/#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a> ở dạng <a target="_blank" href="https://rnaseqcoban.github.io/def/#sparse-matrix" data-tooltip="{{site.data.dict.Sparse_matrix}}"  data-tooltip-location="top">sparse matrix</a> giúp tối ưu việc lưu trữ dữ liệu. Về cơ bản cũng là một <a target="_blank" href="https://rnaseqcoban.github.io/def/#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a> để chúng ta có thể đưa vào <a target="_blank" href="https://rnaseqcoban.github.io/def/#seurat-object" data-tooltip="{{site.data.dict.Seurat_Object}}"  data-tooltip-location="top">Seurat object</a>.

## 4. Tạo object

Đưa dữ liệu PBMCs và tạo <a target="_blank" href="https://rnaseqcoban.github.io/def/#seurat-object" data-tooltip="{{site.data.dict.Seurat_Object}}"  data-tooltip-location="top">Seurat object</a>.

```R
# Tạo Seurat object
pbmc <- CreateSeuratObject(counts = v3, project = "pbmc_1k_v3")
# Lưu lại Seurat object
saveRDS(pbmc, file = "PBMC_1k.RDS")
```

Từ những tutorial sau, chúng ta sẽ sử dụng Seurat object của PBMC để phân tích.

-----------------------------------------

Nguồn:

[https://scrnaseq-course.cog.sanger.ac.uk/website/tabula-muris.html](https://scrnaseq-course.cog.sanger.ac.uk/website/tabula-muris.html)
