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

Nếu bạn muốn sử dụng trình duyệt để download dữ liệu Tabula Muris: [Link Download](https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells/5715040). Sau khi download về, bạn giải nén và đổi tên ...

Nếu bạn sử dụng terminal:

```bash
wget https://ndownloader.figshare.com/files/10038307
unzip -o 10038307
wget https://ndownloader.figshare.com/files/10038310
mv 10038310 FACS_metadata.csv
wget https://ndownloader.figshare.com/files/10039267
mv 10039267 FACS_annotations.csv
```

Bạn có thể đếm số dòng của file:
```bash
wc -l droplet_annotation.csv
```
```console
## 54838 droplet_annotation.csv
```

**Câu hỏi:** Có bao nhiêu tế bào có chú giải từ FACS?

**Đáp án:** FACS: 54,838 tế bào.


# 3. Đọc dữ liệu

Chúng ta có thể đọc dữ liệu ma trận count (count có nghĩa là số lượng reads chồng lên nhau ở mỗi gene khi giải trình tự) từ file .csv. Sau đó quan sát thử dataframe (bảng):
```R
dat = read.delim("FACS/Kidney-counts.csv", sep=",", header=TRUE)
dat[1:5,1:5]
```

```console
##               X A14.MAA000545.3_8_M.1.1 E1.MAA000545.3_8_M.1.1
## 1 0610005C13Rik                       0                      0
## 2 0610007C21Rik                       1                      0
## 3 0610007L01Rik                       0                      0
## 4 0610007N19Rik                       0                      0
## 5 0610007P08Rik                       0                      0
##   M4.MAA000545.3_8_M.1.1 O21.MAA000545.3_8_M.1.1
## 1                      0                       0
## 2                      0                       0
## 3                      0                       0
## 4                      0                       0
## 5                      0                       0
```

Chúng ta có thể thấy cột đầu tiên trong bảng là tên của gene, vì vậy chúng ta cần chuyển cột đầu tiên thành tên dòng (row names) để có được ma trận số của count:

```R
dim(dat)
```

```console
## [1] 23433   866
```
```R
rownames(dat) <- dat[,1]
dat <- dat[,-1]
```

Lưu ý dữ liệu Smartseq2 có thể chứa spike-ins (ký hiệu bắt đầu bằng ERCC-):
```R
 rownames(dat)[grep("^ERCC-", rownames(dat))]
```
```console
##  [1] "ERCC-00002" "ERCC-00003" "ERCC-00004" "ERCC-00009" "ERCC-00012"
##  [6] "ERCC-00013" "ERCC-00014" "ERCC-00016" "ERCC-00017" "ERCC-00019"
## [11] "ERCC-00022" "ERCC-00024" "ERCC-00025" "ERCC-00028" "ERCC-00031"
## [16] "ERCC-00033" "ERCC-00034" "ERCC-00035" "ERCC-00039" "ERCC-00040"
## [21] "ERCC-00041" "ERCC-00042" "ERCC-00043" "ERCC-00044" "ERCC-00046"
## [26] "ERCC-00048" "ERCC-00051" "ERCC-00053" "ERCC-00054" "ERCC-00057"
## [31] "ERCC-00058" "ERCC-00059" "ERCC-00060" "ERCC-00061" "ERCC-00062"
## [36] "ERCC-00067" "ERCC-00069" "ERCC-00071" "ERCC-00073" "ERCC-00074"
## [41] "ERCC-00075" "ERCC-00076" "ERCC-00077" "ERCC-00078" "ERCC-00079"
## [46] "ERCC-00081" "ERCC-00083" "ERCC-00084" "ERCC-00085" "ERCC-00086"
## [51] "ERCC-00092" "ERCC-00095" "ERCC-00096" "ERCC-00097" "ERCC-00098"
## [56] "ERCC-00099" "ERCC-00104" "ERCC-00108" "ERCC-00109" "ERCC-00111"
## [61] "ERCC-00112" "ERCC-00113" "ERCC-00116" "ERCC-00117" "ERCC-00120"
## [66] "ERCC-00123" "ERCC-00126" "ERCC-00130" "ERCC-00131" "ERCC-00134"
## [71] "ERCC-00136" "ERCC-00137" "ERCC-00138" "ERCC-00142" "ERCC-00143"
## [76] "ERCC-00144" "ERCC-00145" "ERCC-00147" "ERCC-00148" "ERCC-00150"
## [81] "ERCC-00154" "ERCC-00156" "ERCC-00157" "ERCC-00158" "ERCC-00160"
## [86] "ERCC-00162" "ERCC-00163" "ERCC-00164" "ERCC-00165" "ERCC-00168"
## [91] "ERCC-00170" "ERCC-00171"
```

Bây giờ, chúng ta có thể tách một số metadata từ tên cột:

```R
cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))
```
Chúng ta có thể kiểm tra phân bố chuột:

```R
summary(factor(Mouse))
```
```console
## 3_10_M 3_11_M 3_38_F 3_39_F  3_8_M  3_9_M 
##    104    196    237    169     82     77
```

Sau đó kiểm tra có lỗi kỹ thuật nào không:

```R
table(Mouse, Plate)
```

```console
##         Plate
## Mouse    B001717 B002775 MAA000545 MAA000752 MAA000801 MAA000922
##   3_10_M       0       0         0       104         0         0
##   3_11_M       0       0         0         0       196         0
##   3_38_F     237       0         0         0         0         0
##   3_39_F       0     169         0         0         0         0
##   3_8_M        0       0        82         0         0         0
##   3_9_M        0       0         0         0         0        77
```

Cuối cùng chúng ta đọc chú giải về loại tế bào từ file annotation và hợp (match) với dữ liệu ma trận biểu hiện gene.

```R
ann <- read.table("FACS_annotations.csv", sep=",", header=TRUE)
ann <- ann[match(cellIDs, ann[,1]),]
celltype <- ann[,3]
```

# 4. Tạo object
Để tính toán hiệu quả hơn thì chúng ta cần đưa các dữ liệu vào các loại object khác nhau và tính toán toàn bộ xung quanh object đó. Ở đây chúng ta sẽ thử 2 loại object khác nhau là SingleCellExperiment và Seurat object (xem thêm ở phần setup).

## 4.1. SingleCellExperiment
Để tạo ra SingleCellExperiment object, chúng ta cần đưa tất cả thông tin chú giải tế bào vào 1 dataframe duy nhất.
```R
library("SingleCellExperiment")
library("scater")
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
rownames(cell_anns) <- colnames(dat)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(dat)), colData=cell_anns)
```


Cuối cùng, chúng ta cần đánh dấu spike-ins:

```R
isSpike(sceset, "ERCC") <- grepl("ERCC-", rownames(sceset))
```
## 4.2. Seurat
