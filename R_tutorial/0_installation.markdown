---
layout: page
title: Hướng dẫn cài đặt thư viện
permalink: /R_tutorial/setup/
nav_order: 2
parent: Hướng dẫn phiên bản R
---

# Hướng dẫn cài đặt

## Cài đặt R

R là một ngôn ngữ lập trình có thế mạnh về tính toán thống kê và đồ họa. Bạn có thể cài đặt R theo các bước sau:

- Bước 1: Truy cập [https://cloud.r-project.org/](https://cloud.r-project.org/). Click vào "Download R for ..." (Linux, OS X hoặc Windows). Sau đó click vào "base". Rồi click vào "Download R ..."

![](../assets/images/install/1.png)

- Bước 2: Cài đặt R trên máy tính của bạn.

![](../assets/images/install/2.png)

## Cài đặt Rstudio

RStudio là một chương trình bao gồm các tiện ích và chức năng giúp người dùng dễ dàng sử dụng ngôn ngữ lập trình R hơn. Bạn có thể cài đặt RStudio theo các bước sau:

- Bước 1: Truy cập [https://www.rstudio.com/](https://www.rstudio.com/). Click vào bản RStudio phù hợp với hệ điều hành máy tính bạn đang sử dụng.

![](../assets/images/install/3.png)

- Bước 2: Cài đặt RStudio trên máy tính của bạn.

![](../assets/images/install/4.png)

## Cài đặt thư viện

Trong khóa học, một số thư viện về xử lý dữ liệu và scRNA-seq cần phải cài đặt:

- `Seurat`: Là một trong những công cụ xử lý dữ liệu scRNA-seq mạnh nhất bây giờ. Bạn có thể xem tại trang chủ [https://satijalab.org/seurat/index.html](https://satijalab.org/seurat/index.html). Cách cài đặt như sau:

```R
# Nhập câu lệnh sau đây trong R (hoặc RStudio, nếu đã được cài đặt)
install.packages('Seurat')
library(Seurat)
```

Nếu bạn gặp cảnh báo dưới đây, thì hãy nhập `y` và enter:

```console
package which is only available in source form, and may need compilation of C/C++/Fortran: 'Seurat'
Do you want to attempt to install these from sources?
y/n:
```

- `SingleCellExperiment`: thư viện này là tùy chọn, nếu bạn muốn thực hành lưu trữ dữ liệu scRNAseq trong SingleCellExperiment thì mới cần cài đặt. Thư viện này cung cấp class dùng để lưu trữ dữ liệu scRNA-seq.

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)
```
