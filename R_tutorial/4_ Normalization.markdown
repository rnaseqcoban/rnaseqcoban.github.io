---
layout: page
title: 4. Chuẩn hoá dữ liệu 
permalink: /R_tutorial/Normalizing data/
nav_order: 7
parent: Hướng dẫn phiên bản R
---

# Chuẩn hoá dữ liệu (Normalizing data)

## 1. Chuẩn hóa [library](https://en.wikipedia.org/wiki/CDNA_library) size tế bào

Một yếu tố góp phần tạo ra sự sai khác đối với dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> là sự khác nhau trong kích thước của Library. <a target="_blank" href="https://rnaseqcoban.github.io/def/#library-size" data-tooltip="{{site.data.dict.Library_size}}"  data-tooltip-location="top">Library size</a> khác nhau vì rất nhiều lý do, bao gồm sự khác nhau tự nhiên giữa kích cỡ các tế bào, khả năng bắt giữ được RNA hay hiệu quả khuếch đại PCR khi sử dụng để tạo ra đủ RNA giúp giải trình tự. Ngoài ra, vì dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> thường được giải trình tự bằng các <a target="_blank" href="https://rnaseqcoban.github.io/def/#multiplex-platform" data-tooltip="{{site.data.dict.Multiplex_platform}}"  data-tooltip-location="top">multiplex platform</a> (giải trình tự gộp nhiều mẫu), do đó tổng số reads được lấy ra từ từng tế bào có thể khác nhau đáng kể giữa các mẫu.

Nhìn vào kết quả, chúng ta thấy rằng thể tích của tế bào có thể cung cấp thông tin về kiểu hình của tế bào đó, tuy nhiên sự khác nhau này lại có thể do các yếu tố kỹ thuật. Vì vậy, dữ liệu của <a target="_blank" href="https://rnaseqcoban.github.io/def/#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> thường được chuẩn hóa để có hàm lượng RNA tương đương giữa các tế bào, giúp loại bỏ các lỗi kỹ thuật hơn là thay đổi về mặt sinh học. Một điều quan trọng cần lưu ý đó là tất cả các lý luận về sự khác biệt giữa các tế bào sau khi quá trình chuẩn hóa này xảy ra, sẽ bị hạn chế trong việc đặt câu hỏi về bài toán sinh học. Có nghĩa là mọi so sánh sẽ chỉ mang tính tương đối, không phải tuyệt đối.

Một số phương pháp định lượng đã kết hợp <a target="_blank" href="https://rnaseqcoban.github.io/def/#library-size" data-tooltip="{{site.data.dict.Library_size}}"  data-tooltip-location="top">library size</a> khi ước tính <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> và không cần đến bước chuẩn hóa này (về cơ bản nó đã tích hợp bước này rồi). Còn lại hầu hết các phương pháp khác bắt buộc phải chuẩn hóa cho phù hợp.

Có hai cách tiếp cận chính để chuẩn hóa. Một nhóm các phương pháp sử dụng <a target="_blank" href="https://rnaseqcoban.github.io/def/#scale" data-tooltip="{{site.data.dict.Scale}}"  data-tooltip-location="top">scale</a> tuyến tính đơn giản (ví dụ như đưa về thang 0-1) để điều chỉnh sao cho mỗi tế bào có cùng <a target="_blank" href="https://rnaseqcoban.github.io/def/#library-size" data-tooltip="{{site.data.dict.Library_size}}"  data-tooltip-location="top">library size</a>. Một nhóm khác thì phức tạp hơn, thường liên quan tới mô hình tham số của dữ liệu expression để thực hiện chuẩn hóa phi tuyến tính. Các phương pháp này hữu ích khi có nhiều nguồn variation không mong muốn (ví dụ: sự không đồng nhất mạnh giữa kích thước các tế bào trong 1 mô).

Ở phần này, chúng ta sẽ sử dụng phương pháp chuẩn hóa theo <a target="_blank" href="https://rnaseqcoban.github.io/def/#global-scaling" data-tooltip="{{site.data.dict.Global_scaling}}"  data-tooltip-location="top">global-scaling</a> `LogNormalize` để chuẩn hóa các phép đo <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> cho mỗi tế bào với tổng lượng <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a>, nhân giá trị này với hệ số tỷ lệ (mặc định là 10.000) và biến đổi hàm log cho kết quả. Có rất nhiều phương pháp để chuẩn hóa dữ liệu, nhưng đây là phương pháp đơn giản và trực quan nhất. Phép chia cho tổng <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> được thực hiện để thay đổi tất cả các số lượng <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> thành một số đo tương đối. Kinh nghiệm cho thấy rằng các yếu tố kỹ thuật (ví dụ tỷ lệ bắt giữ được RNA, hiệu quả của phiên mã ngược) chịu trách nhiệm phần lớn cho sự thay đổi số lượng phân tử RNA trên mỗi tế bào. Ngoài ra, các yếu tố sinh học (ví dụ giai đoạn chu kỳ tế bào, kích thước tế bào) cũng đóng một vai trò không nhỏ. Phép biến đổi log là một phép biến đổi thường được sử dụng do có nhiều đặc tính chúng ta mong muốn, chẳng hạn như sự ổn định về phương sai. Chúng ta sẽ sử dụng hàm `NormalizeData` trong Seurat.

Đầu tiên chúng ta sẽ chuẩn hoá dữ liệu bằng phương pháp "LogNormalize" với scale factor là 1e4

```R
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
```

```console
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```

## 2. Phát hiện các gene biến đổi giữa các tế bào

Seurat tính toán các gen có khả năng thay đổi cao và sẽ tập trung vào những gen này để áp dụng các phân tích về sau. `FindVariableFeatures` là hàm tính <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> trung bình và độ phân tán cho mỗi gen, đặt các gen này vào các [bins](https://docs.tibco.com/pub/spotfire/7.0.1/doc/html/bin/bin_what_is_binning.htm), sau đó tính z-score cho sự phân tán trong mỗi bin. Điều này giúp kiểm soát mối quan hệ giữa độ biến thiên và <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> trung bình. Người dùng nên tự đặt các thông số này để biểu đồ phân tán trở nên trực quan hơn. Nhưng để cài đặt thông số chính xác thì phải dựa trên loại dữ liệu thực tế, tính không đồng nhất giữa các mẫu và chiến lược chuẩn hóa. Các tham số ở đây xác định ~ 2,000 gen biến đổi và thường chúng ta sử dụng tham số điển hình cho dữ liệu UMI được chuẩn hóa với tổng số `1e4` phân tử.

Sau đó chúng ta tìm 2000 gene biến đổi nhiều nhất bằng hàm "FindVariableFeatures"


```R
### seurat <- FindVariableFeatures(object = seurat, selection.method = ?, nfeatures = ?) 
### Task: ?FindVariableFeatures into the console to read about different selection methods
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
```
```console
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```

 Chúng ta có thể kiểm tra tốp 10 gene biến đổi nhiều nhất trong số 2000 gene

```R
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
```

```console
 [1] "GNLY"   "IGLC2"  "IGLC3"  "S100A9" "FCGR3A" "S100A8" "CDKN1C" "GZMB"   "ITM2C"  "LYZ"
```
Chúng ta có thể biểu diễn 2000 gene này bằng đồ thị điểm 

```R
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
```
![](../../assets/images/Part4/plot.png)

## 3. Scale dữ liệu và loại bỏ các nguồn biến đổi/biến thể không mong muốn

Tập dữ liệu tế bào của bạn có thể chứa các nguồn biến đổi/biến thể không mong muốn. Không chỉ bao gồm các yếu tố gây nhiễu về mặt kỹ thuật, mà còn cả [batch effect](https://en.wikipedia.org/wiki/Batch_effect) (hiệu ứng lô/mẫu), hoặc thậm chí là các nguồn biến đổi sinh học (các giai đoạn chu kỳ tế bào khác nhau). Hồi quy những tín hiệu này ra khỏi phân tích có thể cải thiện việc giảm số chiều và phân nhóm sau này. Để giảm thiểu ảnh hưởng của những tín hiệu này, Seurat xây dựng các mô hình tuyến tính để dự đoán <a target="_blank" href="https://rnaseqcoban.github.io/def/#gene-expressionbiểu-hiện-gene" data-tooltip="{{site.data.dict.Gene_expression}}"  data-tooltip-location="top">gene expression</a> dựa trên các biến do người dùng xác định. [Z-score](https://www.investopedia.com/ask/answers/021115/what-difference-between-standard-deviation-and-z-score.asp) (thang chuẩn hóa) residual của các mô hình này được lưu trữ trong `scale.data` và được sử dụng để giảm số chiều và phân nhóm.

Cuối cùng chúng ta scale dữ liệu và loại bỏ biến đổi gây ra do sự khác biệt về gene và phần trăm mitochondrial gene trong từng tế 

```R
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
```

```console
Regressing out nCount_RNA, percent.mito
  |================================================================================================================| 100%
Centering and scaling data matrix
  |================================================================================================================| 100%
```

----------------------------------------------------

Nguồn:

[https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html](https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html)

[https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/02-normalization.html](https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/02-normalization.html)

[https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html](https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html)
