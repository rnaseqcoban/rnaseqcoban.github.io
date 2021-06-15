---
layout: page
title: 1. Giới thiệu về single cell RNA sequencing (scRNA-seq)
permalink: /R_tutorial/intro/
nav_order: 4
parent: Hướng dẫn phiên bản R
---

# Giới thiệu về single cell RNA-seq

## 1. Đặc điểm chính của Bulk RNA-seq và single cell RNA-seq

### a. Bulk RNA-seq

-	Đo **biểu hiện trung bình của một gene** trong một tập hợp tế bào 
-	Hữu dụng trong so sánh hệ phiên mã và định lượng gene biểu hiện đặc hiệu (Ví dụ: so sánh hệ phiên mã của mô gan từ chuột lành và chuột mang bệnh để đánh giá mức độ tăng giảm biểu hiện của gene)
-	Không phù hợp trong nghiên cứu những hệ thống mang tính không đồng nhất (Ví dụ: so sánh sự thay đổi của cơ quan trong quá trình phát triển, hoặc so sánh các mẫu mô có cấu trúc và thành phần phức tạp như mô não)

### b. Single cell RNA-seq (scRNA-seq)

-	Là một công nghệ mới, được đề cập tới đầu tiên trong bài báo [Tang et al. 2009](https://www.nature.com/articles/nmeth.1315)
- Trở nên phổ biến vào khoảng năm 2014 sau khi có nhiều <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#phương-pháp-scrna-seq" data-tooltip="{{site.data.dict.Phuong_phap_scRNA_seq}}"  data-tooltip-location="top">phương pháp scRNA-seq</a> mới ra đời và chi phí giải trình tự giảm xuống
-	Có rất nhiều protocol cho <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> được sử dụng, ví dụ SMART-seq2, CELL-seq, Drop-seq,...
-	Có nhiều hệ thống máy đã được thương mại hoá của <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> như Fluidigm C1, Wafergen ICELL8, và 10X Genomics Chromium,...
-	Nhiều phương pháp phân tích tính toán từ <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#bulk-rna-seq" data-tooltip="{{site.data.dict.Bulk_RNA_seq}}"  data-tooltip-location="top">Bulk RNA-seq</a> có thể được sử dụng cho <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a>
-	Trong phần lớn các trường hợp, phương pháp phân tích tính toán yêu cầu sự thay đổi tương ứng hoặc phát triển những phương pháp mới 
-	Đo **phân bố biểu hiện của từng gene** của từng tế bào một tập hợp tế bào đơn lẻ. 
- Mở ra những câu hỏi sinh học mới trong đó cho phép xác định sự thay đổi nào trong hệ phiên mã của một tế bào là quan trọng nhất(Ví dụ: đó có thể là loại tế bào, sự đáp ứng không đồng nhất của tế bào, tính ngẫu nhiên của gene biểu hiện, hệ quả của điều hoà hoạt động gene...)

## 2. Các bước cơ bản trong thí nghiệm scRNA-seq

1. Phân lập tế bào từ mô 

2. Phân tách RNA từ từng tế bào 

3. Tổng hợp <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#cdna" data-tooltip="{{site.data.dict.cDNA}}" data-tooltip-location="top">cDNA</a>, nhân lên <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#cdna" data-tooltip="{{site.data.dict.cDNA}}" data-tooltip-location="top">cDNA</a> bằng phản ứng <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#pcr" data-tooltip="{{site.data.dict.PCR}}" data-tooltip-location="top">PCR</a> hoặc nhân lên <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#cdna" data-tooltip="{{site.data.dict.cDNA}}" data-tooltip-location="top">cDNA</a> thông qua phản ứng phiên mã trong ống nghiệm (*in vitro* transcription) kết hợp với phiên mã ngược 

4. Tạo thư viện DNA giải trình tự và giải trình tự 

5. Phân tích kết quả 

![RNA-Seq_workflow-5 pdf](https://user-images.githubusercontent.com/59919924/111622576-78713600-883d-11eb-9b16-fcfdcbdb8972.jpg)



## 4. Các bước cơ bản trong phân tích dữ liệu scRNA-seq
Phân tích dữ liệu scRNAseq thường bắt đầu với "<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a>". Trong <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a>, mỗi hàng biểu diễn một gene và mỗi cột biểu diễn một tế bào. Như vậy, mỗi ô biểu diễn mức độ biểu hiện của một gene trong một tế bào. 
Các bước xây dựng một "<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a>" bao gồm: 

1. **Read Quality control**: Đọc và kiểm tra chất lượng đoạn giải trình tự. (Công cụ: FastQC hay Karen).
![image](https://user-images.githubusercontent.com/59919924/112940337-bc870380-9178-11eb-8cb6-21636601b279.png)
2. **Alignment**: Liên kết đoạn giải trình tự với hệ genome tham khảo (reference genome) để tìm gene tương ứng với đoạn giải trình tự. (Công cụ sử dụng: STAR, TopHat)  
![image](https://user-images.githubusercontent.com/59919924/112940242-9b261780-9178-11eb-976c-bb6ec2d56b13.png)
3. **Mapping Quality control**: Kiểm tra độ chính xác của bước liên kết đoạn giải trình tự với reference genome dựa trên các chỉ số như: tỉ lệ giữa <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#ribosome-rna" data-tooltip="{{site.data.dict.Ribosome_RNA}}"  data-tooltip-location="top">ribosome RNA</a> và <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#transposon-rna" data-tooltip="{{site.data.dict.Transposon_RNA}}"  data-tooltip-location="top">transposon RNA</a>, tỉ lệ các đoạn nối đặc hiệu, độ sâu (<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#read-depth" data-tooltip="{{site.data.dict.Read_depth}}"  data-tooltip-location="top">read depth</a>), đoạn đọc trải dài qua các điểm nối (<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#splice-junctionsmối-nối" data-tooltip="{{site.data.dict.Splice_junctions}}"  data-tooltip-location="top">splice junctions</a>).  
![image](https://user-images.githubusercontent.com/59919924/112940214-906b8280-9178-11eb-9ef0-24cddc201055.png)
4. **Reads Quantification**: Tính toán mức độ biểu hiện của mỗi gene trong một tế bào. Với <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#phương-pháp-scrna-seq" data-tooltip="{{site.data.dict.Phuong_phap_scRNA_seq}}"  data-tooltip-location="top">phương pháp scRNA-seq</a> sử dụng kỹ thuật "<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#tag-base" data-tooltip="{{site.data.dict.Tag_base}}"  data-tooltip-location="top">tag-base</a>", <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#umi" data-tooltip="{{site.data.dict.UMI}}"  data-tooltip-location="top">UMI</a> có thể được dùng để tính số lượng tuyệt đối của một phân tử RNA. 

Sau khi được xây dựng, "<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#count-matrixexpression-matrix" data-tooltip="{{site.data.dict.Count_matrix}}"  data-tooltip-location="top">Expression matrix</a>" được dùng trong phân tích tiếp theo thông qua các nền tảng thiết kế riêng cho <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a>. 

Hiện nay, có rất nhiều nền tảng dùng để phân tích dữ liệu <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a>. Tuy nhiên trong tutorials này chúng ta sẽ chủ yếu sử dụng các nền tảng sau:

1. Scanpy trong Python

2. Seurat trong R

3. Bioconductor: nguồn mở, chứa các phần mềm dùng cho phân tích dữ liệu genomics trong R

## 5. Thách thức

Sự khác biệt lớn nhất giữa phân tích <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#bulk-rna-seq" data-tooltip="{{site.data.dict.Bulk_RNA_seq}}"  data-tooltip-location="top">Bulk RNA-seq</a> và <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#scrna-seq" data-tooltip="{{site.data.dict.ScRNA_seq}}"  data-tooltip-location="top">scRNA-seq</a> là mỗi thư viện giải trình tự là đặc trưng của một tế bào thay vì của một quần thể tế bào. Vì vậy chúng ta cần lưu ý khi so sánh giữa các tế bào khác nhau. Những yếu tố có thể dẫn đến sự không tương quan khi so sánh bao gồm:

1. Sự nhân lên thiếu đồng đều của <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#cdna" data-tooltip="{{site.data.dict.cDNA}}" data-tooltip-location="top">cDNA</a> (Bias amplification)

2. Gene “bị loại bỏ” (Gene dropout) là trường hợp gene biểu hiện ở mức độ trung bình hoặc yếu trong một tế bào, nhưng không xuất hiện trong tế bào khác. 

Nguyên nhân của việc này là do số lượng đầu vào thấp của RNA trong một tế bào. Làm thế nào để cải thiện hiệu suất và giảm sự không đồng đều vẫn đang được tiếp tục nghiên cứu

## 6. Phương pháp thí nghiệm 

![image](https://user-images.githubusercontent.com/59919924/111624962-6f359880-8840-11eb-9be9-b05088a116e3.png)


Có rất nhiều phương pháp thí nghiệm, tuy nhiên, các phương pháp này có thể được phân nhóm dựa trên hai đặc điểm sau: (1) Cách định lượng RNA (2) Cách “bắt giữ” từng tế bào

Có hai cách định lượng chính. Một là giải trình tự toàn bộ đoạn RNA (full-length). Hai là dùng “đầu dò” gắn đầu 3’ hoặc 5’, và chỉ giải trình tự đoạn gắn “đầu dò” (<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#tag-base" data-tooltip="{{site.data.dict.Tag_base}}"  data-tooltip-location="top">tag-based</a>).  Nên sử dụng phương pháp nào phụ thuộc vào mục đích tạo dữ liệu, vì cả hai đều có ưu và nhược điểm riêng. Về lí thuyết, "full-length" sẽ bao phổ (<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#coverageđộ-phủ" data-tooltip="{{site.data.dict.Coverage}}"  data-tooltip-location="top">coverage</a>) đồng đều đoạn phiên mã RNA. Tuy nhiên, thực tế cho thấy sự bao phổ (<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#coverageđộ-phủ" data-tooltip="{{site.data.dict.Coverage}}"  data-tooltip-location="top">coverage</a>) là không đồng đều. Nghĩa là có phần của đoạn phiên mã được bao phổ nhiều hơn so với các phần khác. Ưu điểm lớn nhất của phương pháp “<a target="_blank" href="https://rnaseqcoban.github.io/R/def//#tag-base" data-tooltip="{{site.data.dict.Tag_base}}"  data-tooltip-location="top">tag-based</a>” là có thêm trình tự nhận dạng phân tử đặc hiệu (Unique Molecular Identifier - <a target="_blank" href="https://rnaseqcoban.github.io/R/def//#umi" data-tooltip="{{site.data.dict.UMI}}"  data-tooltip-location="top">UMI</a>) trong đoạn mồi. Trình tự này là đặc trưng cho từng loại RNA khác nhau, vậy nên có thể được sử dụng trong định lượng. Tuy nhiên, nhược điểm của phương pháp này là bị giới hạn ở một đầu của RNA, và gây khó khăn trong việc phân biệt các isoform.

Khi nghĩ về phương pháp "bắt giữ" từng tế bào, chúng ta có thể xem xét ba sự lựa chọn sau: microwell, microfluidic, và droplet:

- **Microwell:** sử dụng các phương pháp phân tách như FACS, pippette, laser capturing, để đưa từng tế bào vào trong một giếng có kích thước micro. Ưu điểm của microwell là có thể chụp hình ảnh của tế bào trong một giếng. Nhược điểm là quy mô nhỏ và khối lượng công việc lớn.

- **Microfulidic (Fluidigm C1):** sử dụng nguyên lý thuỷ động lực học để dẫn dắt từng tế bào vào một buồng phản ứng. Ưu điểm là tính tự động hoá và quy mô cao hơn so với microwell. Tuy nhiên, hiệu suất thấp, chỉ có 10% số lượng tế bào được giữ lại ở các buồng phản ứng; và giá thành cho một đĩa “chip” quá đắt đỏ. 

- **Droplet (10X Genomics, in-Drop):** ý tưởng của phương pháp droplet là bao bọc một tế bào trong một giọt dầu đã có sẵn <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#bead" data-tooltip="{{site.data.dict.Bead}}" data-tooltip-location="top">bead</a> gắn đoạn mồi và nguyên liệu cho phản ứng tạo <a target="_blank" href="https://rnaseqcoban.github.io/R/def/#cdna" data-tooltip="{{site.data.dict.cDNA}}" data-tooltip-location="top">cDNA</a>. Đây là phương pháp có hiệu suất và quy mô cao nhất trong ba loại. 


-----------------------------------------------

Nguồn:

[https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html](https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html)
