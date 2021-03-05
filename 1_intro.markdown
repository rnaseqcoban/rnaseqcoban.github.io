---
layout: page
title: 1. Giới thiệu về single-cell RNA sequencing (scRNAseq)
permalink: /intro/
nav_order: 1
---

Source:

[https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html](https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html)

# Giới thiệu về single-cell RNA-seq

## 1. Bulk RNA-seq

-	Đo mức độ biểu hiện trung bình của một gene trong một tập hợp tế bào 
-	Hữu dụng trong so sánh hệ phiên mã (transcriptomic), ví dụ cùng một loại mô của các loài khác nhau
-	Hữu dụng trong định lượng gene đặc hiệu biểu hiện, ví dụ trong nghiên cứu về bệnh
-	Không phù hợp trong nghiên cứu những hệ thống không có tính đồng nhất, vì sụ nghiên cứu về sinh học phát triển, mẫu mô có cấu trúc và thành phần phức tạp như mô não

## 2. scRNAseq

-	Là một công nghệ mới, được đề cập tới đầu tiên trong bài báo của Tang et al. 2009
-	Chỉ bắt đầu trở nên phổ biến vào khoảng năm 2014 khi có nhiều protocols mới ra đời và sự giảm gía thành cho việc giải trình tự 
-	Đo phân bố biểu hiện của một gene trong một tập hợp tế bào 
-	Trả lời câu hỏi sự thay đổi nào trong hệ phiên mã của một tế bào thực sự là quan trọng. Đó có thể là loại tế bào, sự đáp ứng không đồng nhất của tế bào, tính ngẫu nhiên của gene biểu hiện, hệ quả của điều hoà hoạt động gene. 
-	Có rất nhiều protocol cho scRNAseq được sử dụng, ví dụ SMART-seq2, CELL-seq, Drop-seq
-	Nền tảng thương mại hoá của scRNAseq gồm có Fluidigm C1, Wafergen ICELL8, và 10X Genomocs Chromium 
-	Nhiều phương pháp phân tích tính toán từ bulk RNA-seq có thể được sử dụng 
-	Trong phần lớn các trường lợp, phương pháp phân tích tính toán yêu cầu sự thay đổi tương ứng hoặc là phát triển những phương pháp mới. 

## 3. Các bước cơ bản trong thí nghiệm scRNAseq

(1)	Phân lập tế bào từ mô 

(2)	Phân tách RNA từ từng tế bào 

(3)	Tổng hợp cDNA

(4)	Nhân lên cDNA bằng phản ứng PCR hoặc nhân lên cDNA thông qua phản ứng phiên mã trong ống nghiệm (in vitro trascription) kết hợp với phiên mã ngược 

(5)	Tạo thư viện DNA giải trình tự 

(6)	Giải trình tự 

(7)	Phân tích kết quả 

## 4. Các bước cơ bản trong phân tích dữ liệu scRNAseq


(1)	Read QC

(2)	Alignment

(3)	Mapping QC

(4)	Cell QC

(5)	Normalization 

(6)	Cofounders

Hiện nay, có rất nhiều nền tảng dùng để phân tích dữ liệu scRNAseq. Tuy nhiên trong Tutorials này chúng ta sẽ chủ yếu sử dụng các nền tảng sau:

(1)	Scanpy

(2)	Seurat: dùng trong R

(3)	Bioconductor: nguồn mở, chứa các phần mềm dùng cho phân tích dữ liệu genomics

## 5. Thách thức

Sự khác biệt lớn nhất giữa phân tích bulk RNAseq và scRNAseq là mỗi thư viện giải trình tự là đặc trưng của một tế bào thay vì của một quần thể tế bào. Vì vậy chúng ta cần lưu ý khi so sánh giữa các tế bào khác nhau. Những yếu tố có thể dẫn đến sự không tương quan khi so sánh bao gồm:

(1)	Sự nhân lên thiếu đồng đều của cDNA (Amplification)

(2)	Gene “bị loại bỏ” (Gene dropout) là trường hợp gene biểu hiện ở mức độ trung bình hoặc yếu trong một tế bào, nhưng không xuất hiện trong tế bào khác. 

Nguyên nhân của việc này là do số lượng đầu vào thấp  của RNA trong một tế bào. Làm thế nào để cải thiện hiệu suất và giảm sự không đồng đều vẫn đang được tiếp tục nghiên cứu

## 6. Phương pháp thí nghiệm 

Có rất nhiều phương pháp thí nghiệm, tuy nhiên, các phương pháp này có thể được phân nhóm dựa trên hai đặc điểm sau: (1) Cách định lượng RNA (2) Cách “bắt dữ” từng tế bào

Có hai cách định lượng chính, một là giải trình tự toàn bộ đoạn RNA (full-length), hai là dùng “đầu dò” gắn đầu 3’ hoặc 5’, và chỉ giải trình tự đoạn gắn “đầu dò” (tag-based).  Nên sử dụng phương pháp nào phụ thuộc vào mục đích sử dụng dữ liệu vì cả hai đều có ưu và nhược điểm riêng. Về lí thuyết, “full-length” sẽ bao phổ (coverage) đồng đều đoạn phiên mã RNA. Tuy nhiên, thực tế cho thấy sự bao phổ (coverage) là không đồng đều. Nghĩa là có phần của đoạn phiên mã được bao phổ nhiều hơn so với các đoạn khác. Ưu điểm lớn nhất của phương pháp “tag-based” là có thêm trình tự nhận dạng đặc hiệu (Unique Molecular Identifier -UMI) trong đoạn mồi. Trình tự này là đặc trưng cho từng loại RNA khác nhau, vậy nên có thể được sử dụng trong định lượng. Tuy nhiên, nhược điểm của phương pháp này là bị giới hạn ở một đầu của RNA, và gây khó khăn trong việc phân biệt các isoform

Khi nghĩ về phương pháp phân tách từng tế bào, chúng ta có thể xem xét ba sự lựa chọn sau: microwell, microfluidic, và droplet:

- **Microwell:** sử dụng các phương pháp chọn lọc và phân tách, ví dụ như FACS, pippette, laser capturing, để đưa từng tế bào vào trong một giếng có kích thước micro. Ưu điểm là có thể chụp hình ảnh của tế bào trong một giếng. Nhược điểm là quy mô nhỏ và khối lượng công việc lớn.

- **Microfulidic (Fluidigm C1):** sử dụng nguyên lý thuỷ động lực học để đua từng tế bào vào một buồng phản ứng. Ưu điểm là tính tự động hoá và quy mô cao hơn so với Microwell. Tuy nhiên, hiệu suất thấp, chỉ có 10% số lượng tế bào được giữ lại ở các buồng phản ứng, và giá thành cho một đĩa “chip” quá đắt đỏ. 

- **Droplet (10X Genomics, in-Drop):** ý tưởng của phương pháp droplet là bao bọc một tế bào trong một giọt dầu đã có sẵn “bead” gắn đoạn mồi và nguyên liệu cho phản ứng tạo cDNA. Đây là phương pháp có hiệu suất và quy mô cao nhất trong ba loại. 
