---
layout: page
title: Từ điển thuật ngữ
permalink: /def/
nav_order: 3
---

# Thuật ngữ

### Bulk RNA-seq

Nếu đặt lên bàn cân so sánh với scRNA-seq thì bulk RNA-seq Là công nghệ giải trình tự cũ hơn, nó không đạt được đến độ phân giải từng tế bào mà chỉ đo được thông tin transcriptomics trung bình của cả một mẫu sinh học (trung bình của một nhóm các tế bào). Mỗi điểm dữ liệu nhận được tương đương với 1 mẫu sinh học (trung bình tập hợp rất nhiều tế bào).

### ScRNA-seq

scRNA-seq là công nghệ giải trình tự thông lượng cao với khả năng đo được biểu hiện gene với độ phân giải ở mức tế bào. Mỗi điểm dữ liệu thu nhận được tương đương với 1 tế bào.

### Phương pháp scRNA-seq

Phương thức thực hiện giải trình tự. Mỗi phương thức sẽ khác nhau về giá cả, thông lượng, và phẩm chất đầu ra.

### UMI

Viết tắt của Unique molecular identifier, nghĩa là định danh phân tử duy nhất. Nó giúp xác định được phân tử phiên mã bắt nguồn từ read nào. UMI giúp thu gọn lại các bản sao từ PCR.

### Count matrix/Expression matrix. 

Count ở đây chỉ số lượng UMI. Có thể hiểu rằng count định lượng cho biểu hiện gene. Count matrix là một ma trận với cột là tế bào, dòng là gene. Giá trị ở điểm bất kỳ trong ma trận là lượng UMI của một gene ở 1 tế bào.

### Metadata

Các dữ liệu đi kèm như thông tin của tế bào, platform giải trình tự,...

### Tag-base

Là kỹ thuật giải trình tự chỉ sự dụng một đoạn ngắn (tag) ở một vị trí nhất định để xác định phân tử RNA được giải trình tự.

### Thông lượng

Là khả năng tự động hóa thí nghiệm cho việc giải trình tự lặp lại ở quy mô số lượng tế bào lớn trở nên khả thi. Có thể hiểu đơn giản là quy mô số tế bào được giải trình tự.

### Coverage/Độ phủ

Mô tả số lượng read trung bình phủ lên bộ gene tham chiếu. Mức độ bao phủ = (tổng số base được tạo ra) / (kích thước của bộ gene tham chiếu)

### Sparse matrix

Ma trận thưa là ma trận chứa rất ít phần tử khác 0.

### library size

RNA đã được giải trình tự gọi là thư viện RNA. Library size có nghĩa là tổng số reads được giải trình tự.

### Dataframe

Là một bảng hoặc một cấu trúc mảng hai chiều, trong đó mỗi cột chứa các giá trị của một biến số và mỗi hàng chứa một bộ giá trị tương ứng từ mỗi cột.

### cDNA

Viết tắt của Complementary DNA, nghĩa là DNA bổ sung. Nó được tổng hợp từ khuôn mẫu RNA sợi đơn trong một phản ứng được xúc tác bởi enzyme phiên mã ngược.

### PCR

Phản ứng chuỗi polymerase giúp khuếch đại bản sao DNA trong quá trình giải trình tự.

### Single cell RT-qPCR

Single cell reverse transcription quantitative real-time PCR là một kỹ thuật giúp khuếch đại và định lượng mức độ biểu hiện gene với thông lượng thấp.

### Ribosome RNA

Phân tử RNA trong tế bào tạo thành một phần của bào quan tổng hợp protein được gọi là ribosome và được xuất ra tế bào chất để giúp dịch thông tin trong RNA thông tin (mRNA) thành protein.

### Transposon RNA

Là sản phẩm của các trình tự DNA có thể thay đổi vị trí trong bộ gene.

### Read depth

Độ sâu giải trình tự. Mô tả số lần một nucleotide nhất định trong bộ gene được đọc trong khi giải trình tự. Độ sâu càng lớn thì độ tin cậy của kết quả giải trình tự càng cao.

### splice junctions/mối nối

Trong quá trình nối RNA (RNA splicing), mối nối là vị trí của một intron cũ trong mRNA trưởng thành.

### Bead

Bead là các hạt nhựa silica hoặc hạt từ phủ silica sử dụng muối chaotropic để phá vỡ liên kết hydro và liên kết với axit nucleic, tạo điều kiện cho giúp rửa trôi các chất ô nhiễm.

### Seurat Object

Là trung tâm của mỗi phân tích scRNA-seq khi sử dụng Seurat. Nó lưu trữ tất cả thông tin được liên kết với dữ liệu scRNA-seq, bao gồm count matrix, chú thích, phân tích, v.v. Tất cả những gì cần thiết để xây dựng một Seurat Object là count matrix (hàng là gene, cột là tế bào).

### Multiplex platform

Là các platform có khả năng giải trình tự gộp nhiều mẫu cùng một lúc.

### Gene expression/biểu hiện gene

Biểu hiện gene ở đây có thể được định lượng bằng UMI count.

### Scale

Quy mô giá trị hoặc kích thước.

### Global-scaling

Thay đổi quy mô giá trị toàn cục.

### Bins

Cách chia các nhóm dữ liệu thành từng bin/thùng có chiều rộng bằng nhau.

### Batch effect

Trong sinh học phân tử, hiệu ứng lô (batch) xảy ra khi các yếu tố không phải sinh học gây ra những thay đổi trong dữ liệu của một thí nghiệm.

### Confounding effect

Một số hiệu ứng gây nhiễu phổ biến như: read nhiễu, hiện tượng drop-out,... là hậu quả của read depth thấp, có tính chất ngẫu nhiên.

### Z-score residual

Là một thuật ngữ phổ biến cho phần dư chuẩn hóa. Để tính toán phần dư chuẩn hóa của một tập dữ liệu, giá trị trung bình và độ lệch chuẩn của giá trị dữ liệu phải được ước tính.

### Feature

Đặc tính, hay thuật ngữ dùng để thay thế từ gene trong phân tích dữ liệu scRNA-seq.

### Resampling

Lẫy mẫu lại

### jackStraw procedure

Phương pháp jackstraw cho phép xác định một tập hợp các gene được liên kết với bất kỳ PC nhất định nào, hoặc với một tập hợp con của các PC, hoặc sự kết hợp tuyến tính của hai hoặc nhiều PC,...

### Null distribution

Trong kiểm định giả thuyết thống kê, phân phối rỗng là phân phối xác suất của thống kê kiểm định khi giả thuyết rỗng là đúng.

### cutoff

Giá trị dùng để xác định điểm cắt.

### tSNE

t-distributed stochastic neighbor embedding là một phương pháp thống kê để trực quan hóa dữ liệu nhiều chiều bằng cách cung cấp cho mỗi điểm dữ liệu một vị trí trong bản đồ hai hoặc ba chiều sử dụng kỹ thuật Stochastic Neighbor Embedding.

### UMAP

Uniform Manifold Approximation and Projection là một kỹ thuật giảm kích thước chiều dữ liệu được xây dựng từ một khung lý thuyết dựa trên hình học Riemannian và tôpô đại số.

### Random walk

Trong toán học, bước đi ngẫu nhiên là một đối tượng toán học, được gọi là quá trình ngẫu nhiên, mô tả một con đường bao gồm liên tiếp các bước ngẫu nhiên trên một số không gian toán học.

### Nearest-neighbour network

Là một mạng dữ liệu bao gồm các đối tượng liên kết với nhau dựa trên khoảng cách nào đó. Trong scRNA-seq, có thể hiểu là mạng tế bào, trong đó các tế bào có mối quan hệ gần gũi về mặt biểu hiện gene sẽ liên kết với nhau.

### Stochastic

Đề cập đến đặc tính được mô tả bởi phân phối xác suất ngẫu nhiên.

### Seed

Một seed (lấy nghĩa của từ hạt giống) ngẫu nhiên là điểm khởi đầu trong việc tạo ra các số ngẫu nhiên. Khi seed giống nhau, các quá trình ngẫu nhiên sẽ diễn ra giống hệt nhau.

### Random-number generator

Trình tạo ra số ngẫu nhiên, sẽ ảnh hưởng bởi seed ngẫu nhiên.

### Communities

Nghĩa là cộng đồng. Đối với đồ thị, có thể được định nghĩa là một tập hợp con của các nút được kết nối chặt chẽ với nhau và kết nối lỏng lẻo với các nút trong các cộng đồng khác trong cùng một đồ thị.

### Mô-đun

Biểu thị cho gene hoặc nhóm gene.

### Gene marker

Gene được biểu hiện mang tính đại diện cho một nhóm tế bào, một kiểu hình hay một tính chất nào đó của tế bào.
