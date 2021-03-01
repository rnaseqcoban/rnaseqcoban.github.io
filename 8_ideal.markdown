---
layout: page
title: 8. Pipeline chuẩn cho scRNAseq
permalink: /pipeline/
nav_order: 8
---

#  Ideal scRNAseq pipeline (as of Oct 2017)

## Experimental Design
* Avoid confounding biological and batch effects (Figure 10.1)
	* Multiple conditions should be captured on the same chip if possible
	* Perform multiple replicates of each condition where replicates of different conditions should be performed together if possible
	* Statistics cannot correct a completely confounded experiment!

* Unique molecular identifiers
	* Greatly reduce noise in data
	* May reduce gene detection rates (unclear if it is UMIs or other protocol differences)
	* Lose splicing information
	* Use longer UMIs (~10bp)
	* Correct for sequencing errors in UMIs using UMI-tools
* Spike-ins
	Useful for quality control
	May be useful for normalizing read counts
	Can be used to approximate cell-size/RNA content (if relevant to biological question)
	Often exhibit higher noise than endogenous genes (pipetting errors, mixture quality)
	Requires more sequencing to get enough endogenous reads per cell
* Cell number vs Read depth
	Gene detection plateaus starting from 1 million reads per cell
	Transcription factor detection (regulatory networks) require high read depth and most sensitive protocols (i.e. Fluidigm C1)
	Cell clustering & cell-type identification benefits from large number of cells and doesn’t requireas high sequencing depth (~100,000 reads per cell).


##  Processing Reads
* Read QC & Trimming
	* FASTQC, cutadapt
* Mapping
	* Small datasets or UMI datasets: align to genome/transcriptome using STAR
	* Large datasets: pseudo-alignment with Salmon or kallisto
* Quantification
	* Small dataset, no UMIs : featureCounts
	* Large datasets, no UMIs: Salmon, kallisto
	* UMI dataset : UMI-tools + featureCounts

## Preparing Expression Matrix

* Cell QC
	* scater
	* consider: mtRNA, rRNA, spike-ins (if available), number of detected genes per cell, total reads/molecules per cell
* Library Size Normalization
	* scran
* Batch correction (if appropriate)
	* Replicates/Confounded RUVs
	* Unknown or unbalanced biological groups mnnCorrect
	* Balanced design ComBat
## Biological Interpretation
* Feature Selection
	* M3Drop
* Clustering and Marker Gene Identification
	* ≤ 5000 cells : SC3
	* > 5000 cells : Seurat
* Pseudotime
	* distinct timepoints: TSCAN
	* small dataset/unknown number of branches: Monocle2
	* large continuous dataset: destiny
* Differential Expression
	* Small number of cells and few groups : scde
	* Replicates with batch effects : mixture/linear models
	* Balanced batches: edgeR or MAST
	* Large datasets: Kruskal-Wallis test (all groups at once), or Wilcox-test (compare 2-groups at a time).