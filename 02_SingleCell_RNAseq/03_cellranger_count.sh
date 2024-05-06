#!/bin/bash
# 10x cellranger count - align sequences to transcriptome and generates a .cloupe file for visualization
# Author: Kate Castellano
# ----------------------------

##align all samples to genome (ref made with gene_name attribute added and filtered gtf)

#data for each FBS condition separately
/data/app/cellranger-6.1.2/cellranger count --id=Lv_count_3FBS_expectcells \
--fastqs=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/KCastellano_GMGI_10xsc_10Jan2022/Files/Lv_scRNA_fastq/outs/fastq_path/HWJ2HDRXY \
--sample=Lv_3FBS \
--transcriptome=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/cellranger_count/Lv_genome_filteredGenes_addGeneName_10xindx \
--expect-cells=10000 2>&1 | tee -a cellranger_count_3FBS_log.txt

/data/app/cellranger-6.1.2/cellranger count --id=Lv_count_5FBS \
--fastqs=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/KCastellano_GMGI_10xsc_10Jan2022/Files/Lv_scRNA_fastq/outs/fastq_path/HWJ2HDRXY \
--sample=Lv_5FBS \
--transcriptome=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/cellranger_count/Lv_genome_filteredGenes_addGeneName_10xindx \
--expect-cells=10000 2>&1 | tee -a cellranger_count_5FBS_log.txt

/data/app/cellranger-6.1.2/cellranger count --id=Lv_count_10FBS \
--fastqs=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/KCastellano_GMGI_10xsc_10Jan2022/Files/Lv_scRNA_fastq/outs/fastq_path/HWJ2HDRXY \
--sample=Lv_10FBS \
--transcriptome=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/cellranger_count/Lv_genome_filteredGenes_addGeneName_10xindx \
--expect-cells=10000 2>&1 | tee -a cellranger_count_10FBS_log.txt

/data/app/cellranger-6.1.2/cellranger count --id=Lv_count_15FBS \
--fastqs=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/KCastellano_GMGI_10xsc_10Jan2022/Files/Lv_scRNA_fastq/outs/fastq_path/HWJ2HDRXY \
--sample=Lv_15FBS \
--transcriptome=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/cellranger_count/Lv_genome_filteredGenes_addGeneName_10xindx \
--expect-cells=10000 2>&1 | tee -a cellranger_count_15FBS_log.txt

