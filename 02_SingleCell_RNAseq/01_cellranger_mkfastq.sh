#!/bin/bash
#10x single cell sequencing - convert base calls to fastq files
#Author: Kate Castellano
# ----------------------------

/data/app/cellranger-6.1.2/cellranger mkfastq --id=Lv_scRNA --run=/data/prj/urchin/cell-culture/2022Jan_10xscSeq/KCastellano_GMGI_10xsc_10Jan2022/Files \
--samplesheet=cellranger_sample.csv 2>&1 | tee -a cellranger_mkfastq_log.txt
