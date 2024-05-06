#!/bin/bash
#fastqc of raw RNA-seq data

module load fastQC/v0.12.1

fileList="*.gz"

for file in ${fileList}
do

	fastqc ${file}

done
