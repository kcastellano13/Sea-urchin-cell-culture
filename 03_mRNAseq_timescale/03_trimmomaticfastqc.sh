#!/bin/bash
#fastqc of trimmed RNA-seq data
#K. Castellano

module load fastQC/v0.12.1

fileList="*.gz"

for file in ${fileList}
do

	fastqc ${file}

done
