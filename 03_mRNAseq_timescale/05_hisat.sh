#!/bin/bash
#Map trimmed RNA-seq reads to the L. variegatus 3.0 genome
#K. Castellano

module load histat2/v2.2.1
module load samtools/v1.18

fileList="trim_10FBS_D20_S277 trim_10FBS_D313_S283 trim_10FBS_D445_S286 trim_10FBS_D738_S289 trim_15FBS_D182_S281 trim_15FBS_D20_S278 trim_15FBS_D313_S284 trim_15FBS_D445_S287 trim_15FBS_D738_S290 trim_5FBS_D182_S279 trim_5FBS_D20_S276 trim_5FBS_D313_S282 trim_5FBS_D445_S285 trim_5FBS_D738_S288"
path=<path to files>

for file in ${fileList}
do

exec &> ${file}_hisat2.log
hisat2 -p 8 -x Lvar3.0_hisat2index \
	-1 ${path}${file}_R1.fastq.gz \
	-2 ${path}${file}_R2.fastq.gz | \
	samtools view -@ 8 -S -h -u - | \
	samtools sort -@ 8 -T ${file} - > ${file}.bam

done
