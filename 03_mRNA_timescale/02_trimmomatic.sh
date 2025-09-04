#!/bin/bash
#K. Castellano

module load Trimmomatic/v0.39

fileList="5FBS_D20_S276 5FBS_D182_S279 5FBS_D313_S282 5FBS_D445_S285 5FBS_D738_S288 10FBS_D20_S277 10FBS_D182_S280 10FBS_D313_S283 10FBS_D445_S286 10FBS_D738_S289 15FBS_D20_S278 15FBS_D182_S281 15FBS_D313_S284 15FBS_D445_S287 15FBS_D738_S290"
path=<path to files>

for file in ${fileList}
do

java -jar /data/resources/app_modules/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 \
	${path}${file}_L004_R1_001.fastq.gz \
	${path}${file}_L004_R2_001.fastq.gz \
        trim_${file}_R1.fastq.gz singles_trim_${file}_R1.fastq.gz \
        trim_${file}_R2.fastq.gz singles_trim_${file}_R2.fastq.gz \
        ILLUMINACLIP:/data/resources/app_modules/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
	HEADCROP:5 \
	CROP:94 \
        SLIDINGWINDOW:4:25 \
        MINLEN:45 2>&1 | tee -a ${file}_trimmomatic_log.txt

done
