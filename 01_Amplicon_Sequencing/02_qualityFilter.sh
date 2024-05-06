
#!/bin/bash
#Quality filter 18S reads using the Trimmomatic version 0.38 and fastqc of trimmed reads
#Trimmomatic version 0.38 and FastQC v0.11.10.devel
#KCastellano

##Notes on flags used for trimming
#SLIDINGWINDOW:4:25 means that if at any point in a scan of the sequence the average quality drops below 25 in a 4bp span, the rest of the read will be cut off
#MINLEN:45 indicates that the read should be dropped altogether if it is trimmed shorter than 45bp
#CROP: number of bases to keep from the beginning (or number of bases to remove from 3' end)
#HEADCROP: number of bases to remove from the beginning of the read

path=<location of raw fastq reads>

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D22-3FBS_S63_L001_R1_001.fastq \
        ${path}Lv-D22-3FBS_S63_L001_R2_001.fastq \
        Lv-D22-3FBS_S63_L001_R1_001_trim.fastq Lv-D22-3FBS_S63_L001_R1_001_single.fastq \
        Lv-D22-3FBS_S63_L001_R2_001_trim.fastq Lv-D22-3FBS_S63_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D22-5FBS_S64_L001_R1_001.fastq \
        ${path}Lv-D22-5FBS_S64_L001_R2_001.fastq \
        Lv-D22-5FBS_S64_L001_R1_001_trim.fastq Lv-D22-5FBS_S64_L001_R1_001_single.fastq \
        Lv-D22-5FBS_S64_L001_R2_001_trim.fastq Lv-D22-5FBS_S64_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D22-10FBS_S65_L001_R1_001.fastq \
        ${path}Lv-D22-10FBS_S65_L001_R2_001.fastq \
        Lv-D22-10FBS_S65_L001_R1_001_trim.fastq Lv-D22-10FBS_S65_L001_R1_001_single.fastq \
        Lv-D22-10FBS_S65_L001_R2_001_trim.fastq Lv-D22-10FBS_S65_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D22-15FBS_S66_L001_R1_001.fastq \
        ${path}Lv-D22-15FBS_S66_L001_R2_001.fastq \
        Lv-D22-15FBS_S66_L001_R1_001_trim.fastq Lv-D22-15FBS_S66_L001_R1_001_single.fastq \
        Lv-D22-15FBS_S66_L001_R2_001_trim.fastq Lv-D22-15FBS_S66_L001_R1_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D112-3FBS_S67_L001_R1_001.fastq \
        ${path}Lv-D112-3FBS_S67_L001_R2_001.fastq \
        Lv-D22-3FBS_S67_L001_R1_001_trim.fastq Lv-D22-3FBS_S67_L001_R1_001_single.fastq \
        Lv-D22-3FBS_S67_L001_R2_001_trim.fastq Lv-D22-3FBS_S67_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D112-5FBS_S68_L001_R1_001.fastq \
        ${path}Lv-D112-5FBS_S68_L001_R2_001.fastq \
        Lv-D112-5FBS_S68_L001_R1_001_trim.fastq Lv-D22-5FBS_S68_L001_R1_001_single.fastq \
        Lv-D112-5FBS_S68_L001_R2_001_trim.fastq Lv-D22-5FBS_S68_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D112-10FBS_S69_L001_R1_001.fastq \
        ${path}Lv-D112-10FBS_S69_L001_R2_001.fastq \
        Lv-D112-10FBS_S69_L001_R1_001_trim.fastq Lv-D22-10FBS_S69_L001_R1_001_single.fastq \
        Lv-D112-10FBS_S69_L001_R2_001_trim.fastq Lv-D22-10FBS_S69_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D112-15FBS_S70_L001_R1_001.fastq \
        ${path}Lv-D112-15FBS_S70_L001_R2_001.fastq \
        Lv-D112-15FBS_S70_L001_R1_001_trim.fastq Lv-D22-15FBS_S70_L001_R1_001_single.fastq \
        Lv-D112-15FBS_S70_L001_R2_001_trim.fastq Lv-D22-15FBS_S70_L001_R1_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D112-15FBS_S70_L001_R1_001.fastq \
        ${path}Lv-D112-15FBS_S70_L001_R2_001.fastq \
        Lv-D112-15FBS_S70_L001_R1_001_trim.fastq Lv-D22-15FBS_S70_L001_R1_001_single.fastq \
        Lv-D112-15FBS_S70_L001_R2_001_trim.fastq Lv-D22-15FBS_S70_L001_R1_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D70-P3-EB_S71_L001_R1_001.fastq \
        ${path}Lv-D70-P3-EB_S71_L001_R2_001.fastq \
        Lv-D70-P3-EB_S71_L001_R1_001_trim.fastq Lv-D70-P3-EB_S71_L001_R1_001_single.fastq \
        Lv-D70-P3-EB_S71_L001_R2_001_trim.fastq Lv-D70-P3-EB_S71_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 \
        ${path}Lv-D229-P3-MrBig_S72_L001_R1_001.fastq \
        ${path}Lv-D229-P3-MrBig_S72_L001_R2_001.fastq \
        Lv-D229-P3-MrBig_S72_L001_R1_001_trim.fastq Lv-D229-P3-MrBig_S72_L001_R1_001_single.fastq \
        Lv-D229-P3-MrBig_S72_L001_R2_001_trim.fastq Lv-D229-P3-MrBig_S72_L001_R2_001_single.fastq \
        SLIDINGWINDOW:4:25 \
        CROP:230 \
        HEADCROP:9 \
        MINLEN:45 2>&1 | tee -a trimmomatic.log

##run Fastqc of alltrimmed files

fileList="*.trim.fastq"

for file in ${fileList}
do

        fastqc ${file} -t 5 2>&1 | tee -a trimmomatic.log

done