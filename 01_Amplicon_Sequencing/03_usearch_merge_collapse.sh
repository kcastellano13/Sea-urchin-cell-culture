
#!/bin/bash
#merge R1 and R2 reads and then collapse reads so we have on read per region
#usearch version v11.0.667_i86linux64
#K Castellano

#merge forward and reverse reads
usearch -fastq_mergepairs Lv-D22-3FBS_S67_L001_R1_001_trim.fastq \
-reverse Lv-D22-3FBS_S67_L001_R2_001_trim.fastq \
-fastaout Lv-D22-3FBS_S67_merged.fasta

usearch -fastq_mergepairs Lv-D22-5FBS_S64_L001_R1_001_trim.fastq \
-reverse Lv-D22-5FBS_S64_L001_R2_001_trim.fastq \
-fastaout Lv-D22-5FBS_S64_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D22-10FBS_S65_L001_R1_001_trim.fastq \
-reverse Lv-D22-10FBS_S65_L001_R2_001_trim.fastq \
-fastaout Lv-D22-10FBS_S65_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D22-15FBS_S66_L001_R1_001_trim.fastq \
-reverse Lv-D22-15FBS_S66_L001_R2_001_trim.fastq \
-fastaout Lv-D22-15FBS_S66_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D112-3FBS_S67_L001_R1_001_trim.fastq \
-reverse Lv-D112-3FBS_S67_L001_R2_001_trim.fastq \
-fastaout Lv-D112-3FBS_S67_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D112-5FBS_S68_L001_R1_001_trim.fastq \
-reverse Lv-D112-5FBS_S68_L001_R2_001_trim.fastq \
-fastaout Lv-D112-5FBS_S68_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D112-10FBS_S69_L001_R1_001_trim.fastq \
-reverse Lv-D112-10FBS_S69_L001_R2_001_trim.fastq \
-fastaout Lv-D112-10FBS_S69_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D112-15FBS_S70_L001_R1_001_trim.fastq \
-reverse Lv-D112-15FBS_S70_L001_R2_001_trim.fastq \
-fastaout Lv-D112-15FBS_S70_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D70-P3-EB_S71_L001_R1_001_trim.fastq \
-reverse Lv-D70-P3-EB_S71_L001_R2_001_trim.fastq \
-fastaout Lv-D70-P3-EB_S71_merged.fasta 2>&1 | tee -a usearch.log

usearch -fastq_mergepairs Lv-D229-P3-MrBig_S72_L001_R1_001_trim.fastq \
-reverse Lv-D229-P3-MrBig_S72_L001_R2_001_trim.fastq \
-fastaout Lv-D229-P3-MrBig_S72_merged.fasta 2>&1 | tee -a usearch.log


#collapse reads
fileList="*merged.fasta"

for file in ${fileList}
do

    #collapse unique reads
    #size out will add the number of reads that fall into the one unique sequence
    #must be atleast 90% identical to cluster (-id 0.9)
    prefix=$(echo ${file} | cut -d "." -f 1)
    usearch -cluster_fast ${file} -id 0.9 -sizeout -centroids ${prefix}_uniq.fasta 2>&1 | tee -a usearch.log

done