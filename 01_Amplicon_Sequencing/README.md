# Validating the identity of the cells in culture using 18S rDNA targeted amplicon sequencing

## 1. QC raw data using FastQC v0.11.10.devel
- script: 01_fastqc_rawdata.sh

## 2. Trim reads using Trimmomatic version 0.38 and then check the quality of trimmed reads using FastQC v0.11.10.devel
- script: 02_qualityFilter
- important parameters:
    - #SLIDINGWINDOW:4:25 means that if at any point in a scan of the sequence the average quality drops below 25 in a 4bp span, the rest of the read will be cut off
    - #MINLEN:45 indicates that the read should be dropped altogether if it is trimmed shorter than 45bp
    - #CROP: number of bases to keep from the beginning (or number of bases to remove from 3' end)
    - #HEADCROP: number of bases to remove from the beginning of the read

## 3. Merge reads and collapse uniqs using usearch v11.0.667_i86linux64
- program and version: usearch v11.0.667_i86linux64
- input: trimmed forward and reverse reads for each sample from step B above
- script: 03_usearch_merge_collapse

## 4. All final collapsed reads were analyzed using NCBI nucleotide Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
- Results can ve viewed in the excel sheet above
