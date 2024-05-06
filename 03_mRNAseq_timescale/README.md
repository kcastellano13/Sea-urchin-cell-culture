# Differential Expression analysis of L. variegatus Cell Culture time series
### Samples: 5%, 10% and 15% each on Day 20, 182, 313, 445, 738

## 1. Check Data quality 
- program: fastqc/v0.12.1
- script: 01_fastqc.sh

## 2. Trim and quality filter reads
- program: Trimmomatic/v0.39
- script: 02_trimmomatic.sh

## 3. Check Data quality after trimming (fastqc/v0.12.1 nad multiQC/v1.16.
- program: fastqc/v0.12.1
- script: 03_trimmomaticfastqc.sh

## 4. Map reads to the L. variegatus genome
- program: histat2/v2.2.1
    ##### A. Make the L. var v3.0 genome into a hisat index
    - script: 04_hisat_index.sh

    ##### B. Map read to L. var v3.0 genome 
    - script: 05_hisat.sh

    ##### C. Alignment QC - get summary stats to make sure the alignment itself went ok 
    - script: 06_hisat_summaryStats.sh

## 5. Create count files for TPM analysis
- program: htseq/v2.0.4 - referenced based
- script: 07_htseq_count.sh

## 6. Get Gene lengths for TPM analysis
- htseq only counts alignments to exons so use GTF tools to get exon lengths for TPM analysis
- program: gtf tools v0.9.0 (https://www.genemine.org/gtftools.php.

    ##### A. Rename L. var chromsomes because gtf tools only accepts chromosome in human format (so 1-22.
    Edit file with sed - Chromosome number associations are from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018143015.1/.
    ```
    sed -i 's/CM031015.1/1/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031016.1/2/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031017.1/3/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031018.1/4/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031019.1/5/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031020.1/6/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031021.1/7/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031022.1/8/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031023.1/9/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031024.1/10/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031025.1/11/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031026.1/12/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031027.1/13/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031028.1/14/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031029.1/15/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031030.1/16/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031031.1/17/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031032.1/18/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/CM031033.1/19/g' GCF_018143015.1_Lvar_3.0_genomic.gtf
    sed -i 's/MG676468.1/20/g' GCF_018143015.1_Lvar_3.0_genomic.gtf #mitochondrial but just labeled as a number for gtf tools
    ```
    ##### B. Get gene lengths of merged exons via gtf tools
    ###### -l = Calculate gene lengths. Since a gene may have multiple isoforms, there are multiple ways to calculate gene length based on literature. Three simple ways are considering the mean, median and maximum of the lengths of isoforms as the length of the gene. A fourth way is to calculate the length of merged exons of all isoforms (i.e. non-overlapping exonic length.. So, in total, four different types of gene lengths(the mean, median and max of lengths of isoforms of agene, and the length of merged exons of isoforms of a gene. are provided. format. Needed for e.g. calculating FPKM in RNA-seq data analysis, where gene length is required.

```
python GTFtools_0.9.0/gtftools.py -l Lvar3.0_exon_length.txt GCF_018143015.1_Lvar_3.0_genomic.gtf
#make a file with the gene ID and the merged column to use for TPM analysis. I do not care about isoforms in this case and want all counts associated with the gene so that is why I chose the merged exon count column
awk -v OFS="\t" '{print $1,$5}' Lvar3.0_exon_length.txt > Lvar3.0_exonMerged_length.txt
```

## 7. Prepare data for TPM analysis in R v2023.09.1
```{r import and reformat tables}
#in excel manually remove rows at the bottom (alignment not unique, ambiguous, no feature, not aligned, too low quality.
FBS5_D20 <- data.frame(read.table("trim_5FBS_D20_S276.counts", sep="\t"..
colnames(FBS5_D20. <- c("target_id", "FBS5_D20_count".
FBS5_D182 <- data.frame(read.table("trim_5FBS_D182_S279.counts", sep="\t"..
colnames(FBS5_D182. <- c("target_id", "FBS5_D182_count".
FBS5_D313 <- data.frame(read.table("trim_5FBS_D313_S282.counts", sep="\t"..
colnames(FBS5_D313. <- c("target_id", "FBS5_D313_count".
FBS5_D445 <- data.frame(read.table("trim_5FBS_D445_S285.counts", sep="\t"..
colnames(FBS5_D445. <- c("target_id", "FBS5_D445_count".
FBS5_D738 <- data.frame(read.table("trim_5FBS_D738_S288.counts", sep="\t"..
colnames(FBS5_D738. <- c("target_id", "FBS5_D738count".

FBS10_D20 <- data.frame(read.table("trim_10FBS_D20_S277.counts", sep="\t"..
colnames(FBS10_D20. <- c("target_id", "FBS10_D20_count".
FBS10_D182 <- data.frame(read.table("trim_10FBS_D182_S280.counts", sep="\t"..
colnames(FBS10_D182. <- c("target_id", "FBS10_D182_count".
FBS10_D313 <- data.frame(read.table("trim_10FBS_D313_S283.counts", sep="\t"..
colnames(FBS10_D313. <- c("target_id", "FBS10_D313_count".
FBS10_D445 <- data.frame(read.table("trim_10FBS_D445_S286.counts", sep="\t"..
colnames(FBS10_D445. <- c("target_id", "FBS10_D445_count".
FBS10_D738 <- data.frame(read.table("trim_10FBS_D738_S289.counts", sep="\t"..
colnames(FBS10_D738. <- c("target_id", "FBS10_D738count".

FBS15_D20 <- data.frame(read.table("trim_15FBS_D20_S278.counts", sep="\t"..
colnames(FBS15_D20. <- c("target_id", "FBS15_D20_count".
FBS15_D182 <- data.frame(read.table("trim_15FBS_D182_S281.counts", sep="\t"..
colnames(FBS15_D182. <- c("target_id", "FBS15_D182_count".
FBS15_D313 <- data.frame(read.table("trim_15FBS_D313_S284.counts", sep="\t"..
colnames(FBS15_D313. <- c("target_id", "FBS15_D313_count".
FBS15_D445 <- data.frame(read.table("trim_15FBS_D445_S287.counts", sep="\t"..
colnames(FBS15_D445. <- c("target_id", "FBS15_D445_count".
FBS15_D738 <- data.frame(read.table("trim_15FBS_D738_S290.counts", sep="\t"..
colnames(FBS15_D738. <- c("target_id", "FBS15_D738count".
```

```{r merge tables}
#merge table to make count table
### merge subsets -- all stages ###
countsTable <- Reduce(function (x,y. merge(x=x, y=y, by="target_id"., list(FBS5_D20, FBS5_D182, FBS5_D313, FBS5_D445, FBS5_D738, FBS10_D20,FBS10_D182, FBS10_D313, FBS10_D445, FBS10_D738, FBS15_D20, FBS15_D182, FBS15_D313, FBS15_D445, FBS15_D738..


#remove rows 1-5 which have some summary stats (alignment not unique, ambiguous, no feature, not aligned, too low quality.
countsTable <- countsTable[-c(1, 2, 3, 4, 5., ]
#reset row numbers to start with 1
rownames(countsTable. <- 1:nrow(countsTable.

#write.table(countsTable, "FBStimecourse_combined_readCount.tsv", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE.
```

```{r }
## remove any genes that are not covered by any of the samples
## ie, if there are no read counts in any of the samples, remove this gene
## ASC can handle zero count data, but I filter it out anyway
fdta <- countData_m[!apply(countData_m[, 2:4], 1, function(x. all(x==0..,]
length(fdta[,1]. #22131
#original file
length(countData_m[,1]. #26301 - 22131 = 4,170 genes removed with no counts in any   condition
```
#at this point the counts table can be used to make MDS plot (I transfer merged counts table to computer and run in R - uses TCseq to make MDS plot.


```{r }
#transpose rows and columns for TPM analysis to make each column represent a gene and each row represent a sample
countsTable_transform <- t(countsTable. #this creates a transposed matrix
#convert back to data frame
test <- as.data.frame(countsTable_transform.
write.table(test, "FBStimecourseAll_noZeros_Transpose.tsv", sep = "\t", quote = FALSE.

# and without row names for TPM pipeline
write.table(test, "FBStimecourseAll_noZeros_Transpose.tsv", row.names=FALSE, sep = "\t", quote = FALSE.

```
## 8. Convert counts to TPM
- _program:_ python3
- _script:_ 08_convert_to_tpm.py

## 9. transpose the file again to swap rows and columns so that the rows represent genes and columns represent samples, done in R
```{r }
#read table back in
countsTable <- read.table("FBStimecourseAll_noZeros_TPM.csv", sep = "," , header = TRUE.

countsTable_transform <- as.data.frame(t(countsTable.. #the transpose function will convert your data to a matrix so to keep it as a data frame, use as.data.frame

write.table(countsTable_transform, "FBStimecourseAll_noZeros_TPM_Transpose.tsv", sep = "\t", quote = FALSE.
```
