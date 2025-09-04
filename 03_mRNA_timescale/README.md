# Differential Expression analysis of *L. variegatus* Cell Culture time series
## Pre-processing of Data done by Kate Castellano
### Samples: 5%, 10% and 15% each on Day 20, 182, 313, 445, 738 (done in singlicate)

## 1. Check Data quality 
- program: fastqc/v0.12.1
- script: 01_fastqc.sh

## 2. Trim and quality filter reads
- program: Trimmomatic/v0.39
- script: 02_trimmomatic.sh

## 3. Check Data quality after trimming
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
- program: gtf tools v0.9.0 (https://www.genemine.org/gtftools.php)

##### A. Rename L. var chromsomes because gtf tools only accepts chromosome in human format (so 1-22)

Edit file with sed - Chromosome number associations are from NCBI ([GCF_018143015.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018143015.1/))
    
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
#run GTF tools
python GTFtools_0.9.0/gtftools.py -l Lvar3.0_exon_length.txt GCF_018143015.1_Lvar_3.0_genomic.gtf

#make a file with the gene ID and the merged exon column to use for TPM analysis.
#I do not care about isoforms in this case and want all counts associated with the gene so that is why I chose the merged exon count column
awk -v OFS="\t" '{print $1,$5}' Lvar3.0_exon_length.txt > Lvar3.0_exonMerged_length.txt
```

## 7. Prepare data for TPM analysis
- Program: R v2023.09.1
```{r import and reformat tables}
#in excel manually remove rows at the bottom (alignment not unique, ambiguous, no feature, not aligned, too low quality)
FBS5_D20 <- data.frame(read.table("trim_5FBS_D20_S276.counts", sep="\t"))
colnames(FBS5_D20. <- c("target_id", "FBS5_D20_count")
FBS5_D182 <- data.frame(read.table("trim_5FBS_D182_S279.counts", sep="\t"))
colnames(FBS5_D182. <- c("target_id", "FBS5_D182_count")
FBS5_D313 <- data.frame(read.table("trim_5FBS_D313_S282.counts", sep="\t"))
colnames(FBS5_D313. <- c("target_id", "FBS5_D313_count")
FBS5_D445 <- data.frame(read.table("trim_5FBS_D445_S285.counts", sep="\t"))
colnames(FBS5_D445. <- c("target_id", "FBS5_D445_count")
FBS5_D738 <- data.frame(read.table("trim_5FBS_D738_S288.counts", sep="\t"))
colnames(FBS5_D738. <- c("target_id", "FBS5_D738count")

FBS10_D20 <- data.frame(read.table("trim_10FBS_D20_S277.counts", sep="\t"))
colnames(FBS10_D20. <- c("target_id", "FBS10_D20_count")
FBS10_D182 <- data.frame(read.table("trim_10FBS_D182_S280.counts", sep="\t"))
colnames(FBS10_D182. <- c("target_id", "FBS10_D182_count")
FBS10_D313 <- data.frame(read.table("trim_10FBS_D313_S283.counts", sep="\t"))
colnames(FBS10_D313. <- c("target_id", "FBS10_D313_count")
FBS10_D445 <- data.frame(read.table("trim_10FBS_D445_S286.counts", sep="\t"))
colnames(FBS10_D445. <- c("target_id", "FBS10_D445_count")
FBS10_D738 <- data.frame(read.table("trim_10FBS_D738_S289.counts", sep="\t"))
colnames(FBS10_D738. <- c("target_id", "FBS10_D738count")

FBS15_D20 <- data.frame(read.table("trim_15FBS_D20_S278.counts", sep="\t"))
colnames(FBS15_D20. <- c("target_id", "FBS15_D20_count")
FBS15_D182 <- data.frame(read.table("trim_15FBS_D182_S281.counts", sep="\t"))
colnames(FBS15_D182. <- c("target_id", "FBS15_D182_count")
FBS15_D313 <- data.frame(read.table("trim_15FBS_D313_S284.counts", sep="\t"))
colnames(FBS15_D313. <- c("target_id", "FBS15_D313_count")
FBS15_D445 <- data.frame(read.table("trim_15FBS_D445_S287.counts", sep="\t"))
colnames(FBS15_D445. <- c("target_id", "FBS15_D445_count")
FBS15_D738 <- data.frame(read.table("trim_15FBS_D738_S290.counts", sep="\t"))
colnames(FBS15_D738. <- c("target_id", "FBS15_D738count")
```

```{r merge tables}
#merge table from each sample to make one big count table
countsTable <- Reduce(function (x,y. merge(x=x, y=y, by="target_id"., list(FBS5_D20, FBS5_D182, FBS5_D313, FBS5_D445, FBS5_D738, FBS10_D20,FBS10_D182, FBS10_D313, FBS10_D445, FBS10_D738, FBS15_D20, FBS15_D182, FBS15_D313, FBS15_D445, FBS15_D738..

#remove rows 1-5 which have some summary stats (alignment not unique, ambiguous, no feature, not aligned, too low quality)
countsTable <- countsTable[-c(1, 2, 3, 4, 5., )]

#reset row numbers to start with 1
rownames(countsTable. <- 1:nrow(countsTable)

#write.table(countsTable, "FBStimecourse_combined_readCount.tsv", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
```

```{r }
## remove any genes that are not covered by any of the samples
## ie, if there are no read counts in any of the samples, remove this gene
## ASC can handle zero count data, but I filter it out anyway
fdta <- countData_m[!apply(countData_m[, 2:4], 1, function(x. all(x==0..,]
length(fdta[,1]) #22131

#original file
length(countData_m[,1]) #26301

#26301 - 22131 = 4,170 genes removed with no counts in any condition
```

- at this point the counts table can be used to make MDS plot (I transfer merged counts table to computer and run in R - uses TCseq to make MDS plot.

```{r }
#transpose rows and columns for TPM analysis to make each column represent a gene and each row represent a sample
countsTable_transform <- t(countsTable. #this creates a transposed matrix

#convert back to data frame
test <- as.data.frame(countsTable_transform.
write.table(test, "FBStimecourseAll_noZeros_Transpose.tsv", sep = "\t", quote = FALSE.

# also write file without row names for TPM pipeline
write.table(test, "FBStimecourseAll_noZeros_Transpose.tsv", row.names=FALSE, sep = "\t", quote = FALSE.
```

## 8. Convert counts to TPM
- program: python3
- script: 08_convert_to_tpm.py

## 9. Transpose the file again to swap rows and columns so that the rows represent genes and columns represent samples
- program: R v2023.09.1
```{r }
#read table back in
countsTable <- read.table("FBStimecourseAll_noZeros_TPM.csv", sep = "," , header = TRUE.

countsTable_transform <- as.data.frame(t(countsTable.. #the transpose function will convert your data to a matrix so to keep it as a data frame, use as.data.frame

write.table(countsTable_transform, "Lv_FBStimecourseAll_noZeros_TPM.tsv", sep = "\t", quote = FALSE.
```

## Clustering of time-course RNAseq data.
### Analysis by RMK is based on methods in Chille et al 2021: https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/3-WGCNA/WGCNA.Rmd
### Manuscript DOI: https://doi.org/10.1242/jeb.243187

```{r }
#Import necessary libraries
library("tidyverse")
library("genefilter")
# Check the version of the genefilter package
packageVersion("genefilter")
library("DESeq2")
packageVersion("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
packageVersion("ComplexHeatmap")
library("goseq")
library("dplyr")
library("clusterProfiler")
library("simplifyEnrichment")
library("Rmisc")
library("factoextra")
library("extrafont")
library("ggsignif")
library("rmarkdown")
#install.packages("magick")
library("magick")
packageVersion("NbClust")

##Data input, cleaning, and pre-processing

#Read in the TPM count data from the CSV file TPM_COUNTS_LvarCulturedCells.csv
#The sheet was originally titled "Lv_FBStimecourseAll_noZeros_TPM", and I renamed it to TPM_COUNTS_LvarCulturedCells.

#Read in the original data with Annotations column, save for later
original_gene_matrix_annotated <- as.data.frame(read.csv("TPM_COUNTS_LvarCulturedCells.csv", header=TRUE, row.names="FeatureName"))
original_gene_matrix_annotated$gene_id <- rownames(original_gene_matrix_annotated) #add the column "gene_id"

#Read in only gene counts (no annotation column) as "gcount" and proceed with filtering.
gcount <- as.data.frame(read.csv("TPM_COUNTS_LvarCulturedCells.csv", header=TRUE, row.names="FeatureName"))
head(gcount)
#remove "annotation" column
library(dplyr)
gcount <- gcount %>% select(-Annotation)
dim(gcount) #9984 genes x 15 samples
summary(gcount) #quick way to see if any have NAN values


#Read in the treatment info file: a metadata .csv
metadata <- read.csv("metadata_lvar_cellcultureexpt.csv", header = TRUE, sep = ",")
#This is a .csv file with 3 columns: sample_name, condition, and time_point.

# Create the new column 'condition_time_point', which is a merge of condition and time_point
metadata <- metadata %>%
  mutate(condition_time_point = paste(condition, time_point, sep = "_"))
#Treatment is 5, 10, or 15% FBS
head(metadata)

##Quality-filter gene counts
### Check that there are no genes with 0 counts across all samples: if there are, remove those genes

nrow(gcount) #9,984
gcount<-gcount %>%
  mutate(Total = rowSums(.[, 1:15]))%>% #adds a new column Total, which is the sum of values in the first 15 columns for each row
  filter(!Total==0)%>% #filters out rows where total=0
  dplyr::select(!Total) #removes the total column from the dataframe

nrow(gcount) #9,984
#Removed nothing, since zeros had been removed previously.

### Conduct data filtering using pOverA

# The following info taken from Ariana Huffmyer and Erin Chille scripts: 
#   *pOverA*: Specifying the minimum count for a proportion of samples for each gene. Here, we are using a pOverA of 0.07. This is because we have 15 samples with a minimum of n=1 samples per condition per time. Therefore, we will accept genes that are present in 1/15 = 0.07 of the samples because we expect different expression by group. 
# 
# We are further setting the minimum count of genes to 10, such that 10% of the samples must have a gene count of >10 in order for the gene to remain in the data set.  
# 
# Filter in the package "genefilter". Pre-filtering our dataset to reduce the memory size dataframe, increase the speed of the transformation and testing functions, and improve quality of statistical analysis by removing low-coverage counts. Removed counts could represent outliers in the data and removing these improves sensitivity of statistical tests.


filt <- filterfun(pOverA(0.07,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before = 9,984
nrow(gcount_filt) #After = 3,527

# Before filtering, we had 9,984 genes. After filtering for pOverA, we have approximately 3,527 genes. This indicates that there were many genes present in <7% of samples at <10 counts per gene.  
# 
# In order for the DESeq2 algorithms to work, the SampleIDs on the metadata file and count matrices have to match exactly and in the same order. The following R clump will check to make sure that these match. Should return TRUE. 
```

# Display current order of metadata and gene count matrix.  
metadata$sample_name
colnames(gcount_filt)

#Order metadata the same as the column order in the gene matrix.  
list<-colnames(gcount_filt)
list<-as.factor(list)
list

metadata$sample_name<-as.factor(metadata$sample_name)

# Re-order the levels
metadata$sample_name <- factor(as.character(metadata$sample_name), levels=list)
metadata$sample_name

# Re-order the data.frame
metadata_ordered <- metadata[order(metadata$sample_name),]
metadata_ordered$sample_name

#Checking that all row and column names match. Should return "TRUE"
all(rownames(metadata$sample_name) %in% colnames(gcount_filt))
all(rownames(metadata$sample_name) == colnames(gcount_filt)) 


### Plot the sequence counts after filtering sample to confirm that sequence depth is not different between groups
#Visualize distribution of total read counts across different groups to check for sig differences in sequence depth
# Calculate total read counts per sample
total_counts <- colSums(gcount_filt)

# Create a data frame for visualization
sample_ids <- colnames(gcount_filt)
group_info <- metadata$condition_time_point  # Adjust this to match your metadata structure
visualization_data <- data.frame(sample_id = sample_ids, total_counts = total_counts, group = group_info)

# Create boxplot
ggplot(visualization_data, aes(x = group, y = total_counts)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Total Read Counts per Sample by Group", x = "Group", y = "Total Read Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#5FBS_Day_182 and 5FBS_Day_313 have lower total read counts compared to other samples


##################Set DESeq2 design###############################
#NOTE this is not DOING DE; it is simply creating a DESeqDataSet object from the count matrix.
#This object holds the count matrix, the metadata, and the design.
#DESeq2 requires INTEGERs. If you have TPM as input, the values need to be rounded.

gcount_filt_rounded <- round(gcount_filt)

metadata_ordered <- metadata_ordered %>% unite(group, condition, time_point, remove = FALSE) #combines the condition and time_point columns into a new column called group. remove=False keeps the original time_point and treatment columns in the dataframe.

#Set DESeq2 design using group
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt_rounded,
                                colData = metadata_ordered,
                                design = ~group)

# First log-transform the data using a variance stabilizing transformation (VST). This is only for visualization purposes. Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.

# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4.

# Chunk should return TRUE if <4.  

### Estimate size factors for gdds
SF.gdds <- estimateSizeFactors(gdds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst

print(sizeFactors(SF.gdds)) #View size factors
all(sizeFactors(SF.gdds)) < 4 #returned TRUE

#All size factors are less than 4, so we can use VST transformation.  

## Apply VST transformation
gvst <- vst(gdds, blind=TRUE) #apply a variance stabilizing transformation to minimize effects of small counts and normalize wrt library size. Setting to blind=TRUE ensures the transformation is performed without using info about the experimental design, which is done when you dont have replicates.

# NOTE If you set blind=FALSE you get this error due to lack of replicates:
#   Error in estimateDispersionsGeneEst(object.sub, quiet = TRUE) : 
#   the number of samples and the number of model coefficients are equal,
# i.e., there are no replicates to estimate the dispersion.
# use an alternate design formula

head(assay(gvst), 3) #view transformed gene count data for the first three genes in the dataset.

#Convert the VST-transformed counts matrix stored in 'gvst' to a dataframe
gcount_filt_vst <- as.data.frame(assay(gvst)) ### Converts the VST-transformed counts matrix (gvst) to a dataframe.
nrow(gcount_filt) #3,527
nrow(gcount_filt_vst) #same total before and after, 3527 genes.



#### Differential Gene Expression Analysis
#### Run DE analysis using a Wald model. 

#Group model
#DEG <- DESeq(gdds) #run differential expression test by group using the Wald model

#ERROR
#estimating size factors
#estimating dispersions
#Error in checkForExperimentalReplicates(object, modelMatrix) : 
  
#  The design matrix has the same number of samples and coefficients to fit,
#so estimation of dispersion is not possible. Treating samples
#as replicates was deprecated in v1.20 and no longer supported since v1.22.


#It is not possible to de DESeq without replicates.
#But we still find the optimal number of k-means clusters for the VST-transformed data

### Compute optimal number for clustering
library(NbClust)
#NOTE: this failed when trying to use TPM counts (gcount_filt). It succeeded when using VST transformed TPM counts (gcount_filt_vst). 
#It failed when using scaled TPM counts.
##The error thrown when trying to use just TPM normalized values, gcount_filt:
#Error in if ((resCritical[ncB - min_nc + 1, 3] >= alphaBeale) && (!foundBeale)) { : 
#    missing value where TRUE/FALSE needed
#  In addition: Warning message:
#    In pf(beale, pp, df2) : NaNs produced

#Run NbClust on gcount_filt_vst
nb <- NbClust(gcount_filt_vst, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans") 
#use k-means clustering, considering a minimum of 2 clusters and a max of 10
# ** : The Hubert index is a graphical method of determining the number of clusters.
# In the plot of Hubert index, we seek a significant knee that corresponds to a 
# significant increase of the value of the measure i.e the significant peak in Hubert
# index second differences plot. 
# 
# *** : The D index is a graphical method of determining the number of clusters. 
# In the plot of D index, we seek a significant knee (the significant peak in Dindex
#                                                     second differences plot) that corresponds to a significant increase of the value of
# the measure. 
# 
# ******************************************************************* 
#   * Among all indices:                                                
#   * 12 proposed 2 as the best number of clusters 
# * 6 proposed 3 as the best number of clusters 
# * 1 proposed 4 as the best number of clusters 
# * 1 proposed 5 as the best number of clusters 
# * 1 proposed 7 as the best number of clusters 
# * 1 proposed 8 as the best number of clusters 
# * 1 proposed 9 as the best number of clusters 
# * 1 proposed 10 as the best number of clusters 
# 
# ***** Conclusion *****                            
#   
#   * According to the majority rule, the best number of clusters is  2 

#NOTE: nbclust suggests 2

##To visualize the vst-transformed counts in a heatmap, they need to be centered.
#Subtracts the mean expression value of each row from each element in that row (centers the data for each gene around zero)
gcount_filt_vst_centered <- gcount_filt_vst - rowMeans(gcount_filt_vst) 

#Do k-means clustering on the centered VST matrix, specifying 3 clusters.
calc_kmeans <- kmeans(gcount_filt_vst_centered, 3)

#Create a dataframe of the clustering results
cluster_res <- data.frame(gene_id = names(calc_kmeans$cluster), cluster = calc_kmeans$cluster)

# Total number of genes
total_genes <- nrow(cluster_res)
cat("Total number of genes:", total_genes, "\n") #3,527 total genes

# Number of genes in each cluster
cluster_counts <- table(cluster_res$cluster)
print(cluster_counts)

#cluster 1: 1,466 genes
#cluster 2: 922 genes
#cluster 3: 1,139 genes


#Merge the original_gene_matrix_annotated to add the clustering assignments to each gene
library(dplyr)

### Merge cluster_res with detailed gene descriptions using the original_gene_matrix_annotated file of TPM counts
genes_list_with_clusters_and_genedescriptions <- merge(cluster_res, original_gene_matrix_annotated, by = "gene_id", all.x = TRUE)

#Save the result
write.csv(genes_list_with_clusters_and_genedescriptions, file="gene_list_with_clusters.csv")


# First, make sure the cluster assignments are in the same order as the rows in mat_DEG2_centered
gene_clusters <- cluster_res[match(rownames(gcount_filt_vst_centered), cluster_res$gene_id), "cluster"]


# Convert to a factor with desired cluster order and labels
gene_clusters <- factor(gene_clusters, levels = c(1, 2, 3), labels = c("Cluster1", "Cluster2", "Cluster3"))


#Save gene_list_with_clusters_genedescriptions
write.csv(genes_list_with_clusters_and_genedescriptions, file="geneslist_with_clusters_and_genedescriptions.csv")


########### Plot a heatmap of centered, VST transformed TPM counts#######
#Create a new data frame TPM_clust containing only the gene_id and cluster columns.
TPM_clust <- subset(genes_list_with_clusters_and_genedescriptions, select = c(gene_id, cluster))
nrow(TPM_clust) #3257

#Prepare heatmap annotations ('hm_ann') for rows (individual genes per cluster) and columns (condition and time)
hm_ann_row <- unique(TPM_clust) #remove duplicate rows, if any, and assign to hm_ann_row
nrow(hm_ann_row) #3257

rownames(hm_ann_row) <- hm_ann_row$gene_id #set row names of hm_ann_row to the gene_id values
hm_ann_row <- subset(hm_ann_row, select=cluster) #keep only the cluster column

#replaces numeric cluster labels with more descriptive names (Cluster1, Cluster2).
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)  
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)
hm_ann_row$cluster <- gsub(3,"Cluster3",hm_ann_row$cluster)

#converts hm_ann_row to a matrix, ensuring the rows match those in gcount_filt_vst_centered
hm_ann_row <- as.matrix(hm_ann_row[rownames(gcount_filt_vst_centered),])
hmGroup <- colData(gvst)[c("condition", "time_point")] #Extracts the condition and time_point columns from object gvst

hmGroup

hmGroup$time_point <- as.character(hmGroup$time_point)  #Converts the time_point column to character type.
# Ensure correct order of timepoints
hmGroup$time_point <- factor(hmGroup$time_point, levels = c("Day_20", "Day_182", "Day_313", "Day_445", "Day_738"))
# Ensure correct order of conditions
hmGroup$condition <- factor(hmGroup$condition, levels = c("5FBS", "10FBS", "15FBS"))

#To see what the levels of "condition" are
hmGroup@listData[["condition"]]

Conditioncol = c("5FBS"="#1f77b4", "10FBS" ="#ffdd57", "15FBS" = "#d62728") #define color mappings for the condition and time_point annotations.
Timepointcol = c(
  "Day_20" = "#e5f5e0",  # Light green
  "Day_182" = "#a1d99b",  # Light-medium green
  "Day_313" = "#74c476",  # Medium green
  "Day_445" = "#31a354",  # Dark-medium green
  "Day_738" = "#006d2c"  # Dark green
)

hm_ann_col <- HeatmapAnnotation(df=hmGroup, #gp = gpar(col = "grey"),
                                col = list(condition=Conditioncol, time_point=Timepointcol)) # creates a HeatmapAnnotation object using the hmGroup data frame and the defined color mappings. Condition (FBS) treatment on top, Timepoint below.



## Convert the data frame to a matrix
gcount_filt_vst_centered_matrix <- as.matrix(gcount_filt_vst_centered)


#First heatmap try with 3 clusters

#Do NOT let ComplexHeatmap compute it's own k-means clustering-- use the assignments you calculated above using kmeans(). 

DEGheatmap <-  Heatmap(gcount_filt_vst_centered_matrix, column_title = "", 
                       use_raster = TRUE,
                       name = "expression",
                       show_row_names = FALSE, 
                       top_annotation = hm_ann_col, 
                       show_column_names = TRUE, 
                       row_dend_side = "left" ,
                       column_split = 3, column_dend_height = unit(0.2, "in"),
                       cluster_row_slices = FALSE, #disable reordering, ensure cluster1 is on top
                       row_split = gene_clusters, #use this if you want to cluster based on kmeans() assignments above
                       #km=3, #use this if you want ComplexHeatmap() to do it's own kmeans clustering
                       row_km_repeats = 100,
                       row_gap = unit(2.5, "mm"), border = TRUE,
                       column_names_gp =  gpar(fontsize = 10))


print(DEGheatmap)
dev.off()


#How many genes did Heatmap put into each of the clusters?
# Extract the row clustering (k-means) information
row_clusters <- row_order(DEGheatmap)


# Count how many genes are in each cluster
cluster_counts <- sapply(row_clusters, length)
cluster_counts
Cluster1 Cluster2 Cluster3 
1466     922     1139

# Total number of genes across all clusters
total_genes <- sum(cluster_counts)

# View the result
total_genes #3,527 genes in heatmap


### List of individual clusters 
cluster1 <- cluster_res %>% filter(cluster == 1) 
cluster2 <- cluster_res %>% filter(cluster == 2)
cluster3 <- cluster_res %>% filter(cluster == 3)


nrow(cluster1) #1466
nrow(cluster2) #922
nrow(cluster3) #1139

#Save lists of gene IDs assigned to each cluster
#write.csv(cluster1, file = 'Cluster1_ID_list.csv')
#write.csv(cluster2, file = 'Cluster2_ID_list.csv')
#write.csv(cluster3, file = 'Cluster3_ID_list.csv')


## Heatmap WITHOUT cluster2
######### Try heatmap without Cluster3 genes (strange peak at 313 days; 1,139 genes total)
#Need to modify gcount_filt_vst_centered_matrix to drop all cluster3 genes
#Extract cluster3 gene IDs
cluster3_genes <- rownames(cluster3)
length(cluster3_genes)

#Find genes in are in matrix but NOT in cluster3
genes_to_keep <- setdiff(rownames(gcount_filt_vst_centered_matrix), cluster3_genes)
length(genes_to_keep) #should be 2,388 total

#Subset matrix to keep only desired genes (those in cluster 1 and 2)
gcount_filt_vst_centered_matrix_cluster3_removed <- gcount_filt_vst_centered_matrix[genes_to_keep, ]

nrow(gcount_filt_vst_centered_matrix_cluster3_removed) #2388, correct

# Filter out Cluster 3 from your cluster assignments
cluster_res_filtered <- cluster_res[cluster_res$cluster != 3, ]

# Match the filtered cluster assignments to the filtered matrix
gene_clusters_new <- cluster_res_filtered[match(rownames(gcount_filt_vst_centered_matrix_cluster3_removed), 
                                                cluster_res_filtered$gene_id), "cluster"]

#Convert to factor
gene_clusters_new <- factor(gene_clusters_new, levels = c("Cluster1", "Cluster2"))




### Timepoint ordering###

# Combine condition and time_point into a single column for sorting
hmGroup$sample_id <- paste(hmGroup$time_point, hmGroup$condition, sep = "_")

# Create a desired order of sample IDs
desired_order <- with(hmGroup, order(time_point, condition))

# Reorder hmGroup and the matrix columns
hmGroup <- hmGroup[desired_order, ]
gcount_filt_vst_centered_matrix_cluster3_removed <- gcount_filt_vst_centered_matrix_cluster3_removed[, desired_order]

# Recreate the column annotation with the reordered hmGroup
hm_ann_col <- HeatmapAnnotation(df = hmGroup,
                                col = list(condition = Conditioncol, time_point = Timepointcol))



# Plot new heatmap
DEGheatmap_no_cluster3 <- Heatmap(gcount_filt_vst_centered_matrix_cluster3_removed, 
                                  column_title = "", 
                                  use_raster = TRUE,
                                  name = "expression",
                                  show_row_names = FALSE, 
                                  top_annotation = hm_ann_col, 
                                  show_column_names = TRUE, 
                                  row_dend_side = "left",
                                  #column_split = 5, 
                                  cluster_columns = FALSE, #FALSE to remove top dendrogram
                                  column_dend_height = unit(0.2, "in"),
                                  row_split = gene_clusters_new,
                                  cluster_row_slices = FALSE, #preserves factor levels
                                  #km=2,
                                  row_km_repeats = 100, 
                                  row_gap = unit(2.5, "mm"), 
                                  border = TRUE,
                                  column_names_gp = gpar(fontsize = 10))

print(DEGheatmap_no_cluster3)




#How many genes did Heatmap put into each of the clusters?
# Extract the row clustering (k-means) information
row_clusters <- row_order(DEGheatmap_no_cluster3)


# Count how many genes are in each cluster
cluster_counts <- sapply(row_clusters, length)
cluster_counts
Cluster1 Cluster2
1466    922

# Total number of genes across all clusters
total_genes <- sum(cluster_counts)

# View the result
total_genes #2,388 genes in heatmap

#Save lists of gene IDs assigned to each cluster
write.csv(cluster1, file = 'Cluster1_ID_list.csv')
write.csv(cluster2, file = 'Cluster2_ID_list.csv')






#####




########Convert Cluster1 and Cluster2 IDs to protein FASTAS #############################
#3_17_25 Retrieve protein sequences FASTA given an input of Lvariegatus LOC gene IDs.
#Biostrings is used to read/handle FASTA files.

library(Biostrings)

################################## CLUSTER 1 ##########################
#Read the protein FASTA into R
fasta_file <- "Lvar_3.0.pep.all_LOCheaders.fasta"
protein_sequences <- readAAStringSet(fasta_file)

#Original protein FASTA was downloaded from Ensembl http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/lytechinus_variegatus_gca018143015v1/. I processed it through my python script "SeaUrchin Protein FASTA Manipulator convert headers to LOC IDs" to convert headers to LOC IDs. The final FASTA was saved as Lvar_3.0.pep.all_LOCheaders.fasta.

#Load Your List of LOC IDs: Assuming your list of LOC IDs is in a text file or a CSV, you can read it into R using read.table or read.csv.
Cluster1_loc_ids <- read.csv("Cluster1_ID_list.csv", header = TRUE, stringsAsFactors = FALSE)
Cluster1_loc_ids <- Cluster1_loc_ids$gene  # Assuming the LOC IDs are in column named "gene"
Cluster1_loc_ids
summary(Cluster1_loc_ids) #should be 1466 total

#Match LOC IDs and Extract Sequences: Match the LOC IDs with the names in your protein_sequences object and extract the corresponding sequences.
Cluster1_matched_sequences <- protein_sequences[names(protein_sequences) %in% Cluster1_loc_ids]
summary(Cluster1_matched_sequences) #should be 1466 total.
#Find unmatched IDs
Cluster1_unmatched_ids <- Cluster1_loc_ids[!Cluster1_loc_ids %in% names(protein_sequences)]
# Print the unmatched IDs
print(Cluster1_unmatched_ids) #86 of them. why so many? quick check with NCBI: many are ncRNA.
summary(Cluster1_unmatched_ids)

#Input the unmatched IDs to echinobase, try to find protein sequence, add it to the reference protein FASTA. Run matching/protein seq pull again.

#Save the matched sequences AND unmatched sequences to a new FASTA file
writeXStringSet(Cluster1_matched_sequences, "Cluster1_proteinFASTA.fasta")
write.csv(Cluster1_unmatched_ids, file = 'Cluster1_NCRNAs.csv')

################################## CLUSTER 2 ##########################
#Read the protein FASTA into R
fasta_file <- "Lvar_3.0.pep.all_LOCheaders.fasta"
protein_sequences <- readAAStringSet(fasta_file)
#Original protein FASTA was downloaded from Ensembl http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/lytechinus_variegatus_gca018143015v1/. I processed it through my python script "SeaUrchin Protein FASTA Manipulator convert headers to LOC IDs" to convert headers to LOC IDs. The final FASTA was saved as Lvar_3.0.pep.all_LOCheaders.fasta.

#Load Your List of LOC IDs: Assuming your list of LOC IDs is in a text file or a CSV, you can read it into R using read.table or read.csv.
Cluster2_loc_ids <- read.csv("Cluster2_ID_list.csv", header = TRUE, stringsAsFactors = FALSE)
Cluster2_loc_ids <- Cluster2_loc_ids$gene  # Assuming the LOC IDs are in column named "gene"
Cluster2_loc_ids
summary(Cluster2_loc_ids) #should be 922 total

#Match LOC IDs and Extract Sequences: Match the LOC IDs with the names in your protein_sequences object and extract the corresponding sequences.
Cluster2_matched_sequences <- protein_sequences[names(protein_sequences) %in% Cluster2_loc_ids]
summary(Cluster2_matched_sequences) #836 total, so many non-protein coding.
#Find unmatched IDs
Cluster2_unmatched_ids <- Cluster2_loc_ids[!Cluster2_loc_ids %in% names(protein_sequences)]
# Print the unmatched IDs
print(Cluster2_unmatched_ids) #252 of them. why so many? quick check with NCBI: many are ncRNA.
summary(Cluster2_unmatched_ids)

#Input the unmatched IDs to echinobase, try to find protein sequence, add it to the reference protein FASTA. Run matching/protein seq pull again.

#Save the matched sequences AND unmatched sequences to a new FASTA file
writeXStringSet(Cluster2_matched_sequences, "Cluster2_proteinFASTA.fasta")
write.csv(Cluster2_unmatched_ids, file = 'Cluster2_NCRNAs.csv')



####### MDS plot using the VST-transformed TPM count matrix gcount_filt_vst_centered_matrix 
library(edgeR)
packageVersion("edgeR")

#Define shapes to signify condition (5%FBS= black square, 10%FBS= black circle, 15%FBS=black triangle)
Shape<- c(15,15,15,15,15,18,18,18,18,18,17,17,17,17,17) 

#Define colors for timepoints, where red shows day 20 samples, blue shows day 182 samples, green shows day 313 samples, pink shows day 445 samples, and purple shows day 738 samples)
Color<- c("red", "blue", "green", "magenta", "purple","red", "blue", "green", "magenta", "purple", "red", "blue", "green", "magenta", "purple") #colors=days in culture 20, 182, 313,445 and 738


#Plot MDS
plotMDS(gcount_filt_vst_centered_matrix, pch=Shape, col=Color, cex=2, main="MDS Plot of TPM-VST RNAseq Transcript Counts")

#Add legend for FBS Conditions
legend("topleft", inset=c(0,0), pch=c(15,18,17), legend=c("5% FBS","10% FBS", "15% FBS"), col="black", bty="n", pt.cex=1.5)
#pch sets the shape of the point
#inset dictates where the legend appears
```
#pt.cex increases the size of the point

#Add legend for timepoints
legend("topleft", inset=c(0.21,0), pch=16, legend=c("Day 20","Day 182", "Day 313", "Day 445", "Day 738"), col=c("red", "blue", "green", "magenta", "purple"), bty="n", pt.cex=1.5)
