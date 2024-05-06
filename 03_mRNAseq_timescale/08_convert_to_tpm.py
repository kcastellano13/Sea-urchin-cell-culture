## python script to convert htseq counts to TPM
## RPKM,FPKM and TPM explained: https://www.youtube.com/watch?v=TTUrtCY2k-w
## How to calculate TPM with python (code below from this video): https://www.youtube.com/watch?v=WHvajzS9f8M

import numpy as np
import csv
import pandas as pd

#read in gene length data

length_dict = {}
with open("Lvar3.0_exonMerged_length.csv") as csvfile:

    reader = csv.reader(csvfile, delimiter=',')

    for row in reader:
        length_dict[row[0]] = float(row[1])

#read in RNA-seq raw count data
data = [] #empty list
with open("FBStimecourseAll_noZeros_TransposenoRowName.tsv") as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
        data.append(row)

genes = data[0] #says the first row are the gene names
data = data[1:] #says the data is from the second row down

#data is written as string so cast as a float (in numpy array)
data = np.array(data).astype(float)

print(data) #to check it looks good

## TPM Normalization
## Step 1: Normalize for gene length (kilobase) -> this gives reads per kilobase (RPK)
## Step 2: Normalize for sequencing depth (divide by total RPK)

#Step 1 - for loop to iterate through the rows and get rpk


rpk = np.zeros(shape=(data.shape[0],data.shape[1])) #array of zeros that is the same size as our data table

for i in range(data.shape[0]): #iterate through rows
    for j in range(data.shape[1]): #iterate through columns
        gene_length = length_dict[genes[j]] #get gene length

        gene_length_in_kilo = gene_length/1000 #some genes are long so get them down to numbers that are easier to work with (in kilobases which is standard for RPK)

        rpk[i,j] = data[i,j] / gene_length_in_kilo

#Step 2 
tpm = np.zeros(shape=(data.shape[0],data.shape[1])) #array of zeros that is the same size as our data table

for i in range(rpk.shape[0]): #iterate through rows (samples)
    total_rpk = np.sum(rpk[i,:]) #calculate total rpk (aka total sequencing depth)

    scaled_rpk = total_rpk / 1000000 #scale per million (these can both stay in the outer loop bc it will be the same value for all columns)

    for j in range(rpk.shape[1]): #iterate through columns

        if scaled_rpk != 0:
            tpm[i,j] = rpk[i,j] / scaled_rpk #If only this line is used, it will give you an error if any of the values are 0 so this is why we include the if else statement
        else:
            tpm[i,j] = rpk[i,j] #will just print 0 if the value is 0

pd.DataFrame(tpm).to_csv("tpm.csv", header=genes, index=False)

