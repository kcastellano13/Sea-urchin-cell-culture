#!/bin/bash
# 10x cellranger make reference - L variegatus genome, transcriptome and gene annotation file
# Author: Kate Castellano
# ----------------------------

#make v3.0 L variegatus genome into a 10x index
#/data/app/cellranger-6.1.2/cellranger mkref \
#--genome=Lv_genome_10xindx \
#--fasta=GCF_018143015.1_Lvar_3.0_genomic.fna \
#--gene=GCF_018143015.1_Lvar_3.0_genomic.gtf 2>&1 | tee -a cellranger_mkref_log.txt


#make v3.0 L variegatus genome into a 10x index - redo after filtering GTF file
#/data/app/cellranger-6.1.2/cellranger mkref \
#--genome=Lv_genome_filteredGenes_10xindx \
#--fasta=GCF_018143015.1_Lvar_3.0_genomic.fna \
#--gene=GCF_018143015.1_Lvar_3.0_genomic.filter.gtf 2>&1 | tee -a cellranger_mkref_filterGenes__log.txt


#make v3.0 L variegatus genome into a 10x index - redo after adding gene_name attribute and filtering GTF file
/data/app/cellranger-6.1.2/cellranger mkref \
--genome=Lv_genome_filteredGenes_addGeneName_10xindx \
--fasta=GCF_018143015.1_Lvar_3.0_genomic.fna \
--gene=GCF_018143015.1_Lvar_3.0_genomic_addGeneName.filter.gtf \
2>&1 | tee -a cellranger_mkref_filterGenes_addGeneName__log.txt
