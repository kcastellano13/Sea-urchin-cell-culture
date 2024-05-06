#!/bin/bash
#Make the L. variegatus 3.0 genome into a hisat2 index
#K. Castellano
module load histat2/v2.2.1

hisat2-build GCA_018143015.1_Lvar_3.0_genomic.fna Lvar3.0_hisat2index

