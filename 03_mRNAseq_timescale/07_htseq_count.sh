#!/bin/bash
#K. Castellano

module load htseq/v2.0.4

fileList="trim_10FBS_D445_S286 trim_10FBS_D738_S289 trim_15FBS_D182_S281 trim_15FBS_D20_S278 trim_15FBS_D313_S284 trim_15FBS_D445_S287 trim_15FBS_D738_S290 trim_5FBS_D182_S279 trim_5FBS_D20_S276 trim_5FBS_D313_S282 trim_5FBS_D445_S285 trim_5FBS_D738_S288 trim_10FBS_D20_S277 trim_10FBS_D182_S280 trim_10FBS_D313_S283"
INDIR="/data/prj/urchin/cell-culture/2023Lvar_CellCulture_RNAseqTimeCourse/trimmomatic/hisat_map2LvarGenome"
OUTDIR="/data/prj/urchin/cell-culture/2023Lvar_CellCulture_RNAseqTimeCourse/trimmomatic/htseq_counts/2024March_redo"
GTF="GCF_018143015.1_Lvar_3.0_genomic_exonsOnly.gtf"

for file in ${fileList}
do

#-r pos tells htseq-count that our BAM file is coordinate sorted.
#-f bam indicates that our input file is in BAM format.
htseq-count -r pos -f bam \
	${INDIR}/${file}.bam \
        ${GTF} > ${OUTDIR}/${file}.counts

done
