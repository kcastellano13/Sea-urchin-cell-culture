
# Single-cell sequencing

## 1. Demultiplex and generate fastq files from Illumina basecalls (demultiplexing Illumina base call files (BCL))
- program: cellranger-6.1.2
- script: 01_cellranger_mkfastq

## 2. Convert L. variegatus genome into a cellranger reference
- genome: L. variegatus v3.0 (https://www.ncbi.nlm.nih.gov/assembly/GCF_018143015.1)
- program: cellranger-6.1.2
- script: 02_cellranger_mkref.sh

## 3. Count reads per GEM using cellranger count
- This was run on each sample (3% FBS, 5% FBS, 10% FBS and 15% FBS) individually 
    - program:_ cellranger version 6.1.2
    - important flags:_ --expect cells flag = 10,000
    - script:_ 03_cellranger_count.sh

## 4. Differential Expression analyzed using Seurat v4.3.0 in R v.2023.09.1
- R scripts for each sample can be found above:
    - Lv_Seurat_3FBS.Rmd
    - Lv_Seurat_5FBS.Rmd
    - Lv_Seurat_10FBS.Rmd
    - Lv_Seurat_15FBS.Rmd
