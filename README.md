# MS1943
Scripts associated with the paper "Discovery of a first-in-class EZH2 selective degrader", Ma, A. and Stratikopoulos, E. et al., under review, 2019

We produced RNA-seq data for 4 samples ( MDA-MD-468 breast cancer cell lines): 2 under DMSO and 2 under MS1943 degrader.

1. Quality control of Fastq files to reads count 
After quality control of FASTQ files using FASTQC tool (version 0.11.7) (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), we trimmed low quality bases (PHRED < 10) and adapter sequences and then discarded short reads (length < 60nt) using bbduk tool (version 37.53) (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). When either forward or reserve of a paired-end read are discarded, we discarded the complete paired-end read. We quantified the expression of transcripts using cleaned paired-end reads with Salmon tool (version 0.9.1) using human reference transcriptome from The Cancer Genome Atlas (TCGA) (GDC.h38 GENCODE v22). 

Scripts available under the folder "Script"

2. Differential Gene Expression between DMSO treatment and MS1943 degrader
We performed the differential expression between the two MDA-MD-468 cell lines treated with DMSO and the two MDA-MD-468 cell lines treated with EZH2 degrader at gene level using R tximport library (version 1.6.0) on R (version 3.4.3). We pre-filtered genes to keep only genes that have at least 10 reads in total (Number of genes = 24,924) and then performed differential gene expression using R DESeq2 library (version 1.81.1). We identified as differential expressed between two conditions when the p-value adjusted is at 5% and log2 fold change is more than 1. 

Scripts available under the folder "R/DGE"

3. Enrichment Pathways

We next performed gene set enrichment analysis that capture pathways perturbed towards both directions simultaneously using the 24,448 ranked genes identified in our dataset and annotated in ENSEMBL (version 94) against the KEGG pathways using R GAGE (version 2.28.2) and Pathview (version 1.18.2) libraries8 and the hallmarks gene set collection (MSigDB V6.2) using GSEA. 

Scripts available under the folder "R/Pathways"

4. Visualisation 
We produced several heatmaps and plots using R scripts.

Scripts available under the folder "R/visualisation"


### Developed by Dr. Tiphaine Martin - Icahn School of Medicine at Mount Sinai - Tisch Cancer Institute - https://github.com/TiphaineCMartin
