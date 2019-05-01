# MS1943
Script associated with the paper "Discovery of a first-in-class EZH2 selective degrader"

We produced RNA-seq data for 4 samples ( MDA-MD-468 breast cancer cell lines): 2 under DMSO and 2 under MS1943 degrader.

1. Quality control of Fastq files to reads count 
After quality control of FASTQ files using FASTQC tool (version 0.11.7) (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), we trimmed low quality bases (PHRED < 10) and adapter sequences and then discarded short reads (length < 60nt) using bbduk tool (version 37.53) (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). When either forward or reserve of a paired-end read are discarded, we discarded the complete paired-end read. We quantified the expression of transcripts using cleaned paired-end reads with Salmon tool (version 0.9.1) using human reference transcriptome from The Cancer Genome Atlas (TCGA) (GDC.h38 GENCODE v22). 


To run the pipeline, please update the absolute paths and use the command line on LSF scheduler:

bsub < RNAseq_LSFjob.sh

if you have another scheduler, you just need to adapt the script RNAseq_LSFjob.sh to your scheduler.
