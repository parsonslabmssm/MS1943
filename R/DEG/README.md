# MS1943
Script associated with the paper "Discovery of a first-in-class EZH2 selective degrader"

We produced RNA-seq data for 4 samples ( MDA-MD-468 breast cancer cell lines): 2 under DMSO and 2 under MS1943 degrader.

2. Differential Gene Expression between DMSO treatment and MS1943 degrader
We performed the differential expression between the two MDA-MD-468 cell lines treated with DMSO and the two MDA-MD-468 cell lines treated with EZH2 degrader at gene level using R tximport library (version 1.6.0) on R (version 3.4.3). We pre-filtered genes to keep only genes that have at least 10 reads in total (Number of genes = 24,924) and then performed differential gene expression using R DESeq2 library (version 1.81.1). We identified as differential expressed between two conditions when the p-value adjusted is at 5% and log2 fold change is more than 1. 

Scripts available under the folder "R/DGE"

