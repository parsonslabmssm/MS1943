###A & E are DMSO treated and BCD & FGH are degrader treated
##install different package that we will use below
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
# biocLite("readr")
# biocLite("apeglm")
# biocLite("pheatmap")
# biocLite("ReportingTools")
# biocLite("gage")
# biocLite("pathview")
# biocLite("org.Hs.eg.db")
# biocLite("gageData")
# biocLite("goseq")
# biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
# biocLite("DESeq2")
# biocLite("GenomicFeatures")
# biocLite("biomaRt")
# biocLite("vsn")

# doc http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#load package use to analyse reads count from Salmon analysis
library("tximport")
library("readr")
library("GenomicFeatures")

#director where data from Salmon are saved
dir <- "XXXXXX/results/RNAseq"
samples <- read.table(file.path(dir,"meta_EZH2.csv"), header=TRUE,sep=",")
samples$treated <- as.factor(samples$treated)
files <- file.path(dir,"salmon", paste0(samples$run, "_quant.sf"))

#Make A TxDb Object From Annotations Available in GTF File from TCGA 
gffFile <- paste0(dir,"/reference/GDC.GRCh38_gencode.v22.annotation.gtf")

#available.species() ## to know how to write the name of organism
#https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
txdb <- makeTxDbFromGFF(file=gffFile, format="gtf", dataSource="TCGA",
                        organism="Homo sapiens")
#list of chromosome
seqlevels(txdb)
#colnomes of data
columns(txdb)
#list of transcrip
GR <- transcripts(txdb)

#reading this tx2gene data.frame can be accomplished from a TxDb object
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys =k, keytype ="GENEID", columns ="TXNAME")
#import the necessary quantification data for DESeq2 using the tximport function
names(files) <- paste0("sample",1:4)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

#construct a DESeqDataSet from the txi object and sample information in samples
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ treated)
dim(ddsTxi)
##dim: 59979 4 

##Pre-filtering
#a minimal pre-filtering to keep only rows that have at least 10 reads total. 
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
dim(dds)
# now dim: 24924 4 

#using relevel, just specifying the reference level
dds$treated <- relevel(dds$treated, ref = "DMSO")

##Differential expression analysis
dds_desq2 <- DESeq(dds)
res <- results(dds_desq2, alpha=0.05)
res

#How many adjusted p-values were less than 0.05?
res2 <- results(dds_desq2, contrast=c("treated","degrader","DMSO"), alpha=0.05)

#reorder by pvalue (ascendant) and log2FoldChange (descendant)
res2Ordered <- res2[order(res2$pvalue,-abs(res2$log2FoldChange)),]
sum(res2$padj < 0.05, na.rm=TRUE)
summary(res2)

#Exploring and exporting results
#MA-plot
plotMA(res2Ordered)

#identify gene
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

#Plot counts
plotCounts(dds_desq2, gene=which.min(res2Ordered$padj), intgroup="treated")
plotCounts(dds_desq2, gene="ENSG00000116774.10", intgroup="treated")

d <- plotCounts(dds_desq2, gene="ENSG00000116774.10", intgroup="treated", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=treated, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#More information on results columns
mcols(res2Ordered)$description

# find the Gene name of genes in ENSEMBL
#connexion to ENSEMBL
library("biomaRt")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
res2Ordered$ensembl_gene_id<- gsub("\\..*","",rownames(res2Ordered))
res2Ordered$ensembl_gene_idnum<- rownames(res2Ordered)

#Extract the list of hgcn_symbol and gene_biotype for each ensembl_gene_id in my list
hgnc <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype',"entrezgene"),
              filters = 'ensembl_gene_id', values = res2Ordered$ensembl_gene_id,
              mart = mart)
res2Ordered_ensembl <-merge(as.data.frame(res2Ordered),hgnc,by='ensembl_gene_id')
res2Ordered_ensembl <- res2Ordered_ensembl[order(res2Ordered_ensembl$pvalue,-abs(res2Ordered_ensembl$log2FoldChange)),]

##Write the final file
write.csv(as.data.frame(res2Ordered_ensembl), 
          file=paste0(dir,"/MS1943_DGE_results_ENSEMBL.csv"))
