#Data from R/DGE/Salmon2DGE.R

##### data transformation and visualisation
#Extracting transformed values
vsd <- vst(dds_desq2, blind=FALSE)
rld <- rlog(dds_desq2, blind=FALSE)
head(assay(dds_desq2), 3)

#Effects of transformations on the variance
ntd <- normTransform(dds_desq2)
library("vsn")
meanSdPlot(assay(dds_desq2))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


##Heatmap of the count matrix
library("pheatmap")

resSigLFC <- subset(resLFCOrdered, padj < 0.05 )
nrow(resSigLFC) ## number of significant


select <- order(resLFC$pvalue,-abs(resLFC$log2FoldChange))[1:nrow(resSigLFC)]
df <- as.data.frame(colData(dds_desq2)[,c("treated","cellline")])

#Plot heatmap
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
