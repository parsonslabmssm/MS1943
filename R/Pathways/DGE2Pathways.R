#Data from R/DGE/Salmon2DGE.R

###GSEA
deseq2.fc=res2$log2FoldChange
names(deseq2.fc)=rownames(res2)

## convert ENSEMBL 2 ENTREZ
res2_ensembl <- names(deseq2.fc)
res2_ensembl_short<- gsub("\\..*","",res2_ensembl)
length(res2_ensembl_short)
library("org.Hs.eg.db")
# Gets the entrez gene IDs for thea Ensembl IDs
a <- sapply(res2_ensembl_short, function(x) exists(x, org.Hs.egENSEMBL2EG))
res2_ensembl.existed <- res2_ensembl_short[a]
xx <- as.list(org.Hs.egENSEMBL2EG)
list_names<-xx[res2_ensembl.existed]

entrezID_names <- unlist(list_names)
entrezID2Ensembl_names <- cbind(entrezID_names,names(entrezID_names))
colnames(entrezID2Ensembl_names) <- c("entrezId","EnsembID_gene")
dim(entrezID2Ensembl_names)
names_deseq2.fc <- names(deseq2.fc)
names_res2 <- cbind(names_deseq2.fc,res2_ensembl_short,seq(1,length(names_deseq2.fc),1))
colnames(names_res2)<- c("EnsembID_geneVersion","EnsembID_gene","rang")
dim(names_res2)
names_deseq2.fc_entreid <- merge(names_res2, entrezID2Ensembl_names, by="EnsembID_gene",all.x=TRUE)
dim(names_deseq2.fc_entreid)
names_deseq2.fc_entreid$finalname <- names_deseq2.fc_entreid$entrezId
names_deseq2.fc_entreid$finalname <- as.character(names_deseq2.fc_entreid$finalname)
names_deseq2.fc_entreid$finalname[which(is.na(names_deseq2.fc_entreid$finalname))] <- names_deseq2.fc_entreid$EnsembID_gene[which(is.na(names_deseq2.fc_entreid$finalname))]

head(names_deseq2.fc_entreid)
head(deseq2.fc)
names(deseq2.fc) <- names_deseq2.fc_entreid$finalname
head(deseq2.fc)
exp.fc=deseq2.fc
out.suffix="deseq2"

require(gage)
##KEGG
data(kegg.gs)
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.05 &
  + !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.05 &
  + !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]

kegg_goseq_sel <- which(grepl("3030",rownames(fc.kegg.p$greater)))
path.ids.goseq <- rownames(fc.kegg.p$greater)[kegg_goseq_sel]

path.ids2 <- substr(c(path.ids, path.ids.l,path.ids.goseq), 1, 8)
require(pathview)
#view different pathways
pv.out.list <- sapply(path.ids2, function(pid) 
  pathview(gene.data = exp.fc, pathway.id = pid,
           species = "Homo sapiens", out.suffix=out.suffix))

#### GOseq
## significant
library("goseq")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
res2_sign <- res2[which(res2$padj <0.05 & res2$log2FoldChange !=0),]
dim(res2_sign)

genes=as.integer(p.adjust(res2$pvalue[res2$log2FoldChange!=0],
                          method="BH")<.05)

names(genes)=gsub("\\..*","",row.names(res2[res2$log2FoldChange!=0,]))
table(genes)

## KEGG
# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Hs.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg=as.list(org.Hs.egPATH)
# Define a function which gets all unique KEGG IDs
# associated with a set of Entrez IDs
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)
head(kegg)

#Create a table between Ensembl_gene_id and kegg_id
kegg_df <- c("ensembl_gene_id","kegg_id")

for (i in 1:length(kegg)){
  names_ensembl <- names(kegg[i])
  num_kegg <-length(kegg[[i]])
  # cat(num_kegg)
  # cat("\n")
  for( j in 1:num_kegg){
    # cat(names_ensembl)
    #cat(" ")
    name_kegg <- kegg[[i]][j]
    new_elt <- c(names_ensembl,name_kegg)
    #cat(paste(names_ensembl,name_kegg,"\n"))
    kegg_df <- rbind(kegg_df,new_elt)
  }
}
colnames(kegg_df)<-kegg_df[1,]
kegg_df <- kegg_df[-1,]
fileKEGG_df="/Users/tiphainecmartin/Documents/Projects_MSSM/ParsonsLab/Elias/github/parsonslab/EZH2degrader/results/RNAseq/enrichment/EnsemblGeneID2KEGG.csv"
write.table(kegg_df,file=fileKEGG_df,row.names =FALSE,sep=",")

pwf=nullp(genes,"hg38","ensGene")
KEGG=goseq(pwf,gene2cat=kegg)
head(KEGG)
KEGG$padj <- p.adjust(KEGG$over_represented_pvalue,method="BH")
enriched.KEGG=KEGG[KEGG$padj<.05,]
head(enriched.KEGG)
fileEnrichKEGG="/Users/tiphainecmartin/Documents/Projects_MSSM/ParsonsLab/Elias/github/parsonslab/EZH2degrader/results/RNAseq/enrichment/degraderEZH2_EnrichKEGG.csv"
write.table(enriched.KEGG,file=fileEnrichKEGG,row.names =FALSE,sep=",")

##list of genes for category_significant
num_kegg_sign <- nrow(enriched.KEGG)
ensembl_keggsignif<-c("ensembl_gene_id","kegg_id")
for(i in 1:num_kegg_sign){
  category_selected <- enriched.KEGG[i,"category"]
  data <- kegg_df[which(kegg_df[,"kegg_id"] %in% category_selected),]
  ensembl_keggsignif <- rbind(ensembl_keggsignif,data)
}
colnames(ensembl_keggsignif)<- ensembl_keggsignif[1,]
ensembl_keggsignif <- ensembl_keggsignif[-1,]


hgnc_kegg <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype',"entrezgene"),
                   filters = 'ensembl_gene_id', values = ensembl_keggsignif[,"ensembl_gene_id"],
                   mart = mart)

ensembl_keggsignif_annot <-merge(as.data.frame(ensembl_keggsignif),hgnc_kegg,by='ensembl_gene_id')

##ENSEMBL in significantpadj <0.05
ensembl_keggsignif_annot$DESeq2_all <- 0
gene_names <- gsub("\\..*","",row.names(res2[which(res2$padj <0.05 & res2$log2FoldChange!=0),]))
ensembl_keggsignif_annot[which(ensembl_keggsignif_annot$ensembl_gene_id %in% gene_names),
                         "DESeq2_all"] <- 1

ensembl_keggsignif_annot$DESeq2_lfc1 <- 0

##ENSEMBL in significantpadj <0.05 & log2FoldChange >1
gene_names2 <- gsub("\\..*","",row.names(res2[which(res2$padj <0.05 
                                                    & abs(res2$log2FoldChange) > 1),]))
ensembl_keggsignif_annot[which(ensembl_keggsignif_annot$ensembl_gene_id %in% gene_names2),
                         "DESeq2_lfc1"] <- 1
length(which(ensembl_keggsignif_annot$DESeq2_all==1))
length(which(ensembl_keggsignif_annot$DESeq2_lfc1==1))
ensembl_keggsignif_annot[which(ensembl_keggsignif_annot$DESeq2_lfc1==1),]
ensembl_keggsignif_annot[which(ensembl_keggsignif_annot$DESeq2_all==1),]
fileEnsemblEnrichKEGG="XXXX/results/RNAseq/enrichment/degraderEZH2_EnrichKEGG_ensemblannotated.csv"
write.table(ensembl_keggsignif_annot,file=fileEnsemblEnrichKEGG,row.names =FALSE,sep=",")
