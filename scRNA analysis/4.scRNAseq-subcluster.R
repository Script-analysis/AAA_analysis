library(CytoTRACE)
library(Seurat)
library(SeuratData)
library(reticulate)
library(ggplot2)
load("scRNA_annotation.RData")
Idents(scRNA) <- scRNA@meta.data$celltype
scRNA= scRNA[,scRNA@meta.data$celltype %in% c('BCs')]
##############################################Analysis of BCs subcluster
#The dimensionality reduction and annotation of BCs were performed using the same Seurat–Harmony workflow described in Step 2.
save(scRNA, file = "scRNA_BCs_subclusters_annotation.RData")
####################CytoTRACE and Pseudotime Analysis
#CytoTRACE
scRNA$celltype.group <- paste(scRNA$celltype, scRNA$group, sep = "_")
Idents(scRNA) <- "celltype.group"
phe=as.character(phe)
names(phe)<-rownames(scRNA@meta.data)
mat_3k<-as.matrix(scRNA@assays$RNA@counts)
results<-CytoTRACE(mat=mat_3k, ncores = 1)
plotCytoGenes(results, numOfGenes = 10, outputDir = "scRNA_cytotrace")
ggplot2::ggsave('results_plotCytoTRACE.pdf',width=7)
plotCytoTRACE(results, phenotype = phe, gene="FKBP11",outputDir = "FKBP11")
####################Monocle
#The code of Pseudotime Analysis was detailed in https://github.com/MaxMeieran/monocle2
