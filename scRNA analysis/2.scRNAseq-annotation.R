rm(list=ls())
library(dplyr)
library(MAST)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc)
library(biomaRt)
library(Seurat)
library(SingleCellExperiment)
library(ggrepel)
library(graphics)
library(Matrix)
library(lubridate)
library(glue)
library(gridExtra)
library(gdata)
library(slingshot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(gridExtra)
library(SummarizedExperiment)
library(DoubletFinder)
library(monocle)
#####################################
scRNA<- NormalizeData(scRNA)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
library(colourpicker)
scRNA<-ScaleData(scRNA,verbose = FALSE)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object =scRNA),reduction.name = "pca")
scRNA<- RunUMAP(scRNA,reduction = "pca", dims = 1:10, reduction.name = "umap_naive")
library(harmony)
scRNA<- RunHarmony(scRNA,reduction = "pca",group.by.vars = "sample",reduction.save = "harmony")
scRNA<- RunUMAP(scRNA, reduction = "harmony", dims = 1:30,reduction.name = "umap")
scRNA <- FindNeighbors(scRNA,reduction = "harmony", dims = 1:25)
scRNA<- FindClusters(scRNA, resolution = seq(0,1.2,by=0.05))
library(clustree)
clustree(scRNA)
ggsave("clustree.pdf",height=17,width=15)
scRNA<- FindClusters(scRNA,resolution = 0.6)
scRNA = RunTSNE(scRNA, dims = 1:25)
embed_tsne <- Embeddings(scRNA, 'tsne')  
write.csv(embed_tsne,'embed_tsne.csv')
scRNA <- RunUMAP(scRNA, dims =1:25)
embed_umap <- Embeddings(scRNA, 'umap') 
write.csv(embed_umap,'embed_umap.csv') 


library(Seurat)
library(ggplot2)
library(ggunchull)
library(ggrepel) 
library(tidyverse)

pdf("./umap.pdf",width = 5,height = 4.5)
plot4= DimPlot(scRNA, reduction = "umap", pt.size = 1.5,cols=mycolors)
print(plot4)
dev.off()

pdf("./umap_splitbygroup.pdf",width = 10,height = 4.5)
plot4= DimPlot(scRNA, reduction = "umap", split.by="group",pt.size = 1.5,cols=mycolors)
print(plot4)
dev.off()

DefaultAssay(scRNA) <- "RNA"
all.markers  <- FindAllMarkers(scRNA, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75)
significant.markers  <- all.markers [all.markers$p_val_adj < 0.05, ]
write.csv(significant.markers, file = "significant.all.markers.csv")
##############################################
markers <- c('MYH11','TAGLN','ACTA2','DCN','COL1A1','COL3A1',
             'CDH5','PECAM1','FABP4',
             'CD68','CDCA3','C1QB',
             'CD79A','LY6D','CD79B',
             'CD3D','CD3G','CD28',
             "ENPP3","FCER1A",
             'CCR7','IFITM1')

library(ggplot2)
DotPlot(scRNA, features =markers ,assay='RNA') + coord_flip()

#####################
library(openxlsx)
annotation_results <- read.xlsx("annotation.xlsx")
print(head(annotation_results))
scRNA@meta.data$celltype <- "Unknown"  
for (i in seq_len(nrow(annotation_results))) {
  cluster_name <- annotation_results$cluster[i]
  cell_type_name <- annotation_results$celltype[i]

  scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters == cluster_name] <- cell_type_name
}
table(scRNA@meta.data$celltype)
Idents(scRNA)=scRNA@meta.data$celltype

pdf("./umap.pdf",width = 5,height = 4.5)
plot4= DimPlot(scRNA, reduction = "umap", pt.size = 1.5,cols=mycolors)
print(plot4)
dev.off()

pdf("./umap_splitbygroup.pdf",width = 10,height = 4.5)
plot4= DimPlot(scRNA, reduction = "umap", split.by="group",pt.size = 1.5,cols=mycolors)
print(plot4)
dev.off()
save(scRNA, file = "scRNA_annotation.RData")

