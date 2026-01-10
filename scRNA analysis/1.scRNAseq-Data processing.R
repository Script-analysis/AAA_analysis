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
###################
samples=list.files()
samples=samples[1:4]
for (i in seq_along(samples)[1:4]){
  assign(paste0("scRNA", samples[i]), Read10X(data.dir = samples[i]))
}
for (i in seq_along(samples)[1:4]){
  assign(paste0("scRNA", samples[i]), CreateSeuratObject(counts = eval(parse(text = paste0("scRNA", samples[i]))), project = samples[i]))
}
scRNA_AAA <- merge(scRNAAAA1, y = c(scRNAAAA2,scRNAAAA3,scRNAAAA4), add.cell.ids = samples[1:4], project = "AAA")
library(stringr)
meta<-scRNA_AAA@meta.data
scRNA_AAA$orig.ident
orig=as.data.frame((str_split(rownames(scRNA_AAA@meta.data),pattern = '_')))
orig=as.data.frame(t(orig))
scRNA_AAA@meta.data$sample=orig$V1
scRNA_AAA@meta.data$group="AAA"
samples=list.files()
######################################Ctrl
for (i in seq_along(samples)[5:6]){
  assign(paste0("scRNA", samples[i]), Read10X(data.dir = samples[i]))
}
for (i in seq_along(samples)[5:6]){
  assign(paste0("scRNA", samples[i]), CreateSeuratObject(counts = eval(parse(text = paste0("scRNA", samples[i]))), project = samples[i]))
}

scRNA_Ctrl <- merge(scRNACtrl1, y = c(scRNACtrl2), add.cell.ids = samples[5:6], project = "Ctrl")

library(stringr)
meta<-scRNA_Ctrl@meta.data
View(meta)
scRNA_Ctrl$orig.ident
orig=as.data.frame((str_split(rownames(scRNA_Ctrl@meta.data),pattern = '_')))
orig=as.data.frame(t(orig))
scRNA_Ctrl@meta.data$sample=orig$V1
scRNA_Ctrl@meta.data$group="Ctrl"
scRNA <- merge(scRNA_Ctrl, y = c(scRNA_AAA), add.cell.ids = c('Ctrl','AAA') , project = "CtrlvsAAA")
#####################QC
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
mito_genes <- grep("^[M,m][T,t]-", rownames(scRNA), value =TRUE)
scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern ="^[M,m][T,t]-")
scRNA[["percent.hb"]] <- PercentageFeatureSet(scRNA,  pattern = "^HB")
ribo_genes <- grep("^R[P,p][SL,sl]", rownames(scRNA), value =TRUE)
scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern ="^R[P,p][SL,sl]")

VlnPlot(scRNA, group.by = "sample",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"), 
        cols = mycolor2,pt.size = 0, 
        ncol = 5)   + theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
                            legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size = 12),
                            axis.text.x = element_text(color="black",size=12),
                            axis.text.y = element_text(color="black",size=12),
                            axis.title.x = element_text(face="plain", color="black",size=12),
                            axis.title.y = element_text(face="plain", color="black",size=12))

ggsave("vlnplot_before_qc_sample.pdf", width = 20, height = 5) 
#####################
minGene=200
maxGene=4500
pctMT=50
pctRB= 20
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT& percent.ribo < pctRB)
VlnPlot(scRNA, group.by = "sample",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"), 
        cols = mycolor2,
       pt.size = 0, 
        ncol = 5) + 
  theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size = 12),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_text(color="black",size=12),
        axis.title.x = element_text(face="plain", color="black",size=12),
        axis.title.y = element_text(face="plain", color="black",size=12))
ggsave("vlnplot_after_qc_sample.pdf", width = 20, height = 5) 
