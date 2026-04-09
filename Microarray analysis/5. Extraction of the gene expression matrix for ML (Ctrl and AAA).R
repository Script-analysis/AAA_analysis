setwd("C:\\AAA_analysis\\Control_AAA\\step1.ID")
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(sva)
genes="DMR_DEGs_18genes.txt"
genes<-read.table(genes, header=T, sep="\t", check.names=F)
file<-"GSE7084_normalizationt.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
genes <- genes$ID
genes_use <- intersect(genes, rownames(data))
data_sub <- data[genes_use, ,drop = FALSE]
write.table(data_sub,file="GSE7084_18genes.txt",sep="\t",quote=F,col.names=T,row.names = T)
file<-"GSE57691_normalizationt.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
genes <- genes$ID
genes_use <- intersect(genes, rownames(data))
data_sub <- data[genes_use, ,drop = FALSE]
write.table(data_sub,file="GSE57691_18genes.txt",sep="\t",quote=F,col.names=T,row.names = T)