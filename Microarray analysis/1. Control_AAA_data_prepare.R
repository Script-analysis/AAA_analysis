setwd("C:\\AAA_analysis\\Control_AAA\\step1.ID")
library(GEOquery)
library(dplyr)
library(stringr)
library(tidyverse)
gse<-getGEO('GSE7084',destdir = '.',AnnotGPL = T,getGPL = T) 
exprset<-exprs(gse[[1]])%>%as.data.frame()
fdata<-fData(gse[[1]])
GPL<-select(fdata,ID,`Gene symbol`)
exprset$ID=rownames(exprset)
expr<-left_join(GPL,exprset,by='ID')%>%select(-ID)
expr$`Gene symbol`<-data.frame(sapply(expr$`Gene symbol`,
                                      function(x)unlist(strsplit(x,"///"))[1]),
                               stringsAsFactors=F)[,1]
expr<-expr%>%distinct(expr[1],.keep_all =  T)%>%na.omit()%>%remove_rownames()%>%column_to_rownames(var='Gene symbol')
write.table(expr,
            file = "GSE7084_Raw.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
gse<-getGEO('GSE57691',destdir = '.',AnnotGPL = T,getGPL = T) 
exprset<-exprs(gse[[1]])%>%as.data.frame()
fdata<-fData(gse[[1]])
GPL<-select(fdata,ID,`Gene symbol`)
exprset$ID=rownames(exprset)
expr<-left_join(GPL,exprset,by='ID')%>%select(-ID)
expr$`Gene symbol`<-data.frame(sapply(expr$`Gene symbol`,
                                      function(x)unlist(strsplit(x,"///"))[1]),
                               stringsAsFactors=F)[,1]
expr<-expr%>%distinct(expr[1],.keep_all =  T)%>%na.omit()%>%remove_rownames()%>%column_to_rownames(var='Gene symbol')
write.table(expr,
            file = "GSE57691_Raw.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
#Extraction of expression matrix and clinical data
setwd("C:\\AAA_analysis\\Control_AAA\\step2.ID")
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(sva)
library(openxlsx)
clin<-read.xlsx("Meta.xlsx")
file<-"GSE7084_Raw.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
samples <- clin$Sample
samples_use <- intersect(samples, colnames(data))
data_sub <- data[, samples_use, drop = FALSE]
clin_sub <- clin[clin$Sample %in% samples_use, ]
clin_sub <- clin_sub[match(colnames(data_sub), clin_sub$Sample), ]
all(clin_sub$Sample == colnames(data_sub))
write.table(data_sub,file="Rawdata_GSE7084.txt",sep="\t",quote=F,col.names=T,row.names = T)
write.table(clin_sub,file="clin_GSE7084.txt",sep="\t",quote=F,col.names=T,row.names = T)
clin<-read.xlsx("Meta.xlsx")
file<-"GSE57691_Raw.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
samples <- clin$Sample
samples_use <- intersect(samples, colnames(data))
data_sub <- data[, samples_use, drop = FALSE]
clin_sub <- clin[clin$Sample %in% samples_use, ]
clin_sub <- clin_sub[match(colnames(data_sub), clin_sub$Sample), ]
all(clin_sub$Sample == colnames(data_sub))
write.table(data_sub,file="Rawdata_GSE57691.txt",sep="\t",quote=F,col.names=T,row.names = T)
write.table(clin_sub,file="clin_GSE57691.txt",sep="\t",quote=F,col.names=T,row.names = T)

#log2 transformation
file1="Rawdata_GSE7084.txt"
file2="Rawdata_GSE57691.txt"
clin1="clin_GSE7084.txt"
clin2="clin_GSE57691.txt"
expr1 <- read.table(file1, header=T, sep="\t", check.names=F,row.names = 1)
expr2 <- read.table(file2, header=T, sep="\t", check.names=F,row.names = 1)
meta1 <-read.table(clin1, header=T, sep="\t", check.names=F,row.names = 1)
meta2 <-read.table(clin2, header=T, sep="\t", check.names=F,row.names = 1)
meta1$dataset <- "array1"
meta2$dataset <- "array2"
sum(colnames(expr1) %in% meta1$Sample)
sum(colnames(expr2) %in% meta2$Sample)
check_need_log2 <- function(mat) {
  qx <- as.numeric(quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  print(qx)
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    mat[mat <= 0] <- NA
    mat <- log2(mat)
    message("log2 transform applied")
  } else {
    message("log2 transform not needed")
  }
  return(mat)
}
expr1 <- check_need_log2(expr1)
expr2 <- check_need_log2(expr2)
expr1 <- expr1[rowSums(is.na(expr1)) == 0, ]
expr2 <- expr2[rowSums(is.na(expr2)) == 0, ]
expr1_norm <- normalizeBetweenArrays(expr1, method = "quantile")
expr2_norm <- normalizeBetweenArrays(expr2, method = "quantile")
meta1 <- meta1[match(colnames(expr1_norm), meta1$Sample), ]
meta2 <- meta2[match(colnames(expr2_norm), meta2$Sample), ]
pdf("Normalization.pdf",height=3,width=12)
par(mfrow = c(1, 2)) 
boxplot(expr1_norm, outline = FALSE, las = 2, main = "Array1 normalized")
boxplot(expr2_norm, outline = FALSE, las = 2, main = "Array2 normalized")
dev.off()
write.table(expr1_norm,file="GSE7084_normalization.txt",sep="\t",quote=F,col.names=T,row.names = T)
write.table(expr2_norm,file="GSE57691_normalization.txt",sep="\t",quote=F,col.names=T,row.names = T)
