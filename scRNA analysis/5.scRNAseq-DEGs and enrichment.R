rm(list=ls())
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GOplot)
library(Seurat)
library(dplyr)
load("scRNA_BCs_subclusters_annotation.RData")
Idents(scRNA)=scRNA@meta.data$celltype
deg.cluster1<- FindMarkers(scRNA,ident.1 = 'PC',ident.2 = 'Bmem',verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
library(ggplot2)
library(ggrepel)  # 用于更智能地标记文本
deg.cluster1$symbol <- rownames(deg.cluster1)

allDiff<-deg.cluster1
allDiff$gene <- rownames(deg.cluster1)
logFoldChange=0
adjustP=0.05
allDiff$type<- ifelse(allDiff$avg_log2FC> logFoldChange & allDiff$p_val_adj < adjustP, "Up",
                      ifelse(allDiff$avg_log2FC < -logFoldChange & allDiff$p_val_adj < adjustP, "Down", "not-sig")
)
deg.cluster1$Gene <- rownames(deg.cluster1)
gene_all=subset(deg.cluster1, p_val_adj<0.05)
gene_up=subset(deg.cluster1, p_val_adj<0.05&avg_log2FC>0)
gene_down=subset(deg.cluster1, p_val_adj<0.05&avg_log2FC<0)
gene_all_gene = rownames(gene_all)  
gene_up_gene = rownames(gene_up) 
gene_down_gene = rownames(gene_down) 
gene_all =gene_all%>%rownames_to_column('Gene')%>% filter(p_val_adj < 0.05)
gene_up =gene_up%>%rownames_to_column('Gene')%>% filter(p_val_adj<0.05&avg_log2FC>0)
gene_down =gene_down%>%rownames_to_column('Gene')%>% filter(p_val_adj<0.05&avg_log2FC<0)
OrgDb = "org.hsa.eg.db"
gene_all_convert<- bitr(gene_all$Gene, fromType = "SYMBOL",
                        toType = c( "ENTREZID"),
                        OrgDb =OrgDb)
gene_all=gene_all%>%inner_join(gene_all_convert,by=c("Gene"="SYMBOL"))
gene_all = gene_all_convert$ENTREZID


## GO analysis of all differentially expressed genes
kk=enrichGO(gene=gene_all, OrgDb=OrgDb, pvalueCutoff=0.05, qvalueCutoff=0.05, ont="all", readable=T)
GO=as.data.frame(kk)
pvalueFilter=0.05      
qvalueFilter=0.05     
colorSel="qvalue"
if(qvalueFilter>0.05){colorSel="pvalue"}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
go.all_all<- kk
go.all_all=DOSE::setReadable(go.all_all,OrgDb=OrgDb,keyType='ENTREZID')
GO=go.all_all@result
p1<-dotplot(go.all_all,showCategory =10,split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale='free')+scale_color_continuous(low='#bc3c29', high='#dadada')
pdf(file = 'GO_all_1.pdf',width=10,height=8)
print(p1)
dev.off()
library(GOplot)
library(circlize)
library(ComplexHeatmap)
data=GO[order(GO$p.adjust),]
datasig=data[data$p.adjust<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
pdf("GO.circlize_all.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
main.legend = Legend(
  labels = c("Biological Process", "Molecular Function","Cellular Component"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()
##KEGG analysis of all differentially expressed genes
kk.all <- enrichKEGG(gene = gene_all,
                     organism ="hsa",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

kk.all2=DOSE::setReadable(kk.all,OrgDb=OrgDb,keyType='ENTREZID')
kk.all2=kk.all2@result
p1<-dotplot(kk.all,showCategory = 15,label_format=50)+scale_color_continuous(values<-c(low='purple', high='green'), guide = guide_legend(reverse = TRUE))
pdf(file = 'KEGG_all_1.pdf',width=8,height=5)
print(p1)
dev.off()
#GO and KEGG enrichment analyses of significantly upregulated and downregulated genes were performed using the same approach as for all DEGs.
