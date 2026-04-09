setwd("C:\\AAA_analysis\\Control_AAA\\step4.Venn")
library(ggVennDiagram)
library(readr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(openxlsx)
GSE7084   <- read.delim("DEGs_GSE7084_sig.txt", row.names = 1, check.names = FALSE)
GSE57691  <- read.delim("DEGs_GSE57691_sig.txt", row.names = 1, check.names = FALSE)
Merge <- read.delim("DEGs_Merge_sig.txt", row.names = 1, check.names = FALSE)
DMG <- read.xlsx("diff_peak_genes.xlsx")
library(VennDiagram)
library(grid)
genes1 <- unique(na.omit(rownames(GSE7084)))
genes2 <- unique(na.omit(rownames(GSE57691)))
genes3 <- unique(na.omit(DMG$Geneid))
genes4 <- unique(na.omit(rownames(Merge)))
pdf("VennDiagram.pdf", width = 8, height = 8)
venn.plot <- venn.diagram(
  x = list(
    GSE7084   = genes1,
    GSE57691  = genes2,
    DMG = genes3,
    Merge = genes4
  ),
  filename = NULL,
  fill = c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488CC"),
  col = "white",
  lwd = 2,
  cex = 1.2,
  cat.cex = 1.2,
  cat.fontface = "bold",
  fontface = "bold",
  cat.col = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"),
  margin = 0.08
)

grid.newpage()
grid.draw(venn.plot)
dev.off()
common_genes <- Reduce(intersect, list(genes1, genes2, genes3, genes4))
write.table(
  common_genes,
  file = "intersect_DEG_DMG.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
