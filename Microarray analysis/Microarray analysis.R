####################################
library(dplyr)
library(tidyverse)
library(data.table)
library(limma) 
library(ggplot2)
library(ggrepel)
#################################DEGs
rt=read.table("data.txt",sep="\t",header=T,check.names=F,row.names = 1)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
boxplot(data)
data<-normalizeBetweenArrays(data)
group<-"group.txt"
clin_info<-read.table(group,sep="\t",header=T,check.names=F)
Group<-clin_info$Group
clin_info <- clin_info[match(colnames(data), rownames(clin_info)), , drop = FALSE]
all(colnames(data) == rownames(clin_info))
design <- model.matrix(~0+factor(clin_info$Group))
colnames(design) <- levels(factor(clin_info$Group))
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(AAA-Ctrl,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
logFoldChange=1
adjustP=0.05
allDiff$type<- ifelse(allDiff$logFC> logFoldChange & allDiff$adj.P.Val < adjustP, "Up",
                      ifelse(allDiff$logFC < -logFoldChange & allDiff$adj.P.Val < adjustP, "Down", "not-sig"))
write.table(allDiff, "DEG_all.txt", sep = "\t", quote = FALSE, row.names = TRUE)
gene_all=subset(allDiff, adj.P.Val<0.05 & (logFC>logFoldChange|logFC< -logFoldChange))
gene_up=subset(allDiff, adj.P.Val<0.05&logFC>logFoldChange)
gene_down=subset(allDiff, adj.P.Val<0.05&logFC< -logFoldChange)
write.table(gene_all,file="Sig_all.txt",sep="\t",quote=F,row.names=T)
write.table(gene_up ,file="Sig_up.txt",sep="\t",quote=F,row.names=T)
write.table(gene_down,file="Sig_down.txt",sep="\t",quote=F,row.names=T)

####################################Visualization
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
expFile="data_foucus_genes.txt"    
data=read.table(expFile, header=T, sep="\t", check.names=F,row.names = 1)
group<-"group.txt"
clin_info<-read.table(group,sep="\t",header=T,check.names=F)
clin_info<- clin_info [colnames(data), ]
clin_info <- clin_info[match(colnames(data), rownames(clin_info)), , drop = FALSE]
Type<-clin_info$Group
expr=as.data.frame(t(data))
exp=cbind(expr, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")
data$Gene<- factor(data$Gene, levels = unique(data$Gene)) 
data$Type<-factor(data$Type, levels = c("Ctrl","AAA"))
#boxplot
p=ggplot(data, aes(x=Gene, y=Expression))+
  geom_boxplot(aes(color=Type) ,width=1)+
  theme_classic()+
  theme(panel.grid = element_blank(),axis.line.y = element_line(),axis.ticks.y = element_line())+
  facet_wrap(~Gene,scales='free')
p1<-p+stat_compare_means(aes(group=Type),
                         method="wilcox.test",
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),colour='black',size = 0.1,
                         label = "p.signif",label.x=1.5)

pdf(file="boxplot.pdf", width=7.5, height=5)
print(p1)
dev.off()
#half_violin
colors <- c("#EA8379", "#83c5db")  
group_labels <- c("AAA_Female","AAA_Male")
create_violin_plot <- function(data, colors, group_col, group_labels) {
data_long <- data %>%
    tidyr::pivot_longer(
      cols = -all_of(group_col),
      names_to = "Genes",
      values_to = "Values"
    ) %>%
    mutate(!!sym(group_col) := factor(!!sym(group_col), levels = group_labels))
  y_min <- min(data_long$Values, na.rm = TRUE)
  y_max <- max(data_long$Values, na.rm = TRUE)
  y_range <- y_max - y_min
  y_upper <- y_max + y_range * 0.1
  p <- ggplot(data_long, aes(x = Genes, y = Values)) +
    geom_half_violin(
      data = . %>% filter(!!sym(group_col) == group_labels[1]),
      aes(fill = !!sym(group_col)), 
      side = "l", 
      alpha = 0.7,
      trim = FALSE,
      scale = "width"
    ) +
    geom_half_violin(
      data = . %>% filter(!!sym(group_col) == group_labels[2]),
      aes(fill = !!sym(group_col)), 
      side = "r", 
      alpha = 0.7,
      trim = FALSE,
      scale = "width"
    ) +
    stat_summary(
      aes(group = !!sym(group_col)),#####
      fun = mean, 
      geom = "point", 
      shape = 21,
      size = 3,
      color = "black",
      position = position_dodge(0.8)
    ) +
    stat_summary(
      aes(group = !!sym(group_col)),
      fun.min = function(z) {quantile(z, 0.25)},
      fun.max = function(z) {quantile(z, 0.75)},
      geom = "errorbar", 
      width = 0.15,
      size = 0.8,
      color = "black",
      position = position_dodge(0.8)
    ) +
    stat_compare_means(
      aes(group = !!sym(group_col)),
      method = "wilcox.test",
      label = "p.signif",group.by="Genes",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),
      label.y = y_upper,
      size = 2,
      vjust = 0.5,
      hide.ns = FALSE
    ) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "Expression Value") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 45, 
        hjust = 1,
        size = 12,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      axis.title.y = element_text(
        size = 14,
        face = "bold",
        margin = margin(r = 15)
      ),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  return(p)
}
violin_plot <- create_violin_plot(
  data = data, 
  colors = colors, 
  group_col = "group", 
  group_labels = group_labels
)
print(violin_plot)
ggsave("Violin_plot.pdf",width=10,height=4)

####################################CIBERSORT
#The code was obtained from https://github.com/Moonerss/CIBERSORT

####################################Correlation analysis
data <- read.table("data.txt", header = TRUE, sep = "\t", quote = "\"",
                   fill = TRUE, comment.char = "", check.names = FALSE, row.names = 1,
                   stringsAsFactors = FALSE)
data1<- read.table("CIBERSORT_result.txt", header = TRUE, sep = "\t", quote = "\"",
                   fill = TRUE, comment.char = "", check.names = FALSE, row.names = 1,
                   stringsAsFactors = FALSE)
exp<-rbind(data,result1)
suppressPackageStartupMessages({
  library(Hmisc)
  library(corrplot)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
group<-"group.txt"
clin_info<-read.table(group,sep="\t",header=T,check.names=F,row.names = 1)
keep_samples <- rownames(clin_info)[clin_info$Group %in% c("AAA_Female")]
clin_info <- clin_info[keep_samples, , drop=FALSE]
data <- exp[, keep_samples, drop=FALSE]
method <- "pearson"
method <- "spearman"
stopifnot("FKBP11" %in% rownames(data))
td <- t(data)
gene_order_data <- intersect(rownames(data), colnames(td))
gene_order <- c("FKBP11", setdiff(gene_order_data, "FKBP11"))
n_samp <- nrow(td)
cor_test_matrix <- function(mat, method = method){
  g <- colnames(mat); n <- length(g)
  rmat <- matrix(NA_real_, n, n, dimnames = list(g, g))
  pmat <- matrix(NA_real_, n, n, dimnames = list(g, g))
  for(i in seq_len(n)){
    for(j in i:n){
      x <- mat[, i]; y <- mat[, j]
      ok <- complete.cases(x, y)
      if(sum(ok) >= 3){
        if (method == "spearman"){
          ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = TRUE))
        } else {
          ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
        }
        rmat[i, j] <- rmat[j, i] <- unname(ct$estimate)
        pmat[i, j] <- pmat[j, i] <- ct$p.value
      } else {
        r <- suppressWarnings(cor(x, y, method = method, use = "pairwise.complete.obs"))
        rmat[i, j] <- rmat[j, i] <- r
        pmat[i, j] <- pmat[j, i] <- NA_real_
      }
    }
  }
  list(r = rmat, p = pmat)
}
if (n_samp >= 5){
  rc <- Hmisc::rcorr(as.matrix(td), type = method)
  rmat <- rc$r; pmat <- rc$P
} else {
  res <- cor_test_matrix(td, method = method)
  rmat <- res$r; pmat <- res$p
}
show_sig <- TRUE
r_fk <- rmat[gene_order, "FKBP11"]
p_fk <- pmat[gene_order, "FKBP11"]
res_one <- data.frame(
  gene = gene_order,
  r = as.numeric(r_fk),
  p = as.numeric(p_fk),
  stringsAsFactors = FALSE
)
write.csv(res_one, sprintf("FKBP11_vs_all_genes_cor_byDataOrder_%s.csv", method), row.names = FALSE)
write.csv(data.frame(gene = gene_order), "FKBP11_gene_order_from_data.csv", row.names = FALSE)
plot_df <- res_one[res_one$gene != "FKBP11", ]
plot_df$gene <- factor(plot_df$gene, levels = plot_df$gene)
label_sig_fun <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}
pdf(sprintf("FKBP11_corr_oneRow_heatmap_byDataOrder_%s.pdf", method), width = 3, height = 5)
ggplot(plot_df, aes(x = gene, y = "FKBP11", fill = r)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), low = "blue", mid = "white", high = "red") +
  geom_text(aes(label = label_sig_fun(p)), size = 3) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = sprintf("%s r", tools::toTitleCase(method)),
       title = sprintf("FKBP11 vs genes (order = rows of data, n=%d)", n_samp)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))
dev.off()
