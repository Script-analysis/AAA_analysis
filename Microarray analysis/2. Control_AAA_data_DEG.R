#DEGs
setwd("C:\\AAA_analysis\\Control_AAA\\step3.DEGs")
file1="GSE7084_normalization.txt"
file2="clin_GSE7084.txt"
data<- read.table(file1, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<- read.table(file2, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
design <- model.matrix(~0+factor(clin_info$Group))
colnames(design) <- levels(factor(clin_info$Group))
rownames(design) <- colnames(data)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(AAA-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
logFCfilter=1        
adj.P.Val.Filter=0.05   
allDiff$type <- ifelse(allDiff$logFC > logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Up",
                       ifelse(allDiff$logFC < -logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Down", "not-sig"))
write.table(allDiff, file="DEGs_GSE7084.txt", sep="\t", quote=F, row.names=TRUE, col.names=TRUE)
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
write.table(diffSig, file="DEGs_GSE7084_sig.txt", sep="\t", quote=F, col.names=T)
diffUp=diffSig[diffSig$logFC>logFCfilter,]
diffDown=diffSig[diffSig$logFC< -logFCfilter,]
write.table(diffUp, file="DEGs_GSE7084_sig_up.txt", sep="\t", quote=F, col.names=T)
write.table(diffDown, file="DEGs_GSE7084_sig_down.txt", sep="\t", quote=F, col.names=T)

library(ggplot2)
library(ggrepel)
volcano_df <- allDiff
volcano_df$gene <- rownames(volcano_df)
volcano_df$negLogP <- -log10(volcano_df$adj.P.Val)
volcano_colors <- c("Up" = "#E64B35FF", "Down" = "#4DBBD5FF","not-sig" = "grey80"
)
library(ggrepel)
top_up <- volcano_df %>% filter(type=="Up") %>% arrange(adj.P.Val) %>% head(5)
top_down <- volcano_df %>% filter(type=="Down") %>% arrange(adj.P.Val) %>% head(5)
top_genes <- rbind(top_up, top_down)
plot <- ggplot(volcano_df, aes(logFC, negLogP)) +
  geom_point(data = subset(volcano_df, type=="not-sig"),
             color = "grey85", size = 1) +
  geom_point(data = subset(volcano_df, type=="Up"),
             color = "#E64B35FF", size = 1.5) +
  geom_point(data = subset(volcano_df, type=="Down"),
             color = "#4DBBD5FF", size = 1.5) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey50") +
  geom_label_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    color = "black",             
    fill = "white",               
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    label.size = 0.2,            
    max.overlaps = Inf
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "AAA vs Control",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
ggsave("Volcano Plot_top10_GSE7084.pdf",plot=plot,width=5,height=4)
file1="GSE57691_normalization.txt"
file2="clin_GSE57691.txt"
data<- read.table(file1, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<- read.table(file2, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
design <- model.matrix(~0+factor(clin_info$Group))
colnames(design) <- levels(factor(clin_info$Group))
rownames(design) <- colnames(data)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(AAA-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
logFCfilter=1        
adj.P.Val.Filter=0.05   
allDiff$type <- ifelse(allDiff$logFC > logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Up",
                       ifelse(allDiff$logFC < -logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Down", "not-sig"))
write.table(allDiff, file="DEGs_GSE57691.txt", sep="\t", quote=F, row.names=TRUE, col.names=TRUE)
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
write.table(diffSig, file="DEGs_GSE57691_sig.txt", sep="\t", quote=F, col.names=T)
diffUp=diffSig[diffSig$logFC>logFCfilter,]
diffDown=diffSig[diffSig$logFC< -logFCfilter,]
write.table(diffUp, file="DEGs_GSE57691_sig_up.txt", sep="\t", quote=F, col.names=T)
write.table(diffDown, file="DEGs_GSE57691_sig_down.txt", sep="\t", quote=F, col.names=T)
library(ggplot2)
library(ggrepel)
volcano_df <- allDiff
volcano_df$gene <- rownames(volcano_df)
volcano_df$negLogP <- -log10(volcano_df$adj.P.Val)
volcano_colors <- c(
  "Up" = "#E64B35FF",      
  "Down" = "#4DBBD5FF",    
  "not-sig" = "grey80"
)
library(ggrepel)
top_up <- volcano_df %>% filter(type=="Up") %>% arrange(adj.P.Val) %>% head(5)
top_down <- volcano_df %>% filter(type=="Down") %>% arrange(adj.P.Val) %>% head(5)
top_genes <- rbind(top_up, top_down)
plot <- ggplot(volcano_df, aes(logFC, negLogP)) +
  geom_point(data = subset(volcano_df, type=="not-sig"),
             color = "grey85", size = 1) +
  geom_point(data = subset(volcano_df, type=="Up"),
             color = "#E64B35FF", size = 1.5) +
  geom_point(data = subset(volcano_df, type=="Down"),
             color = "#4DBBD5FF", size = 1.5) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey50") +
  geom_label_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    color = "black",             
    fill = "white",               
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    label.size = 0.2,            
    max.overlaps = Inf
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "AAA vs Control",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
ggsave("Volcano Plot_top10_GSE7084.pdf",plot=plot,width=5,height=4)

##batch effect correction
file1<-"GSE57691_normalization.txt"
file2<-"GSE7084_normalization.txt"
expr1_norm<-read.table(file1, header=T, sep="\t", check.names=F,row.names = 1)
expr2_norm<-read.table(file2, header=T, sep="\t", check.names=F,row.names = 1)
common_genes <- intersect(rownames(expr1_norm), rownames(expr2_norm))
length(common_genes)
expr1_common <- expr1_norm[common_genes, ]
expr2_common <- expr2_norm[common_genes, ]
clin1="clin_GSE98278.txt"
clin2="clin_GSE7084.txt"
meta1 <-read.table(clin1, header=T, sep="\t", check.names=F,row.names = 1)
meta2 <-read.table(clin2, header=T, sep="\t", check.names=F,row.names = 1)
expr_merge <- cbind(expr1_common, expr2_common)
dim(expr_merge)
write.table(expr_merge,file="expr_merge_Raw.txt",sep="\t",quote=F,col.names=T,row.names = T)

meta_all <- rbind(meta1, meta2)
meta_all <- meta_all[match(colnames(expr_merge), meta_all$Sample), ]
write.table(meta_all,file="clin_all.txt",sep="\t",quote=F,col.names=T,row.names = T)
meta_all$Cohort <- factor(meta_all$Cohort)
meta_all$Group <- factor(meta_all$Group, levels = c("Control", "AAA"))
data<-expr_merge
clin_info<-meta_all
boxplot(data)
library(tinyarray)
draw_pca(exp = data, group_list = factor(clin_info$Group))
ggsave("pre_PCA_Group.pdf",width=5,height=4)
draw_pca(exp = data, group_list = factor(clin_info$Cohort))
ggsave("pre_PCA_Cohort.pdf",width=5,height=4)
library(sva)
expr_combat <- ComBat(dat = data, batch = clin_info$Cohort)
expr_combat[1:4,1:4]
expr_combat <- as.data.frame(expr_combat)
draw_pca(exp = expr_combat, group_list = factor(clin_info$Group))
ggsave("post_expr_combat_PCA_Group.pdf",width=5,height=4)
write.table(expr_combat, file="expr_combat.txt", sep="\t", quote=F, col.names=T)
file1="expr_combat.txt"
file2="clin_all.txt"
data<- read.table(file1, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<- read.table(file2, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
design <- model.matrix(~0+factor(clin_info$Group))
colnames(design) <- levels(factor(clin_info$Group))
rownames(design) <- colnames(data)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(AAA-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
logFCfilter=1        
adj.P.Val.Filter=0.05   
allDiff$type <- ifelse(allDiff$logFC > logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Up",
                       ifelse(allDiff$logFC < -logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, "Down", "not-sig"))
write.table(allDiff, file="DEGs_Merge.txt", sep="\t", quote=F, row.names=TRUE, col.names=TRUE)
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
write.table(diffSig, file="DEGs_Merge_sig.txt", sep="\t", quote=F, col.names=T)
diffUp=diffSig[diffSig$logFC>logFCfilter,]
diffDown=diffSig[diffSig$logFC< -logFCfilter,]
write.table(diffUp, file="DEGs_Merge_sig_up.txt", sep="\t", quote=F, col.names=T)
write.table(diffDown, file="DEGs_Merge_sig_down.txt", sep="\t", quote=F, col.names=T)
library(ggplot2)
library(ggrepel)
volcano_df <- allDiff
volcano_df$gene <- rownames(volcano_df)
volcano_df$negLogP <- -log10(volcano_df$adj.P.Val)
volcano_colors <- c(
  "Up" = "#E64B35FF",      
  "Down" = "#4DBBD5FF",    
  "not-sig" = "grey80"
)
library(ggrepel)
top_up <- volcano_df %>% filter(type=="Up") %>% arrange(adj.P.Val) %>% head(5)
top_down <- volcano_df %>% filter(type=="Down") %>% arrange(adj.P.Val) %>% head(5)
top_genes <- rbind(top_up, top_down)
plot <- ggplot(volcano_df, aes(logFC, negLogP)) +
  geom_point(data = subset(volcano_df, type=="not-sig"),
             color = "grey85", size = 1) +
  geom_point(data = subset(volcano_df, type=="Up"),
             color = "#E64B35FF", size = 1.5) +
  geom_point(data = subset(volcano_df, type=="Down"),
             color = "#4DBBD5FF", size = 1.5) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey50") +
  geom_label_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    color = "black",             
    fill = "white",               
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    label.size = 0.2,            
    max.overlaps = Inf
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "AAA vs Control",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
ggsave("Volcano Plot_top10_Merge.pdf",plot=plot,width=5,height=4)
