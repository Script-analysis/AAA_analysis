setwd("C:\\AAA_analysis\\sAAA_lAAA\\step2.DMR and DEGs in sAAA and lAAA")
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(sva)
genes="DMR_DEGs_18genes.txt"
genes<-read.table(genes, header=T, sep="\t", check.names=F)
file<-"GSE98278_normalizationt.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
genes <- genes$ID
genes_use <- intersect(genes, rownames(data))
data_sub <- data[genes_use, ,drop = FALSE]
write.table(data_sub,file="GSE98278_18genes.txt",sep="\t",quote=F,col.names=T,row.names = T)
file<-"GSE57691_normalizationt.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
genes <- genes$ID
genes_use <- intersect(genes, rownames(data))
data_sub <- data[genes_use, ,drop = FALSE]
write.table(data_sub,file="GSE57691_18genes.txt",sep="\t",quote=F,col.names=T,row.names = T)

library(reshape2)
library(ggplot2)
library(ggpubr)
file<-"GSE98278_18genes.txt"
clin<-"clin_GSE98278.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<-read.table(clin, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
Type<-clin_info$group
plot_mat_df <- as.data.frame(data)
plot_mat_df$gene <- rownames(plot_mat_df)
library(reshape2)
plot_df <- melt(plot_mat_df, id.vars = "gene")
head(plot_df)
colnames(plot_df) <- c("Gene", "Sample", "Expression")
names(Type) <- colnames(data)
plot_df$Type <- Type[plot_df$Sample]
plot_df$Type <- factor(plot_df$Type, levels = unique(Type))
plot_df$Type<-factor(plot_df$Type,
                     levels = c("sAAA",
                                "lAAA"))
levels(plot_df$Type)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gghalves)
p_df <- plot_df %>%
  group_by(Gene) %>%
  dplyr::summarise(
    p_value = tryCatch(
      wilcox.test(Expression ~ Type)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(Expression, na.rm = TRUE) + 0.35,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_label = ifelse(is.na(p_value), "P=NA", paste0("P=", signif(p_value, 2)))
  )
my_cols <- c("sAAA" = "#decbe4", "lAAA" ="#ccebc5")
pd <- position_dodge(width = 0.65)

p <- ggplot(plot_df, aes(x = Gene, y = Expression, fill = Type)) +
  geom_boxplot(
    position = pd,
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.9,
    color = "black",
    linewidth = 0.4
  ) +
  geom_point(
    aes(fill = Type),
    shape = 21,
    color = "grey30",   
    stroke = 0.4,      
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.65
    ),
    size = 1.8,
    alpha = 0.8
  )+
  geom_text(
    data = p_df,
    aes(x = Gene, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "italic"
  ) +
  scale_fill_manual(values = my_cols) +
  scale_color_manual(values = my_cols) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
p
ggsave("boxplot_with_P_GSE98278.pdf", plot = p, width = 10, height = 5)
library(reshape2)
library(ggplot2)
library(ggpubr)
file<-"GSE57691_18genes.txt"
clin<-"clin_GSE57691.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<-read.table(clin, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
Type<-clin_info$group
plot_mat_df <- as.data.frame(data)
plot_mat_df$gene <- rownames(plot_mat_df)
library(reshape2)
plot_df <- melt(plot_mat_df, id.vars = "gene")
head(plot_df)
colnames(plot_df) <- c("Gene", "Sample", "Expression")
names(Type) <- colnames(data)
plot_df$Type <- Type[plot_df$Sample]
plot_df$Type <- factor(plot_df$Type, levels = unique(Type))
plot_df$Type<-factor(plot_df$Type,
                     levels = c("sAAA",
                                "lAAA"))
levels(plot_df$Type)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gghalves)
p_df <- plot_df %>%
  group_by(Gene) %>%
  dplyr::summarise(
    p_value = tryCatch(
      wilcox.test(Expression ~ Type)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(Expression, na.rm = TRUE) + 0.35,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_label = ifelse(is.na(p_value), "P=NA", paste0("P=", signif(p_value, 2)))
  )
my_cols <- c("sAAA" = "#decbe4", "lAAA" ="#ccebc5")
pd <- position_dodge(width = 0.65)

p <- ggplot(plot_df, aes(x = Gene, y = Expression, fill = Type)) +
  geom_boxplot(
    position = pd,
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.9,
    color = "black",
    linewidth = 0.4
  ) +
  geom_point(
    aes(fill = Type),
    shape = 21,
    color = "grey30",   
    stroke = 0.4,      
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.65
    ),
    size = 1.8,
    alpha = 0.8
  )+
  geom_text(
    data = p_df,
    aes(x = Gene, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "italic"
  ) +
  scale_fill_manual(values = my_cols) +
  scale_color_manual(values = my_cols) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
p
ggsave("boxplot_with_P_GSE57691.pdf", plot = p, width = 10, height = 5)