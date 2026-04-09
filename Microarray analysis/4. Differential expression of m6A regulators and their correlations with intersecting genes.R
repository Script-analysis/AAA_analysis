setwd("C:\\AAA_analysis\\Control_AAA\\step5.Venn")
library(reshape2)
library(ggplot2)
library(ggpubr)
file<-"expr_combat.txt"
clin<-"clin_all.txt"
data<-read.table(file, header=T, sep="\t", check.names=F,row.names = 1)
clin_info<-read.table(clin, header=T, sep="\t", check.names=F,row.names = 1)
all(colnames(data) == clin_info$Sample)
clin_info <- clin_info[match(colnames(data), clin_info$Sample), , drop = FALSE]
Type<-clin_info$Group
genes_to_plot <- c(
  "ALKBH5", "FTO", "HNRNPA2B1", "HNRNPC", "IGF2BP1", "IGF2BP2", "IGF2BP3",
  "FMR1", "METTL14", "METTL3", "RBM15", "RBM15B", "CBLL1", "WTAP",
  "YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", "ZC3H13", "ELAVL1"
)
genes_use <- intersect(genes_to_plot, rownames(data))
genes_use
setdiff(genes_to_plot, rownames(data))
plot_mat <- data[genes_use, , drop = FALSE]
head(plot_mat)
plot_mat_df <- as.data.frame(plot_mat)
plot_mat_df$gene <- rownames(plot_mat_df)
library(reshape2)
plot_df <- melt(plot_mat_df, id.vars = "gene")
head(plot_df)
colnames(plot_df) <- c("Gene", "Sample", "Expression")
names(Type) <- colnames(data)
plot_df$Type <- Type[plot_df$Sample]
plot_df$Type <- factor(plot_df$Type, levels = unique(Type))
plot_df$Type<-factor(plot_df$Type,
                     levels = c("AAA",
                                "Control"))
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
my_cols <- c("Control" = "#e5d8bd", "AAA" ="#b3cde3")
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
ggsave("m6A_regulator_boxplot_with_P.pdf", plot = p, width = 10, height = 5)
setwd("C:\\AAA_analysis\\Control_AAA\\step6.Correlation")
genes_m6A_sig <- p_df %>%
  filter(!is.na(p_value) & p_value < 0.05) %>%
  pull(Gene) %>%
  as.character()
genes_m6A_sig
length(genes_m6A_sig)
file<-"intersect_DEG_DMG.txt"
intersect_genes<-read.table(file, header=T, sep="\t", check.names=F)
intersect_genes<-intersect_genes$ID
genes_target <- intersect_genes
genes_m6A_sig<- intersect(genes_m6A_sig, rownames(data))
genes_target <- intersect(genes_target, rownames(data))
length(genes_m6A_sig)
length(genes_target)
expr_m6A <- t(data[genes_m6A_sig, , drop = FALSE])
expr_target <- t(data[genes_target, , drop = FALSE])
cor_mat <- cor(expr_m6A, expr_target, method = "pearson")
cor_p <- matrix(NA, nrow = ncol(expr_m6A), ncol = ncol(expr_target))
rownames(cor_p) <- colnames(expr_m6A)
colnames(cor_p) <- colnames(expr_target)

for(i in 1:ncol(expr_m6A)){
  for(j in 1:ncol(expr_target)){
    test <- cor.test(expr_m6A[, i], expr_target[, j])
    cor_p[i, j] <- test$p.value
  }
}
sig_mat <- ifelse(cor_p < 0.05, "*", "")
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(dplyr)
p_star <- matrix(
  "",
  nrow = nrow(cor_p),
  ncol = ncol(cor_p),
  dimnames = dimnames(cor_p)
)
p_star[cor_p < 0.05]  <- "*"
p_star[cor_p < 0.01]  <- "**"
p_star[cor_p < 0.001] <- "***"
cor_df <- melt(cor_mat)
p_df   <- melt(cor_p)
colnames(cor_df) <- c("m6A_gene", "Target_gene", "Correlation")
colnames(p_df)   <- c("m6A_gene", "Target_gene", "P_value")
result_df <- merge(cor_df, p_df, by = c("m6A_gene", "Target_gene"))
result_df <- result_df %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH"),
    Significance = case_when(
      is.na(P_value) ~ "",
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
pdf("cor_pearson_significant_only.pdf", height = 5, width = 10)
ht <- Heatmap(
  cor_mat,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "grey85", lwd = 0.5),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(!is.na(cor_p[i, j]) && cor_p[i, j] < 0.05){
      grid.text(
        sprintf("%.2f", cor_mat[i, j]),
        x = x,
        y = y + h * 0.18,
        gp = gpar(fontsize = 8, col = "black")
      )
      grid.text(
        p_star[i, j],
        x = x,
        y = y - h * 0.18,
        gp = gpar(fontsize = 10, col = "black", fontface = "bold")
      )
    }
  }
)
draw(ht)
dev.off()
