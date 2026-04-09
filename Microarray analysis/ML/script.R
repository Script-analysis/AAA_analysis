set.seed(777)
work.path <- "C:\\AAA_analysis\\Control_AAA\\step8.ML"
setwd(work.path) 
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")
for (p in c(code.path, data.path, res.path, fig.path)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
suppressPackageStartupMessages({
  library(openxlsx)
  library(plyr)
  library(randomForestSRC)
  library(glmnet)
  library(plsRglm)
  library(gbm)
  library(caret)
  library(mboost)
  library(e1071)
  library(BART)
  library(MASS)
  library(xgboost)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pROC)
})
source(file.path(code.path, "ML.R"))
FinalModel <- c("panML", "multiLogistic")[1]
Train_expr <- read.table(file.path(data.path, "GSE57691_18genes.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
Train_expr <- as.data.frame(Train_expr)
Train_class <- read.table(file.path(data.path, "clin_GSE57691.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
colnames(Train_class)[colnames(Train_class) == "Group"] <- "outcome"
Train_class$outcome <- ifelse(Train_class$outcome == "Control", 0,
                              ifelse(Train_class$outcome == "AAA", 1, NA))
Train_class$outcome <- as.numeric(Train_class$outcome)
Test_expr <- read.table(file.path(data.path, "GSE7084_combat_18genes.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
Test_expr <- as.data.frame(Test_expr)
Test_class <- read.table(file.path(data.path, "clin_GSE7084.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
colnames(Test_class)[colnames(Test_class) == "Group"] <- "outcome"
Test_class$outcome <- ifelse(Test_class$outcome == "Control", 0,
                             ifelse(Test_class$outcome == "AAA", 1, NA))
Test_class$outcome <- as.numeric(Test_class$outcome)
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = FALSE]
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = FALSE]
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,])
Test_expr <- t(Test_expr[comgene,])
Train_set <- as.matrix(Train_expr)
Test_set  <- as.matrix(Test_expr)
common_features <- intersect(colnames(Train_set), colnames(Test_set))
Train_set <- Train_set[, common_features, drop = FALSE]
Test_set  <- Test_set[, common_features, drop = FALSE]
train_mean <- apply(Train_set, 2, mean, na.rm = TRUE)
train_sd   <- apply(Train_set, 2, sd, na.rm = TRUE)
train_sd[is.na(train_sd) | train_sd == 0] <- 1
Train_set <- scale(Train_set, center = train_mean, scale = train_sd)
Test_set  <- scale(Test_set,  center = train_mean, scale = train_sd)
Train_set <- as.data.frame(Train_set)
Test_set  <- as.data.frame(Test_set)
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
classVar = "outcome"
min.selected.var = 5
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(777)
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,
                                 Train_set = Train_set,
                                 Train_label = Train_class,
                                 mode = "Variable",
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

model <- list()
set.seed(777)
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable = preTrain.var[[method[1]]]
  if (is.null(Variable) || length(Variable) < 1) next
  Train_set = Train_set_bk[, Variable, drop=FALSE]
  Train_label = Train_class
  fit <- RunML(method = method[2],
               Train_set = Train_set,
               Train_label = Train_label,
               mode = "Model",
               classVar = classVar)
  model[[method_name]] <- fit
}
Train_set = Train_set_bk; rm(Train_set_bk)
saveRDS(model, file.path(res.path, "model.rds")) 
if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){
    vars <- ExtractVar(fit); if (!length(vars)) vars <- fit$subFeature
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, vars, drop=FALSE]))
    tmp$subFeature <- vars
    pr <- as.numeric(predict(tmp, type="response"))
    tmp$threshold <- {
      y <- factor(Train_class[[classVar]])
      rocobj <- suppressMessages(roc(y, pr, quiet=TRUE, direction="<"))
      as.numeric(coords(rocobj, "best", ret="threshold", best.method="youden"))
    }
    return(tmp)
  })
  saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds"))
}

use_model <- readRDS(file.path(res.path, "model.rds"))
methodsValid <- names(use_model)

RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = use_model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = use_model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(Class_mat, file.path(res.path, "Class_mat.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

fea_list <- lapply(use_model, ExtractVar)
fea_df <- lapply(use_model, function(fit){ data.frame(ExtractVar(fit)) })
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

AUC_list <- list()
AUC_CI_low <- list(); AUC_CI_up <- list()
for (method in methodsValid){
  auc_vec <- RunEval(fit = use_model[[method]],
                     Test_set = Test_set,
                     Test_label = Test_class,
                     Train_set = Train_set,
                     Train_label = Train_class,
                     Train_name = "Train",
                     cohortVar = "Cohort",
                     classVar = classVar)
  AUC_list[[method]] <- auc_vec
  ci_list <- attr(auc_vec, "ci")
  AUC_CI_low[[method]] <- sapply(ci_list, function(ci) as.numeric(ci[1]))
  AUC_CI_up[[method]]  <- sapply(ci_list, function(ci) as.numeric(ci[3]))
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(do.call(rbind, AUC_CI_low), file.path(res.path, "AUC_CI_low.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(do.call(rbind, AUC_CI_up),  file.path(res.path, "AUC_CI_up.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), , drop=FALSE]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

if(ncol(AUC_mat) < 3) {  CohortCol <- c("#D55E00", "#0072B2") } else { CohortCol <- brewer.pal(n = max(3, ncol(AUC_mat)), name = "Paired") }
names(CohortCol) <- colnames(AUC_mat)

hm <- SimpleHeatmap(AUC_mat,
                    avg_AUC,
                    CohortCol, "steelblue", 
                    cellwidth = 1, cellheight = 0.5,
                    cluster_columns = FALSE, cluster_rows = FALSE)

pdf(file.path(fig.path, "AUC_1.pdf"), width = 8, height = 22)
draw(hm)
invisible(dev.off())
coef_mat <- as.matrix(model$Ridge$beta)
coef_df <- data.frame(
  Gene = rownames(coef_mat),
  Coefficient = coef_mat[,1]
)
coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
head(coef_df)
write.table(coef_df,  file.path(res.path, "coef_df.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
importance_df<-coef_df
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$Coefficient, decreasing = TRUE), ]
top_vars <- head(importance_df, 20)
library(ggplot2)
ggplot(top_vars, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_segment(aes(x = reorder(Gene, Coefficient),
                   xend = reorder(Gene, Coefficient),
                   y = 0, yend = Coefficient),
               color = "grey70", size = 1) +
  geom_point(size = 4, color = "steelblue") +
  coord_flip() +
  labs(title = "",
       x = "Variables",
       y = "Coefficient") +
  theme_bw(base_size = 14)
ggsave("Importance_lollipop.pdf",height=5,width=5)

#################################################
model <- readRDS(file.path(res.path, "model.rds"))
library(fastshap)
library(ggplot2)
library(reshape2)  
rf_fit  <- model[["Ridge"]]
stopifnot(!is.null(rf_fit))

feat <- rf_fit$subFeature 
X_train <- as.data.frame(Train_set[, feat, drop = FALSE])

# randomForest 
pred_wrapper_rfsrc <- function(object, newdata) {
  newdata <- as.data.frame(newdata)[, object$subFeature, drop = FALSE]
  pr <- predict(object, newdata)$predicted
  pos_col <- if ("1" %in% colnames(pr)) "1" else colnames(pr)[ncol(pr)]
  as.numeric(pr[, pos_col])
}
#Ridge
rf_fit <- model[["Ridge"]]
stopifnot(!is.null(rf_fit))
feat <- rf_fit$subFeature
X_train <- as.data.frame(Train_set[, feat, drop = FALSE])
pred_wrapper_glmnet <- function(object, newdata) {
  newdata <- as.matrix(newdata[, object$subFeature, drop = FALSE])
  pr <- predict(object, newx = newdata, type = "response")
  as.numeric(pr)
}
set.seed(2025)
shap_vals <- fastshap::explain(
  object = rf_fit,
  X = X_train,
  pred_wrapper = pred_wrapper_glmnet,
  nsim = 256
)
mean_abs <- sort(colMeans(abs(shap_vals)), decreasing = TRUE)
topK <- 20L
imp_df <- data.frame(
  Variable = names(mean_abs)[seq_len(min(topK, length(mean_abs)))],
  MeanAbsSHAP = as.numeric(mean_abs[seq_len(min(topK, length(mean_abs)))])
)

ggplot(imp_df, aes(x = reorder(Variable, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_segment(aes(xend = Variable, y = 0, yend = MeanAbsSHAP), color = "grey70", linewidth = 1) +
  geom_point(size = 3) +
  coord_flip() +
  labs(title = "Global SHAP Importance",
       x = "Variables", y = "Mean |SHAP|") +
  theme_bw(base_size = 14)

ggsave("SHAP Importance_1.pdf",height=5,width=5)  
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(scales)
shap_long <- melt(shap_vals, varnames = c("sample", "feature"), value.name = "phi")
colnames(shap_long)
stopifnot(all(colnames(shap_vals) %in% colnames(X_train)))
for (f in unique(shap_long$feature)) {
  shap_long$raw_value[shap_long$feature == f] <- X_train[[f]]
}
bin_by_quantile <- function(x) {
  qs <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE)
  cut(x,
      breaks = c(-Inf, qs[1], qs[2], Inf),
      labels = c("Low", "Mid", "High"),
      include.lowest = TRUE, right = TRUE)
}
shap_long <- shap_long %>%
  group_by(feature) %>%
  mutate(value_group = bin_by_quantile(raw_value)) %>%
  ungroup()
shap_long$group_num <- dplyr::recode(shap_long$value_group,
                                     "Low" = 0, "Mid" = 0.5, "High" = 1)
feature_order <- shap_long %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_abs_phi = mean(abs(phi), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_abs_phi)) %>%
  dplyr::pull(feature)

shap_long$feature <- factor(shap_long$feature, levels = rev(feature_order)) 
p_beeswarm <- ggplot(shap_long, aes(x = phi, y = feature, color = group_num)) +
  geom_quasirandom(groupOnX = FALSE, alpha = 0.8, size = 1.0) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name   = "Feature value",
                        breaks = c(0, 0.5, 1),
                        labels = c("Low", "Mid", "High")) +
  theme_bw(base_size = 14) +
  labs(
    title = "SHAP Value Distribution",
    x = "SHAP value",
    y = NULL
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
p_beeswarm
ggsave("beeswarm.pdf", plot = p_beeswarm, width = 7, height = 5)

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(tidyr)
topK <- 10                
grid_cols <- 3           
limit_x_by_quantile <- TRUE   
qx <- c(0.01, 0.99)      
pt_alpha <- 0.7
pt_size  <- 2
loess_span <- 0.8
pdf_file_grid <- "dependency_top10_grid.pdf"
pdf_file_col1 <- "dependency_top10_col.pdf"
shap_all_df <- do.call(rbind, lapply(seq_len(ncol(shap_vals)), function(j){
  data.frame(
    feature = colnames(shap_vals)[j],
    phi     = shap_vals[, j],
    value   = X_train[, colnames(shap_vals)[j]],
    stringsAsFactors = FALSE
  )
}))
shap_all_df <- tidyr::drop_na(shap_all_df, phi, value)
top_features <- shap_all_df %>%
  group_by(feature) %>%
  dplyr::summarise(mean_abs_phi = mean(abs(phi), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_abs_phi)) %>%
  slice_head(n = topK) %>%
  pull(feature)

print(top_features)
plots <- lapply(top_features, function(fea) {
  df_plot <- shap_all_df %>% filter(feature == fea)
  if (limit_x_by_quantile && nrow(df_plot) > 5) {
    xlim_q <- quantile(df_plot$value, probs = qx, na.rm = TRUE, names = FALSE)
  } else {
    xlim_q <- range(df_plot$value, na.rm = TRUE)
  }
  
  ggplot(df_plot, aes(x = value, y = phi, color = value)) +
    geom_point(alpha = pt_alpha, size = pt_size) +
    geom_smooth(method = "loess", se = FALSE, color = "black",
                linewidth = 1, span = loess_span) +
    scale_color_viridis(option = "plasma", direction = -1, name = fea) +
    coord_cartesian(xlim = xlim_q) +
    theme_bw(base_size = 14) +
    labs(
      title = paste("SHAP Dependency Plot:", fea),
      x = paste(fea, "(original value)"),
      y = "SHAP value"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    )
})
p_col1 <- wrap_plots(plots, ncol = 3)
ggsave(pdf_file_col1, p_col1, width = 15, height = 18)
p_grid <- wrap_plots(plots, ncol = grid_cols)
ggsave(pdf_file_grid, p_grid, width = 12, height = 8)