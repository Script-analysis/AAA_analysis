set.seed(777)
work.path <- "D:\\ML"
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
Train_expr <- read.table(file.path(data.path, "Training_expr_AAA.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
Train_expr <- as.data.frame(Train_expr)
Train_class <- read.table(file.path(data.path, "Training_class_AAA.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
Test_expr <- read.table(file.path(data.path, "Testing_expr_AAA.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
Test_expr <- as.data.frame(Test_expr)
Test_class <- read.table(file.path(data.path, "Testing_class_AAA.txt"), header = TRUE, sep = "\t", row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
# Align samples
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = FALSE]
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = FALSE]
# Align genes
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,])
Test_expr <- t(Test_expr[comgene,])
# Standardize by cohort (no leakage)
Train_set = scaleData(data = Train_expr, centerFlags = FALSE, scaleFlags = FALSE) 
invisible(names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort)))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = FALSE, scaleFlags = FALSE)
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
classVar = "outcome"
min.selected.var = 5
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.var <- list()
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
## =========================================================
final_method <- "Stepglm[both]+RF"     
nested_auc_vec <- NestedCV_single(
  method     = final_method,
  Train_set  = Train_set,
  Train_class = Train_class,
  classVar   = classVar,
  K_outer    = 5,
  seed       = 777
)
nested_auc_summary <- data.frame(
  Model   = final_method,
  MeanAUC = mean(nested_auc_vec, na.rm = TRUE),
  SDAUC   = sd(nested_auc_vec,   na.rm = TRUE)
)
write.table(t(nested_auc_vec),
            file = file.path(res.path, "NestedCV_AMRMS_Folds.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(nested_auc_summary,
            file = file.path(res.path, "NestedCV_AMRMS_Summary.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
## =========================================================
# --- Evaluate on training+external cohorts ---
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
                     Train_name = "GSE57691",
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
# --- Plot heatmap ---
AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), , drop=FALSE]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

if(ncol(AUC_mat) < 3) { CohortCol <- c("red","blue") } else { CohortCol <- brewer.pal(n = max(3, ncol(AUC_mat)), name = "Paired") }
names(CohortCol) <- colnames(AUC_mat)

hm <- SimpleHeatmap(AUC_mat,
                    avg_AUC,
                    CohortCol, "steelblue",
                    cellwidth = 1, cellheight = 0.5,
                    cluster_columns = FALSE, cluster_rows = FALSE)

pdf(file.path(fig.path, "AUC.pdf"), width = 8, height = 22)
draw(hm)
invisible(dev.off())
########################
AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),
                      sep = "\t", row.names = 1, header = TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE)
avg_AUC <- apply(AUC_mat, 1, mean, na.rm = TRUE)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
top10_names <- names(avg_AUC)[1:min(10, length(avg_AUC))]
AUC_mat <- AUC_mat[top10_names, , drop = FALSE]
avg_AUC_show <- as.numeric(format(avg_AUC[top10_names], digits = 3, nsmall = 3))
selected_models <- intersect(top10_names, names(use_model))
if (length(selected_models) < 3) {
  selected_models <- names(use_model)
}
library(pROC)
library(dplyr)
library(openxlsx)
selected_models <- intersect(top10_names, names(use_model))
if (length(selected_models) < 3) {
  selected_models <- names(use_model)
}
print(selected_models)
summarise_one_model <- function(m) {
  fit <- use_model[[m]]
  if (is.null(fit)) return(NULL)
  p_tr <- CalPredictScore(fit, Train_set)
  p_te <- CalPredictScore(fit, Test_set)
  
  # 对齐顺序
  p_tr <- p_tr[rownames(Train_class)]
  p_te <- p_te[rownames(Test_class)]
 y_tr <- as.integer(Train_class[[classVar]])
  y_te <- as.integer(Test_class[[classVar]])
  metrics_one <- function(prob, y, dataset_label) {
    rocobj <- pROC::roc(y, prob, quiet = TRUE, direction = "<")
    ci     <- pROC::ci.auc(rocobj)
    youden <- pROC::coords(
      rocobj,
      x = "best",
      best.method = "youden",
      ret = c("threshold","sensitivity","specificity")
    )
    data.frame(
      Model       = m,
      Dataset     = dataset_label,          
      AUC         = as.numeric(pROC::auc(rocobj)),
      AUC_L       = as.numeric(ci[1]),
      AUC_U       = as.numeric(ci[3]),
      Threshold   = as.numeric(youden["threshold"]),
      Sensitivity = as.numeric(youden["sensitivity"]),
      Specificity = as.numeric(youden["specificity"]),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    metrics_one(p_tr, y_tr, "Training"),
    metrics_one(p_te, y_te, "External")
  )
}
metric_list <- lapply(selected_models, summarise_one_model)
ML_performance <- bind_rows(metric_list)
num_cols <- c("AUC","AUC_L","AUC_U","Threshold","Sensitivity","Specificity")
ML_performance[num_cols] <- lapply(ML_performance[num_cols],
                             function(x) round(x, digits = 3))

write.xlsx(ML_performance,file = file.path(res.path, "ML_performance.xlsx"),
  sheetName = "ML_performance",
  overwrite = TRUE
)


model <- readRDS(file.path(res.path, "model.rds"))
library(fastshap)
library(ggplot2)
library(reshape2)   
rf_fit  <- model[["Stepglm[both]+RF"]]
stopifnot(!is.null(rf_fit))
feat <- rf_fit$subFeature
X_train <- as.data.frame(Train_set[, feat, drop = FALSE])
pred_wrapper_rfsrc <- function(object, newdata) {
  newdata <- as.data.frame(newdata)[, object$subFeature, drop = FALSE]
  pr <- predict(object, newdata)$predicted
  pos_col <- if ("1" %in% colnames(pr)) "1" else colnames(pr)[ncol(pr)]
  as.numeric(pr[, pos_col])
}
set.seed(2025)
shap_vals <- fastshap::explain(
  object = rf_fit,
  X = X_train,
  pred_wrapper = pred_wrapper_rfsrc,
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
  labs(title = "Global SHAP Importance (Stepglm + RF)",
       x = "Variables", y = "Mean |SHAP|") +
  theme_minimal(base_size = 14)
ggsave("SHAP Importance.pdf",height=5,width=5)  
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(scales)
shap_long <- melt(shap_vals, varnames = c("sample", "feature"), value.name = "phi")
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
  group_by(feature) %>%
  summarise(mean_abs_phi = mean(abs(phi), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_phi)) %>%
  pull(feature)

shap_long$feature <- factor(shap_long$feature, levels = rev(feature_order)) 
p_beeswarm <- ggplot(shap_long, aes(x = phi, y = feature, color = group_num)) +
  geom_quasirandom(groupOnX = FALSE, alpha = 0.8, size = 1.0) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name   = "Feature value",
                        breaks = c(0, 0.5, 1),
                        labels = c("Low", "Mid", "High")) +
  theme_classic(base_size = 14) +
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
#####
library(pROC)
library(ggplot2)
roc_aplnr_tr <- roc(as.integer(Train_class$outcome),
                    as.numeric(Train_set[, "APLNR"]), direction = "<")
roc_aplnr_te <- roc(as.integer(Test_class$outcome),
                    as.numeric(Test_set[, "APLNR"]), direction = "<")
roc_fkbp11_tr <- roc(as.integer(Train_class$outcome),
                     as.numeric(Train_set[, "FKBP11"]), direction = "<")
roc_fkbp11_te <- roc(as.integer(Test_class$outcome),
                     as.numeric(Test_set[, "FKBP11"]), direction = "<")
df_plot <- rbind(
  data.frame(specificity = rev(roc_aplnr_tr$specificities),
             sensitivity = rev(roc_aplnr_tr$sensitivities),
             Gene = "APLNR", Dataset = "Train"),
  data.frame(specificity = rev(roc_aplnr_te$specificities),
             sensitivity = rev(roc_aplnr_te$sensitivities),
             Gene = "APLNR", Dataset = "Test"),
  data.frame(specificity = rev(roc_fkbp11_tr$specificities),
             sensitivity = rev(roc_fkbp11_tr$sensitivities),
             Gene = "FKBP11", Dataset = "Train"),
  data.frame(specificity = rev(roc_fkbp11_te$specificities),
             sensitivity = rev(roc_fkbp11_te$sensitivities),
             Gene = "FKBP11", Dataset = "Test")
)
auc_labels <- data.frame(
  Gene    = c("APLNR", "APLNR", "FKBP11", "FKBP11"),
  Dataset = c("Train", "Test", "Train", "Test"),
  label   = c(
    sprintf("APLNR Train AUC=%.3f", auc(roc_aplnr_tr)),
    sprintf("APLNR Test  AUC=%.3f", auc(roc_aplnr_te)),
    sprintf("FKBP11 Train AUC=%.3f", auc(roc_fkbp11_tr)),
    sprintf("FKBP11 Test  AUC=%.3f", auc(roc_fkbp11_te))
  )
)
p <- ggplot(df_plot, aes(x = 1 - specificity, y = sensitivity,
                         color = Gene, linetype = Dataset)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_bw(base_size = 14) +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "Univariate ROC curves: APLNR vs FKBP11") +
  scale_color_manual(values = c("APLNR" = "#1f77b4", "FKBP11" = "#d62728")) +
  scale_linetype_manual(values = c("Train" = "solid", "Test" = "dashed")) +
  theme(legend.position = "bottom")

print(p)

ggsave("APLNR_vs_FKBP11_ROC.pdf", p, width =5, height = 5)