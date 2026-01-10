suppressPackageStartupMessages({
  library(glmnet)
  library(mboost)
  library(plsRglm)
  library(randomForestSRC)
  library(gbm)
  library(caret)
  library(e1071)
  library(xgboost)
  library(pROC)
  library(ComplexHeatmap)
  library(grid)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){ sink(tempfile()); on.exit(sink()) }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

standarize.fun <- function(indata, centerFlag, scaleFlag) {  
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }
  if (is.null(centerFlags)){ centerFlags = FALSE }
  if (length(centerFlags)==1){ centerFlags = rep(centerFlags, length(data)) }
  if (is.null(names(centerFlags))){ names(centerFlags) <- names(data) }
  if (is.null(scaleFlags)){ scaleFlags = FALSE }
  if (length(scaleFlags)==1){ scaleFlags = rep(scaleFlags, length(data)) }
  if (is.null(names(scaleFlags))){ names(scaleFlags) <- names(data) }
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = FALSE)
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, , drop=FALSE]
  return(outdata)
}

.get_pos_col <- function(colnames_vec){
  if ("1" %in% colnames_vec) "1" else colnames_vec[length(colnames_vec)]
}
.get_pos_level <- function(y){
  y <- as.factor(y)
  if ("1" %in% levels(y)) "1" else tail(levels(y), 1)
}

.safe_glmnet_cv_args <- function(y){
  y <- as.factor(y)
  if (length(levels(y)) < 2) return(list(nfolds=3, type.measure="deviance"))
  n0 <- sum(y==levels(y)[1]); n1 <- sum(y==levels(y)[2])
  max_folds <- max(3, min(10, n0, n1))
  type.measure <- if (min(n0, n1) >= 10 * max_folds) "auc" else "deviance"
  list(nfolds=max_folds, type.measure=type.measure)
}

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
    "glm"          = names(coef(fit)),
    "svm.formula"  = fit$subFeature,
    "train"        = fit$coefnames,
    "glmboost"     = names(coef(fit)[abs(coef(fit))>0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0],
    "rfsrc"        = {
      imp <- suppressWarnings(na.omit(fit$importance))
      if (length(imp) && any(imp > 0)) names(imp)[imp > 0] else fit$subFeature
    },
    "gbm"          = {
      sm <- summary.gbm(fit, n.trees = fit$best_iter %||% fit$n.trees, plotit = FALSE)
      rownames(sm)[sm$rel.inf>0]
    },
    "xgb.Booster"  = fit$subFeature,
    "naiveBayes"   = fit$subFeature
  ))
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(unique(Feature))
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature, drop=FALSE]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = as.numeric(predict(fit, type = 'response', newx = as.matrix(new_data))),
    "glm"          = as.numeric(predict(fit, type = 'response', newdata = as.data.frame(new_data))),
    "svm.formula"  = { pr <- predict(fit, as.data.frame(new_data), probability = TRUE); prob <- attr(pr, "probabilities"); as.numeric(prob[, .get_pos_col(colnames(prob))]) },
    "train"        = { pr <- predict(fit, new_data, type = "prob"); as.numeric(pr[[.get_pos_col(colnames(pr))]]) },
    "glmboost"     = as.numeric(predict(fit, type = "response", newdata = as.data.frame(new_data))),
    "plsRglmmodel" = as.numeric(predict(fit, type = "response", newdata = as.data.frame(new_data))),
    "rfsrc"        = { pr <- predict(fit, as.data.frame(new_data))$predicted; as.numeric(pr[, .get_pos_col(colnames(pr))]) },
    "gbm"          = as.numeric(predict(fit, newdata = as.data.frame(new_data), type = "response", n.trees = fit$best_iter %||% fit$n.trees)),
    "xgb.Booster"  = as.numeric(predict(fit, newdata = as.matrix(new_data))),
    "naiveBayes"   = { pr <- predict(fit, type = "raw", newdata = new_data); as.numeric(pr[, .get_pos_col(colnames(pr))]) }
  ))
  names(RS) <- rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  scr <- CalPredictScore(fit, new_data)
  thr <- fit$threshold %||% 0.5
  lab <- ifelse(scr > thr, "1", "0")
  names(lab) <- names(scr)
  as.character(lab)
}

.learn_threshold <- function(y, p){
  y <- as.factor(y)
  if (length(levels(y)) < 2) return(0.5)
  pos <- .get_pos_level(y)
  y <- factor(ifelse(y==pos, "1","0"), levels=c("0","1"))
  rocobj <- suppressMessages(roc(y, p, quiet = TRUE, direction = "<"))
  as.numeric(coords(rocobj, "best", ret = "threshold", best.method = "youden"))
}

RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  method <- gsub(" ", "", method)
  m <- regexec("^(\\w+)(?:\\[(.+)\\])?$", method)
  parts <- regmatches(method, m)[[1]]
  method_name <- parts[2]
  raw_param   <- parts[3]
  method_param <- switch(
    method_name,
    "Enet"    = if (!is.na(raw_param)) list(alpha = as.numeric(sub("^alpha=", "", raw_param))) else list(alpha = 0.5),
    "Stepglm" = list(direction = if (!is.na(raw_param)) raw_param else "both"),
    NULL
  )
  message(sprintf("Run %s for %s; params=%s; using %d variables",
                  method_name, mode,
                  if(length(method_param)) paste0(names(method_param),"=",unlist(method_param),collapse=",") else "none",
                  ncol(Train_set)))
  args <- c(list(Train_set = Train_set, Train_label = Train_label, mode = mode, classVar = classVar), method_param)
  obj <- do.call(what = paste0("Run", method_name), args = args) 
  if(mode == "Variable") message(length(obj), " Variables retained;\n") else message("\n")
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  y <- as.numeric(Train_label[[classVar]])
  cv_args <- .safe_glmnet_cv_args(y)
  cv.fit <- cv.glmnet(x = as.matrix(Train_set), y = y,
                      family = "binomial", alpha = alpha,
                      nfolds = cv_args$nfolds, type.measure = cv_args$type.measure)
  fit <- glmnet(x = as.matrix(Train_set), y = y,
                family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature <- colnames(Train_set)
  fit$threshold <- .learn_threshold(Train_label[[classVar]], as.numeric(predict(fit, type='response', newx=as.matrix(Train_set))))
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}
RunLasso <- function(Train_set, Train_label, mode, classVar) RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
RunRidge <- function(Train_set, Train_label, mode, classVar) RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction = "both"){
  dat <- cbind(as.data.frame(Train_set), ..y = as.numeric(Train_label[[classVar]]))
  fit <- step(glm(formula = ..y ~ ., family = "binomial", data = dat),
              direction = direction, trace = 0)
  fit$subFeature <- setdiff(colnames(model.matrix(fit)), "(Intercept)")
  fit$threshold <- .learn_threshold(Train_label[[classVar]], as.numeric(predict(fit, type='response')))
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- svm(formula = eval(parse(text = paste(classVar, "~."))), data= data, probability = TRUE)
  fit$subFeature = colnames(Train_set)
  pr <- predict(fit, data, probability = TRUE); prob <- attr(pr, "probabilities")
  fit$threshold <- .learn_threshold(Train_label[[classVar]], prob[, .get_pos_col(colnames(prob))])
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = train(eval(parse(text = paste(classVar, "~."))), 
              data = data, method="lda",
              trControl = trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  probs <- predict(fit, data, type="prob")
  fit$threshold <- .learn_threshold(Train_label[[classVar]], probs[[.get_pos_col(colnames(probs))]])
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(as.data.frame(Train_set), y = as.factor(Train_label[[classVar]]))
  fit <- glmboost(y ~ ., data = data, family = Binomial())
  k <- min(10, max(3, floor(nrow(data)/5)))
  cvm <- cvrisk(fit, folds = cv(model.weights(fit), type = "kfold", B = k), papply = lapply)
  fit <- fit[mstop(cvm)]
  fit$subFeature = colnames(Train_set)
  pr <- as.numeric(predict(fit, type="response"))
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  X <- as.data.frame(Train_set)
  Y <- as.numeric(Train_label[[classVar]])
  best_nt <- 2
  try({
    cvres <- cv.plsRglm(formula = I(Y) ~ ., data = cbind(X, Y=Y), nt = 10, verbose = FALSE)
    if (!is.null(cvres$cvres)) best_nt <- which.min(cvres$cvres)
  }, silent = TRUE)
  best_nt <- max(1, min(10, as.integer(best_nt)))
  fit <- try(plsRglm(Y, X, nt = best_nt, modele = "pls-glm-logistic", verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    dat <- cbind(X, ..y = Y)
    assign(".__pls_dat", dat, envir = .GlobalEnv)
    on.exit({ if (exists(".__pls_dat", envir=.GlobalEnv)) rm(".__pls_dat", envir=.GlobalEnv) }, add = TRUE)
    fit <- plsRglm(formula = ..y ~ ., data = get(".__pls_dat", envir=.GlobalEnv),
                   nt = best_nt, modele = "pls-glm-logistic", verbose = FALSE)
  }
  fit$subFeature <- colnames(Train_set)
  pr <- as.numeric(predict(fit, type="response", newdata=as.data.frame(Train_set)))
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  # If only used for variable preselection, return all features to avoid 0-length
  if (mode == "Variable") return(colnames(Train_set))
  rf_nodesize = 5
  y <- as.factor(Train_label[[classVar]])
  dat <- cbind(as.data.frame(Train_set), y = y)
  fit <- rfsrc(y ~ ., data = dat, ntree = 1000, nodesize = rf_nodesize,
               importance = TRUE, proximity = FALSE, forest = TRUE)
  fit$subFeature = colnames(Train_set)
  pr <- predict(fit, as.data.frame(Train_set))$predicted
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr[, .get_pos_col(colnames(pr))])
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_set, Train_label, mode, classVar){
  dat <- cbind(as.data.frame(Train_set), ..y = as.numeric(Train_label[[classVar]]))
  fit0 <- gbm(..y ~ ., data = dat, distribution = 'bernoulli',
              n.trees = 10000, interaction.depth = 3,
              n.minobsinnode = 10, shrinkage = 0.01,
              cv.folds = 10)
  best_iter <- gbm.perf(fit0, method="cv", plot.it=FALSE)
  fit <- fit0; fit$best_iter <- best_iter
  fit$subFeature = colnames(Train_set)
  pr <- as.numeric(predict(fit, type='response', newdata=as.data.frame(Train_set), n.trees = best_iter))
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr)
  if (mode == "Model") return(fit)
  if (mode == "Variable"){
    sm <- summary.gbm(fit, n.trees = best_iter, plotit = FALSE)
    return(rownames(sm)[sm$rel.inf>0])
  }
}

RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  y <- as.numeric(Train_label[[classVar]])
  idxs = createFolds(Train_label[[classVar]], k = 5, list=TRUE)
  CV <- vapply(idxs, function(pt){
    dtrain = xgb.DMatrix(data = as.matrix(Train_set[-pt, , drop=FALSE]), label = y[-pt])
    dtest  = xgb.DMatrix(data = as.matrix(Train_set[ pt, , drop=FALSE]), label = y[ pt])
    watchlist <- list(train=dtrain, test=dtest)
    bst <- xgb.train(data=dtrain, max.depth=2, eta=0.1, nrounds=400, nthread=2,
                     watchlist=watchlist, objective="binary:logistic", eval_metric="logloss", verbose=0)
    which.min(bst$evaluation_log$test_logloss)
  }, numeric(1))
  nround <- as.numeric(names(which.max(table(CV))))
  fit <- xgboost(data = as.matrix(Train_set), label = y,
                 max.depth = 2, eta = 0.1, nrounds = nround, nthread = 2,
                 objective = "binary:logistic", verbose = 0)
  fit$subFeature = colnames(Train_set)
  pr <- as.numeric(predict(fit, newdata=as.matrix(Train_set)))
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(fit$subFeature)
}

RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(as.data.frame(Train_set), y = as.factor(Train_label[[classVar]]))
  fit <- naiveBayes(y ~ ., data = data)
  fit$subFeature = colnames(Train_set)
  pr <- predict(object = fit, type = "raw", newdata = Train_set)
  fit$threshold <- .learn_threshold(Train_label[[classVar]], pr[, .get_pos_col(colnames(pr))])
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunEval <- function(fit, 
                    Test_set = NULL, 
                    Test_label = NULL, 
                    Train_set = NULL, 
                    Train_label = NULL, 
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  } 
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature, drop=FALSE],
                                 Test_set[,  fit$subFeature, drop=FALSE])
    if(!is.null(Train_name)) Train_label[[cohortVar]] <- Train_name else Train_label[[cohortVar]] <- "Training"
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1], 
                             levels = c(unique(Train_label[,cohortVar]),
                                        setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature, drop=FALSE]
  }
  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  Predict.out <- Test_label; Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])

  auc_vals <- vapply(Predict.out, function(data){
    y <- factor(data[[classVar]])
    pos <- .get_pos_level(y); y <- factor(ifelse(y==pos,"1","0"), levels=c("0","1"))
    as.numeric(auc(suppressMessages(roc(y, data$RS, quiet=TRUE, direction = "<"))))
  }, numeric(1))

  attr(auc_vals, "ci") <- lapply(Predict.out, function(data){
    y <- factor(data[[classVar]])
    pos <- .get_pos_level(y); y <- factor(ifelse(y==pos,"1","0"), levels=c("0","1"))
    suppressMessages(ci.auc(roc(y, data$RS, quiet=TRUE, direction = "<")))
  })
  return(auc_vals)
}
NestedCV_single <- function(method,
                            Train_set,
                            Train_class,
                            classVar = "outcome",
                            K_outer = 5,
                            seed = 777) {
  set.seed(seed)
  folds <- caret::createFolds(Train_class[[classVar]],
                              k = K_outer,
                              list = TRUE,
                              returnTrain = FALSE)
  
  auc_vec <- rep(NA_real_, length(folds))
  names(auc_vec) <- paste0("Fold", seq_along(folds))
  
  for (k in seq_along(folds)) {
    message("Nested CV - outer fold ", k, " / ", length(folds), " for ", method)
    
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(nrow(Train_set)), test_idx)
    
    X_tr <- Train_set[train_idx, , drop = FALSE]
    Y_tr <- Train_class[train_idx, , drop = FALSE]
    X_te <- Train_set[test_idx,  , drop = FALSE]
    Y_te <- Train_class[test_idx, , drop = FALSE]
    parts <- strsplit(method, "\\+")[[1]]
    if (length(parts) == 1) parts <- c("simple", parts)
    var_method   <- parts[1]
    model_method <- parts[2]
    
    if (var_method == "simple") {
      vars <- colnames(X_tr)
    } else {
      vars <- tryCatch(
        RunML(method     = var_method,
              Train_set   = X_tr,
              Train_label = Y_tr,
              mode        = "Variable",
              classVar    = classVar),
        error = function(e) NULL
      )
      if (is.null(vars) || length(vars) < 1) {
        vars <- colnames(X_tr)
      }
    }
    
    X_tr_sub <- X_tr[, vars, drop = FALSE]
    X_te_sub <- X_te[, vars, drop = FALSE]
    fit <- tryCatch(
      RunML(method     = model_method,
            Train_set   = X_tr_sub,
            Train_label = Y_tr,
            mode        = "Model",
            classVar    = classVar),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    p_te <- tryCatch(
      CalPredictScore(fit, X_te_sub),
      error = function(e) NULL
    )
    if (is.null(p_te)) next
    
    y_te <- factor(Y_te[[classVar]])
    pos  <- .get_pos_level(y_te)
    y_bin <- factor(ifelse(y_te == pos, "1", "0"), levels = c("0", "1"))
    
    rocobj <- suppressMessages(pROC::roc(y_bin, p_te, quiet = TRUE, direction = "<"))
    auc_vec[k] <- as.numeric(pROC::auc(rocobj))
  }
  
  return(auc_vec)
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex, 
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = FALSE)
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = TRUE, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = FALSE)
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = c("#4195C1", "#FFFFFF", "#CB5746"),
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_columns = cluster_columns, cluster_rows = cluster_rows,
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) {
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}
