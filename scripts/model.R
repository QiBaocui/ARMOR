#!/usr/bin/env Rscript

# ==============================================================================
# File: model.R
# Description: Machine learning model construction functions for AMR prediction
#              Based on H2O framework with multiple algorithms:
#              - GBM (Gradient Boosting Machine)
#              - GLM (Generalized Linear Model)
#              - RF (Random Forest)
#              - DL (Deep Learning)
#              - Stacked Ensemble
#
# Functions:
#   splitSample()      - Split data into training and testing sets
#   get_outputperformance() - Extract optimal performance metrics at Youden cutoff
#   get_predict()      - Generate predictions, ROC curves, and variable importance
#   modelFun()         - Train multiple models and save results
#
# Author: ARMOR Pipeline
# Reference: Manuscript.docx - Methods section
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: splitSample
# Description: Split samples into training and testing sets
#
# Parameters:
#   group       - Data frame with sample IDs and group labels
#   split_type  - Splitting method: "proportional" or "sample_size"
#   ratio       - Training set proportion (for proportional method)
#   num         - Total sample size (for sample_size method)
#
# Returns:
#   List containing train_group and test_group data frames
# ------------------------------------------------------------------------------
splitSample <- function(group, split_type = "proportional", ratio = 0.6, num = 100) {
  # Ensure proper column names
  colnames(group) <- c("ID", "Group")

  if (split_type == "proportional") {
    # Proportional split: maintain class distribution in both sets
    library(caret)
    train_index <- createDataPartition(group[, 2], p = ratio, list = FALSE)
    train_group <- group[train_index, ]
    test_group <- group[-train_index, ]

  } else if (split_type == "sample_size") {
    # Fixed sample size split: equal number from each group
    per_group_num <- floor(num / 2)  # Ensure equal numbers per group

    group_levels <- unique(group$Group)
    if (length(group_levels) != 2) {
      stop("split_type='sample_size' only supports 2 groups! Current: ",
           length(group_levels))
    }

    train_list <- list()
    for (i in 1:length(group_levels)) {
      # Get samples from current group
      sub_group <- group[group$Group == group_levels[i], ]

      # Check if enough samples
      if (nrow(sub_group) < per_group_num) {
        stop("Group '", group_levels[i], "' has insufficient samples (",
             nrow(sub_group), " < ", per_group_num, ")")
      }

      # Random sampling without replacement
      sub_train_index <- sample(1:nrow(sub_group), size = per_group_num, replace = FALSE)
      train_list[[i]] <- sub_group[sub_train_index, ]
    }

    # Combine training samples
    train_group <- do.call(rbind, train_list)

    # Test set: samples not in training
    test_group <- group[!group$ID %in% train_group$ID, ]

    # Reset row names
    rownames(train_group) <- NULL
    rownames(test_group) <- NULL

  } else {
    stop("Invalid split_type. Use 'proportional' or 'sample_size'")
  }

  df <- list(train_group = train_group, test_group = test_group)
  return(df)
}


# ------------------------------------------------------------------------------
# Function: get_outputperformance
# Description: Extract performance metrics at optimal cutoff (Youden index)
#
# Parameters:
#   perf - H2O performance object
#
# Returns:
#   Data frame with performance metrics at optimal threshold
#
# Method:
#   Youden Index = Sensitivity + Specificity - 1 = TPR - FPR
#   Optimal cutoff maximizes Youden Index
# ------------------------------------------------------------------------------
get_outputperformance <- function(perf) {
  # Extract thresholds and metrics from H2O performance object
  outputperformance <- perf %>%
    slot("@metrics") %>%
    extract2("thresholds_and_metric_scores")

  outputperformance.R <- as.data.frame(outputperformance)

  # Calculate Youden Index for each threshold
  youdenIndex <- outputperformance.R$tpr - outputperformance.R$fpr

  # Find optimal cutoff (maximum Youden Index)
  cutoff <- max(youdenIndex)
  site <- which(youdenIndex == cutoff)

  # Get performance metrics at optimal cutoff
  cutoffpermance <- outputperformance.R[site, ]

  return(cutoffpermance)
}


# ------------------------------------------------------------------------------
# Function: get_predict
# Description: Generate predictions, ROC curves, and save results
#
# Parameters:
#   model      - Trained H2O model
#   modelName  - Model name/identifier
#   test.R     - Test set as R data frame
#   test.h2o   - Test set as H2O frame
#   outdir     - Output directory
#   feature    - Feature type (Gene/SNP/Kmer/All)
#
# Outputs:
#   - ROC curve plot (PNG and PDF)
#   - Predictions TSV
#   - Performance metrics TSV
#   - ROC data TSV
#   - Variable importance TSV
# ------------------------------------------------------------------------------
get_predict <- function(model, modelName, test.R, test.h2o, outdir, feature) {
  # Generate predictions
  pred_matri <- h2o.predict(object = model, newdata = test.h2o) %>%
    as.data.frame()

  # Add sample info
  pred_matri <- cbind(Group = test.R[, 1], pred_matri)
  pred_matri <- cbind(sample_code = rownames(test.R), pred_matri)

  # Calculate performance metrics
  test_perf <- h2o.performance(model = model, newdata = test.h2o)
  cutoffpermance <- get_outputperformance(test_perf)

  # Add AUC
  cutoffpermance$auc <- h2o.auc(test_perf)
  cutoffpermance$model <- modelName

  # Extract ROC data
  rocData <- test_perf %>%
    slot("@metrics") %>%
    extract2("thresholds_and_metric_scores") %>%
    extract(c("tpr", "fpr")) %>%
    add_row(tpr = 0, fpr = 0, .before = TRUE)  # Add origin point

  # Sort by FPR for proper ROC curve
  rocData <- rocData[order(rocData$fpr), ]
  rocData$label <- paste0(modelName, "_AUC: ", signif(h2o.auc(test_perf), digits = 3))

  # Create ROC curve plot
  rocCurve <- ggplot(rocData, aes(x = fpr, y = tpr, color = label)) +
    geom_line(size = 1.2) +
    theme_bw() +
    xlab("False positive rate") +
    ylab("True positive rate") +
    labs(color = '') +
    theme(
      legend.position = c(.8, .2),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +  # Diagonal reference
    ggtitle('ROC Curve')

  # Calculate variable importance
  if (modelName == "stacked") {
    # For stacked ensemble, use permutation importance
    var <- h2o.permutationImportance(model, newdata = test.h2o) %>%
      as.data.frame()
  } else {
    # For other models, use built-in variable importance
    var <- h2o.varimp(model) %>%
      as.data.frame()
  }

  # Set bitmap type for PNG output
  options(bitmapType = 'cairo')

  # Save ROC curve plot
  png(paste0(outdir, "/", feature, "_", modelName, "_ROC.png"),
      width = 800, height = 800)
  print(rocCurve)  # Must explicitly print in device context
  dev.off()

  # Save PDF version
  ggsave(paste0(outdir, "/", feature, "_", modelName, "_ROC.pdf"),
         rocCurve, width = 8, height = 8)

  # Save predictions
  write.table(pred_matri,
              paste0(outdir, "/", feature, "_", modelName, "_pred.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")

  # Save performance metrics at optimal cutoff
  write.table(cutoffpermance,
              paste0(outdir, "/", feature, "_", modelName, "_cutoffpermance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")

  # Save ROC data
  write.table(rocData,
              paste0(outdir, "/", feature, "_", modelName, "_rocData.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")

  # Save variable importance
  write.table(var,
              paste0(outdir, "/", feature, "_", modelName, "_varimp.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")
}


# ------------------------------------------------------------------------------
# Function: modelFun
# Description: Train multiple ML models and generate predictions
#
# Parameters:
#   train.R           - Training set as R data frame
#   test.R            - Test set as R data frame
#   outdir            - Output directory
#   feature           - Feature type name (Gene/SNP/Kmer/All)
#   feature_type_list - List of algorithms to train
#   nfolds            - Number of cross-validation folds
#
# Supported algorithms:
#   - GBM: Gradient Boosting Machine
#   - GLM: Generalized Linear Model (Elastic Net)
#   - RF: Random Forest
#   - DL: Deep Learning (Neural Network)
#   - xGboost: XGBoost
#   - Bayes: Naive Bayes
#   - stacked: Stacked Ensemble (combines all above)
# ------------------------------------------------------------------------------
modelFun <- function(train.R, test.R, outdir, feature,
                     feature_type_list = c("GBM", "GLM", "RF", "DL"),
                     nfolds = 10) {

  # Initialize H2O cluster (should be done externally, but check)
  # h2o.init()

  # Convert to H2O frames
  train.h2o <- as.h2o(train.R)
  train.h2o$Group <- h2o.asfactor(train.h2o$Group)  # Convert to factor

  test.h2o <- as.h2o(test.R)
  test.h2o$Group <- h2o.asfactor(test.h2o$Group)

  # Define response and predictor variables
  y <- "Group"
  x <- setdiff(names(train.R), y)

  model_list <- list()

  # Train GBM model
  if ("GBM" %in% feature_type_list) {
    message(sprintf("[%s] Training GBM model...", feature))
    model_list[["GBM"]] <- h2o.gbm(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      keep_cross_validation_fold_assignment = TRUE,
      seed = 1
    )
    get_predict(model = model_list[["GBM"]], modelName = "GBM",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train GLM model (Elastic Net)
  if ("GLM" %in% feature_type_list) {
    message(sprintf("[%s] Training GLM model...", feature))
    model_list[["GLM"]] <- h2o.glm(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      keep_cross_validation_fold_assignment = TRUE,
      seed = 1
    )
    get_predict(model = model_list[["GLM"]], modelName = "GLM",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train Random Forest
  if ("RF" %in% feature_type_list) {
    message(sprintf("[%s] Training Random Forest model...", feature))
    model_list[["RF"]] <- h2o.randomForest(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      seed = 1
    )
    get_predict(model = model_list[["RF"]], modelName = "RF",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train Deep Learning model
  if ("DL" %in% feature_type_list) {
    message(sprintf("[%s] Training Deep Learning model...", feature))
    model_list[["DL"]] <- h2o.deeplearning(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      seed = 1
    )
    get_predict(model = model_list[["DL"]], modelName = "DL",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train XGBoost model
  if ("xGboost" %in% feature_type_list) {
    message(sprintf("[%s] Training XGBoost model...", feature))
    model_list[["xGboost"]] <- h2o.xgboost(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      seed = 1
    )
    get_predict(model = model_list[["xGboost"]], modelName = "xGboost",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train Naive Bayes model
  if ("Bayes" %in% feature_type_list) {
    message(sprintf("[%s] Training Naive Bayes model...", feature))
    model_list[["Bayes"]] <- h2o.naiveBayes(
      x = x, y = y,
      training_frame = train.h2o,
      nfolds = nfolds,
      keep_cross_validation_predictions = TRUE,
      fold_assignment = "Stratified",
      seed = 1
    )
    get_predict(model = model_list[["Bayes"]], modelName = "Bayes",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)
  }

  # Train Stacked Ensemble
  if (length(model_list) >= 2) {
    message(sprintf("[%s] Training Stacked Ensemble model...", feature))
    stacked <- h2o.stackedEnsemble(
      x = x, y = y,
      training_frame = train.h2o,
      base_models = model_list
    )
    get_predict(model = stacked, modelName = "stacked",
                test.R = test.R, test.h2o = test.h2o,
                outdir = outdir, feature = feature)

    # Save stacked model
    h2o.saveModel(object = stacked,
                  path = outdir,
                  filename = paste0(feature, "_stacked_model"),
                  force = TRUE)
  }

  # Save all base models
  for (i in names(model_list)) {
    h2o.saveModel(object = model_list[[i]],
                  path = outdir,
                  filename = paste0(feature, "_", i, "_model"),
                  force = TRUE)
  }

  message(sprintf("[%s] Model training complete!", feature))
}
