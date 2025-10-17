###########--------------------------------------------------------------------###########
#          Random Forest Validation of WGCNA Module Assignments
#           Salt Stress Response in Chlamydomonas reinhardtii
#
#  Script for: Integrating co-expression network analysis and machine learning to reveal 
#  the regulatory landscape of GPD genes in Chlamydomonas reinhardtii under salinity stress
#
#  This script validates WGCNA module assignments using Random Forest classification
#  with comprehensive performance metrics including ROC curves, confusion matrices,
#  and UMAP dimensionality reduction visualizations.
#
#
#  Authors: Tzec-Interián et al.
#  Date: October 2025
#
###########--------------------------------------------------------------------###########

condition_prefix <- "salt"
date_export <- format(Sys.Date(), "%Y%m%d")
ml_prefix <- "ML_RandomForest"

###########----------------------- 0. Package Installation and Loading ---------------

required_packages <- c(
  "randomForest", "pROC", "caret", "dplyr", "ggplot2", 
  "reshape2", "umap"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(randomForest)
library(pROC)
library(caret)
library(dplyr)
library(ggplot2)
library(reshape2)
library(umap)

###########----------------------- 1. Data Loading and Preparation ---------------

data_vst <- read.csv(paste0(condition_prefix, "_VST_corrected_fulldata.csv"), 
                     header = TRUE, row.names = 1)

moduleColors_df <- read.csv(paste0(condition_prefix, "_gene_modules.csv"), 
                            header = TRUE, row.names = 1)

# Merge expression data with module assignments
# Format: genes as observations (rows), samples as features (columns)
data_vst$Gene_ID <- rownames(data_vst)

moduleColors_df$Gene_ID <- rownames(moduleColors_df)
data_for_ml <- merge(data_vst, moduleColors_df, by = "Gene_ID")
colnames(data_for_ml)[colnames(data_for_ml) == "module"] <- "Module"

rownames(data_for_ml) <- data_for_ml$Gene_ID
gene_ids_original <- data_for_ml$Gene_ID
data_for_ml$Gene_ID_Tracker <- gene_ids_original
data_for_ml$Module <- as.factor(data_for_ml$Module)

data_for_ml$Gene_ID <- NULL

message("Data loaded successfully.")
message("Total genes: ", nrow(data_for_ml))
message("Total samples: ", ncol(data_vst))
message("Number of modules: ", length(unique(data_for_ml$Module)))

###########----------------------- 2. Data Splitting and Class Balancing ---------------

# Split data into training (80%) and test (20%) sets
# Apply upsampling to balance class distribution in training set

set.seed(123)

trainIndex <- createDataPartition(data_for_ml$Module, p = 0.8, list = FALSE)
data_train <- data_for_ml[trainIndex, ]
data_test <- data_for_ml[-trainIndex, ]

gene_ids_in_test_set <- data_test$Gene_ID_Tracker

message("Training set size: ", nrow(data_train), " genes")
message("Test set size: ", nrow(data_test), " genes")

predictors_train <- data_train[, !(names(data_train) %in% c("Module"))]
class_labels_train <- data_train$Module

data_balanced <- caret::upSample(
  x = predictors_train,
  y = class_labels_train,
  yname = "Module"
)

train_df <- data_balanced
train_model_df <- train_df[, !(names(train_df) == "Gene_ID_Tracker")]
test_df <- data_test

message("Balanced training set size: ", nrow(train_df), " genes")
message("\nModule distribution in balanced training set:")
print(train_df %>% count(Module))

###########----------------------- 3. Random Forest Model Training ---------------

# Train Random Forest classifier with 5-fold cross-validation
# Model parameters: 1500 trees, default mtry optimized via grid search

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  verboseIter = TRUE,
  savePredictions = "final"
)

set.seed(123)
modelo_rf <- train(
  Module ~ ., 
  data = train_model_df, 
  method = "rf", 
  trControl = ctrl, 
  metric = "Accuracy",
  ntree = 1500
)

message("Random Forest model trained successfully.")
print(modelo_rf)

###########----------------------- 4. Model Performance Evaluation ---------------

predictions_prob <- predict(modelo_rf, newdata = test_df, type = "prob")
predicted_classes <- predict(modelo_rf, newdata = test_df, type = "raw")

cm <- confusionMatrix(predicted_classes, test_df$Module)
print(cm)

# Calculate Area Under the Curve (AUC) for each module using one-vs-rest approach
auc_per_module <- sapply(colnames(predictions_prob), function(module) {
  binary_labels <- as.numeric(test_df$Module == module)
  if (length(unique(binary_labels)) == 2) {
    as.numeric(auc(roc(binary_labels, predictions_prob[[module]], quiet = TRUE)))
  } else {
    NA
  }
})

message("\nAUC per module:")
print(round(auc_per_module, 3))

colnames(predictions_prob) <- levels(test_df$Module)
roc_curves <- multiclass.roc(test_df$Module, predictions_prob)
message("Global AUC: ", round(auc(roc_curves), 3))

###########----------------------- 5. ROC Curve Analysis ---------------

predictions_prob_df <- as.data.frame(predictions_prob)
valid_modules <- intersect(colnames(predictions_prob_df), unique(test_df$Module))

roc_data <- data.frame()

for (module in valid_modules) {
  if (!(module %in% test_df$Module)) next
  true_labels <- as.numeric(test_df$Module == module)
  if (length(unique(true_labels)) < 2) next
  pred_scores <- predictions_prob_df[[module]]
  if (all(is.na(pred_scores))) next
  roc_obj <- roc(response = true_labels, predictor = pred_scores, quiet = TRUE)
  df_temp <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    module = module
  )
  roc_data <- rbind(roc_data, df_temp)
}

moduleColors <- setNames(moduleColors_df$module, rownames(moduleColors_df))
module_color_palette <- structure(unique(moduleColors), names = unique(moduleColors))
roc_curve_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = module)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = module_color_palette) +
  labs(
    x = "False Positive Rate (FPR)",
    y = "True Positive Rate (TPR)",
    color = "Module"
  ) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(roc_curve_plot)

ggsave(
  filename = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_ROC_curves.pdf"), 
  plot = roc_curve_plot, 
  width = 8, 
  height = 6
)

###########----------------------- 6. Confusion Matrix Visualization ---------------

cm_table <- as.data.frame(cm$table)
colnames(cm_table) <- c("Prediction", "Reference", "Freq")

heatmap_cm_plot <- ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = Freq), 
    size = 3.5, 
    color = ifelse(cm_table$Freq > max(cm_table$Freq)/2, "white", "black")
  ) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(
    x = "True Class",
    y = "Predicted Class",
    fill = "Frequency"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(heatmap_cm_plot)

ggsave(
  filename = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_confusion_matrix.pdf"), 
  plot = heatmap_cm_plot, 
  width = 8, 
  height = 8
)

###########----------------------- 7. Data Preparation for UMAP Analysis ---------------

# Prepare scaled features for dimensionality reduction
# Remove duplicates and standardize expression values

features_test <- test_df[, -which(names(test_df) == "Module")]
labels_test <- test_df$Module
preds_test <- predicted_classes

sample_columns <- 1:22
features_test_scaled <- scale(features_test[, sample_columns])
features_test_scaled[is.na(features_test_scaled)] <- 0

unique_feature_indices <- !duplicated(features_test[, sample_columns])
features_test_unique_scaled <- features_test_scaled[unique_feature_indices, , drop = FALSE]
labels_test_unique <- labels_test[unique_feature_indices]
preds_test_unique <- preds_test[unique_feature_indices]
item_ids_unique <- gene_ids_in_test_set[unique_feature_indices]

predicted_correctly_factor <- factor(
  ifelse(preds_test_unique == labels_test_unique, "Correct", "Incorrect"),
  levels = c("Incorrect", "Correct")
)

message("Correct predictions: ", sum(predicted_correctly_factor == "Correct"))
message("Incorrect predictions: ", sum(predicted_correctly_factor == "Incorrect"))

unique_cols <- unique(as.character(moduleColors))
paleta_modulos <- structure(unique_cols, names = unique_cols)
niveles_test_module <- levels(test_df$Module)
paleta_modulos_filtrada <- paleta_modulos[names(paleta_modulos) %in% niveles_test_module]

###########----------------------- 8. UMAP: All Modules ---------------

set.seed(42)
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_result <- umap(features_test_unique_scaled, config = umap_config)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Real_Module = labels_test_unique,
  Predicted_Correctly = predicted_correctly_factor
)
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Real_Module, shape = Predicted_Correctly)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = paleta_modulos_filtrada, name = "True Module") +
  scale_shape_manual(
    values = c("Incorrect" = 4, "Correct" = 16), 
    labels = c("Incorrect", "Correct"),
    name = "Prediction"
  ) +
  labs(
    x = "UMAP Dimension 1", 
    y = "UMAP Dimension 2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(umap_plot)

ggsave(
  filename = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_UMAP_all_modules.pdf"), 
  plot = umap_plot, 
  width = 8, 
  height = 6
)

###########----------------------- 9. UMAP: GPD Modules ---------------

# Focus on modules containing glycerol-3-phosphate dehydrogenase (GPD) genes
# These modules are of primary biological interest in salt stress response

gpd_modules <- c("black", "green", "turquoise", "magenta")
gpd_modules_in_data <- intersect(gpd_modules, levels(labels_test_unique))

message("GPD modules found in test set: ", paste(gpd_modules_in_data, collapse = ", "))

gpd_indices <- labels_test_unique %in% gpd_modules_in_data
features_gpd_scaled <- features_test_unique_scaled[gpd_indices, , drop = FALSE]
labels_gpd <- factor(
  as.character(labels_test_unique[gpd_indices]),
  levels = gpd_modules_in_data
)
preds_gpd <- preds_test_unique[gpd_indices]

predicted_correctly_gpd <- factor(
  as.character(preds_gpd) == as.character(labels_gpd),
  levels = c(FALSE, TRUE), 
  labels = c("Incorrect", "Correct")
)

paleta_gpd <- paleta_modulos_filtrada[names(paleta_modulos_filtrada) %in% gpd_modules_in_data]
paleta_gpd <- paleta_gpd[gpd_modules_in_data]

set.seed(42)
umap_config_gpd <- umap.defaults
umap_config_gpd$n_neighbors <- 5
umap_gpd_result <- umap(features_gpd_scaled, config = umap_config_gpd)

umap_df_gpd <- data.frame(
  UMAP1 = umap_gpd_result$layout[, 1],
  UMAP2 = umap_gpd_result$layout[, 2],
  Real_Module = labels_gpd,
  Predicted_Correctly = predicted_correctly_gpd
)
umap_plot_gpd <- ggplot(umap_df_gpd, aes(x = UMAP1, y = UMAP2, color = Real_Module, shape = Predicted_Correctly)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = paleta_gpd, name = "True Module (GPD)") +
  scale_shape_manual(
    values = c("Incorrect" = 4, "Correct" = 16), 
    labels = c("Incorrect", "Correct"),
    name = "Prediction"
  ) +
  labs(
    x = "UMAP Dimension 1", 
    y = "UMAP Dimension 2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(umap_plot_gpd)

ggsave(
  filename = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_UMAP_GPD_modules.pdf"), 
  plot = umap_plot_gpd, 
  width = 8, 
  height = 6
)

###########----------------------- 10. Misclassification Analysis ---------------

# Identify and characterize misclassified genes
# Focus on GPD-relevant modules (magenta and black)

misclassified_mask <- preds_test_unique != labels_test_unique
genes_misclassified_df <- data.frame(
  Gene_ID = item_ids_unique[misclassified_mask],
  Real_Module = labels_test_unique[misclassified_mask],
  Predicted_Module = preds_test_unique[misclassified_mask]
)

total_genes_test <- length(labels_test_unique)
total_misclassified <- sum(misclassified_mask)
pct_misclassified <- round(total_misclassified/total_genes_test*100, 2)

message("\nMisclassification Summary:")
message("Total genes in test set: ", total_genes_test)
message("Total misclassified genes: ", total_misclassified)
message("Misclassification rate: ", pct_misclassified, "%")

misclassified_magenta_mask <- (labels_test_unique == "magenta") & (preds_test_unique != "magenta")
genes_misclassified_magenta <- data.frame(
  Gene_ID = item_ids_unique[misclassified_magenta_mask],
  Real_Module = labels_test_unique[misclassified_magenta_mask],
  Predicted_Module = preds_test_unique[misclassified_magenta_mask]
)

misclassified_black_mask <- (labels_test_unique == "black") & (preds_test_unique != "black")
genes_misclassified_black <- data.frame(
  Gene_ID = item_ids_unique[misclassified_black_mask],
  Real_Module = labels_test_unique[misclassified_black_mask],
  Predicted_Module = preds_test_unique[misclassified_black_mask]
)

magenta_total <- sum(labels_test_unique == "magenta")
black_total <- sum(labels_test_unique == "black")

message("\nGPD Module-Specific Misclassification:")
message("Magenta module:")
message("  Total genes: ", magenta_total)
message("  Misclassified: ", nrow(genes_misclassified_magenta))
message("  Accuracy: ", round((1 - nrow(genes_misclassified_magenta)/magenta_total)*100, 2), "%")

message("Black module:")
message("  Total genes: ", black_total)
message("  Misclassified: ", nrow(genes_misclassified_black))
message("  Accuracy: ", round((1 - nrow(genes_misclassified_black)/black_total)*100, 2), "%")
write.csv(
  genes_misclassified_df, 
  file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_misclassified_genes.csv"),
  row.names = FALSE
)

write.csv(
  genes_misclassified_magenta, 
  file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_misclassified_magenta.csv"),
  row.names = FALSE
)

write.csv(
  genes_misclassified_black, 
  file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_misclassified_black.csv"),
  row.names = FALSE
)

###########----------------------- 11. Model Performance Summary ---------------

# Export comprehensive performance metrics for publication
# Includes overall and per-module statistics

performance_summary <- data.frame(
  Metric = c("Overall Accuracy", "Kappa", "Mean AUC", "Sensitivity", "Specificity"),
  Value = c(
    round(cm$overall["Accuracy"], 4),
    round(cm$overall["Kappa"], 4),
    round(mean(auc_per_module, na.rm = TRUE), 4),
    round(mean(cm$byClass[, "Sensitivity"], na.rm = TRUE), 4),
    round(mean(cm$byClass[, "Specificity"], na.rm = TRUE), 4)
  )
)

write.csv(
  performance_summary, 
  file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_performance_summary.csv"),
  row.names = FALSE
)

auc_df <- data.frame(
  Module = names(auc_per_module),
  AUC = round(auc_per_module, 4)
)

write.csv(
  auc_df, 
  file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_AUC_per_module.csv"),
  row.names = FALSE
)

message("\nPerformance metrics exported successfully.")

###########----------------------- 12. Save R Session ---------------

save.image(file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_session.RData"))

###########----------------------- 13. Session Information Export ---------------

# Export computational environment details for reproducibility
session_info_output <- sessionInfo()

sink(paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_SessionInfo.txt"))
print(session_info_output)
sink()

saveRDS(session_info_output, 
        file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_SessionInfo.rds"))

loaded_packages <- session_info_output$otherPkgs
base_packages <- session_info_output$basePkgs

package_versions_df <- data.frame(
  Package = c(names(loaded_packages), base_packages),
  Version = c(sapply(loaded_packages, function(x) as.character(x$Version)), 
              rep("base", length(base_packages))),
  stringsAsFactors = FALSE
)

write.csv(package_versions_df, 
          file = paste0(date_export, "_", condition_prefix, "_", ml_prefix, "_PackageVersions.csv"),
          row.names = FALSE)

message("Session information exported for reproducibility.")
message("Random Forest validation analysis completed.")


