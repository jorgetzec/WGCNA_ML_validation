###########--------------------------------------------------------------------###########
#           Differential Expression Analysis (DESeq2)
#          Salt Stress Response in Chlamydomonas reinhardtii

#  Script for: Integrating co-expression network analysis and machine learning to reveal 
#  the regulatory landscape of GPD genes in Chlamydomonas reinhardtii under salinity stress
#
#  This script performs differential expression analysis to identify genes
#  differentially expressed under salt stress conditions (200 mM NaCl) across 
#  a time-course experiment using DESeq2.
#
#  Dependencies: Requires WGCNA analysis output (salt_gene_modules.csv) for 
#  module assignment integration.
#
#  Authors: Tzec-Interián et al.
#  Date: October 2025
#
###########--------------------------------------------------------------------###########

condition_prefix <- "salt"
date_export <- format(Sys.Date(), "%Y%m%d")

###########----------------------- 0. Package Installation and Loading ---------------

cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "readr", "reshape2",
  "ggplot2", "ggrepel", "patchwork", "gridExtra", "pheatmap",
  "VennDiagram", "UpSetR", "eulerr", "ggvenn", "stringr", "purrr"
)

bioc_pkgs <- c("DESeq2", "EnhancedVolcano")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(pheatmap)
library(stringr)
library(purrr)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)
library(UpSetR)
library(eulerr)
library(ggvenn)

###########----------------------- 1. Data Loading ---------------

# Load raw count data (columns 5-28 contain sample data)
counts_full_df <- read.delim("Cre_rawCounts.csv", header = TRUE, sep = ",", row.names = 1)
countData_raw <- counts_full_df[, 5:28]

###########----------------------- 2. Sample Metadata Preparation ---------------

sample_ids <- colnames(countData_raw)

colData <- data.frame(
  id = sample_ids,
  condition = ifelse(grepl("^ck", sample_ids), "untreated", "treated"),
  replicates = ifelse(grepl("^ck", sample_ids), "control",
                     gsub("salt(\\d+)h_.*", "Salt_\\1h", sample_ids)),
  time = ifelse(grepl("^ck", sample_ids), "0",
                gsub("salt(\\d+)h_.*", "\\1", sample_ids)),
  stringsAsFactors = FALSE
)

colData$time <- as.numeric(colData$time)
colData$replicates <- as.factor(colData$replicates)
colData$condition <- as.factor(colData$condition)

rownames(colData) <- colData$id

if(!all(rownames(colData) == colnames(countData_raw))) {
  stop("Mismatch between colData rownames and countData_raw colnames.")
}

###########----------------------- 3. Data Filtering and Sample Removal ---------------

# Remove low-expressed genes (mean count > 10) and outlier samples
countData_filtered <- countData_raw %>%
  filter(rowMeans(.) > 10)

# Remove outlier samples identified in WGCNA analysis
samples_to_remove <- c("salt2h_1", "salt72h_1")
countData_filtered <- countData_filtered[, !colnames(countData_filtered) %in% samples_to_remove]
colData_filtered <- colData[!rownames(colData) %in% samples_to_remove, ]

if(!all(colnames(countData_filtered) == rownames(colData_filtered))) {
  stop("Mismatch between filtered count data colnames and filtered colData rownames.")
}

###########----------------------- 4. Multi-factor Analysis (Control vs Treatment) ---------------

# DESeq2 Model for multi-factor analysis: control vs salt treatment
dds <- DESeqDataSetFromMatrix(countData = countData_filtered, 
                              colData = colData_filtered, 
                              design = ~ condition)

# Run DESeq analysis
dds$condition <- as.factor(dds$condition) 
dds$replicates <- as.factor(dds$replicates) 
dds <- DESeq(dds)

# Extract results for control vs treatment comparison
# Parameters: alpha = 0.05 (FDR threshold), lfcThreshold = 2 (log2 fold change threshold)
res_condition <- results(dds,
                        contrast = c("condition", "treated", "untreated"), 
                        alpha = 0.05,   
                        lfcThreshold = 2)

res_condition_df <- as.data.frame(res_condition)

###########----------------------- 5. Annotation and Symbol Mapping ---------------

# Load annotation files from proteome data
proteome_annot <- read.csv("Cre_proteome_uniprot.csv", header = TRUE, sep = ",")

# Create gene symbol mapping from proteome annotations
symbol_ndup2 <- data.frame(
  id = proteome_annot$Gene_id,
  symbol = proteome_annot$Gene_names,
  stringsAsFactors = FALSE
)
symbol_ndup2 <- symbol_ndup2[!duplicated(symbol_ndup2$id), ]

# Convert results to dataframe with gene IDs
if (!"id" %in% colnames(res_condition_df)) {
  res_condition_df <- tibble::rownames_to_column(res_condition_df, "id")
}

# Join with symbol information
res_condition_df <- left_join(res_condition_df, symbol_ndup2, by = "id", unmatched = "drop")

# Add expression status column (up/down/none based on lfc and padj thresholds)
res_condition_df <- res_condition_df %>%
  mutate(expression_status = case_when(
    log2FoldChange >= 2 & padj < 0.05 ~ "up",
    log2FoldChange <= -2 & padj < 0.05 ~ "down",
    TRUE ~ NA_character_
  ))

# Load module assignments (output from WGCNA analysis)
modules <- read.csv("salt_gene_modules.csv", header = TRUE, sep = ",")
colnames(modules)[1] <- "id"

res_condition_df <- res_condition_df %>%
  left_join(modules, by = "id")

# Keep Phytozome IDs if symbol is NA
res_condition_df$symbol[is.na(res_condition_df$symbol)] <- res_condition_df$id[is.na(res_condition_df$symbol)]

# Prepare data for visualization
res_condition_df2 <- res_condition_df[,-1]
rownames(res_condition_df2) <- res_condition_df$symbol
res_condition_df2 <- res_condition_df2[,-8]
res_condition_df2 <- res_condition_df2[,-7]

###########----------------------- 6. Volcano Plot Generation ---------------

# Enhanced volcano plot highlighting GPD genes
EnhancedVolcano(res_condition_df2,
                lab = rownames(res_condition_df2),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('GPD1','GPD2','GPD3','GPD4', 'Cre09.g387763'),
                xlab = bquote(~Log[2]~ 'fold change'),
                boxedLabels = TRUE,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

ggsave(filename = paste0(date_export, "_", condition_prefix, "_GPD_DEGvolcano.pdf"), 
       width = 8, height = 8)

###########----------------------- 7. Time-course Contrast Analysis ---------------

# DESeq2 Model for time-course analysis
dds_time <- DESeqDataSetFromMatrix(countData = countData_filtered, 
                                   colData = colData_filtered, 
                                   design = ~ replicates)

dds_time$condition <- as.factor(dds_time$condition) 
dds_time$replicates <- as.factor(dds_time$replicates) 
dds_time <- DESeq(dds_time)

# Extract results for each time point vs control
# Parameters: alpha = 0.05 (FDR threshold), lfcThreshold = 2 (log2 fold change threshold)
time_points <- c("Salt_2h", "Salt_4h", "Salt_8h", "Salt_12h", "Salt_24h", "Salt_48h", "Salt_72h")
results_list <- list()

for (time_point in time_points) {
  res_temp <- results(dds_time, 
                     contrast = c("replicates", time_point, "control"), 
                     alpha = 0.05, 
                     lfcThreshold = 2)
  results_list[[time_point]] <- res_temp
}

###########----------------------- 8. Wide Format Results (Matrix Style) ---------------

# Create dataframes for each time point in wide format
results_dfs <- list()
for (i in seq_along(time_points)) {
  time_point <- time_points[i]
  res_temp <- results_list[[time_point]]
  
  res_df <- as.data.frame(res_temp)
  colnames(res_df) <- c(paste0(time_point, "_baseMean"),
                        paste0(time_point, "_log2FoldChange"),
                        paste0(time_point, "_lfcSE"),
                        paste0(time_point, "_stat"),
                        paste0(time_point, "_pvalue"),
                        paste0(time_point, "_padj"))
  
  res_df <- cbind(id = rownames(res_df), res_df)
  rownames(res_df) <- 1:nrow(res_df)
  results_dfs[[time_point]] <- res_df
}


# Combine all time point results in wide format using base R
DESeq2_results_wide <- Reduce(function(x, y) inner_join(x, y, by = 'id'), results_dfs)
write.csv(DESeq2_results_wide, paste0(date_export, "_", condition_prefix, "_DESEq2_results_wide_format.csv"))

###########----------------------- 9. Long Format Results (Analysis Style) ---------------

# Process results for each time point in long format
# Classify genes as up/down/none based on padj < 0.05 and |log2FoldChange| > 2
procesar_resultados <- function(res, time) {
  res_df <- as.data.frame(res)
  res_df$time <- time
  
  res_df <- res_df %>%
    mutate(expression = case_when(
      padj < 0.05 & log2FoldChange > 2 ~ "up",
      padj < 0.05 & log2FoldChange < -2 ~ "down",
      TRUE ~ "none"
    ))
  
  res_df <- res_df %>%
    select(time, expression, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    rownames_to_column(var = "gene_id")
  
  return(res_df)
}

# Process all time points
time_labels <- c("2h", "4h", "8h", "12h", "24h", "48h", "72h")
resultados_completos <- list()

for (i in seq_along(time_points)) {
  time_point <- time_points[i]
  time_label <- time_labels[i]
  res_temp <- results_list[[time_point]]
  resultados_completos[[i]] <- procesar_resultados(res_temp, time = time_label)
}

# Combine all results in long format
resultados_completos_df <- bind_rows(resultados_completos)

# Change gene_id column name to id
colnames(resultados_completos_df)[colnames(resultados_completos_df) == "gene_id"] <- "id"

# Join with module information (from WGCNA analysis)
resultados_completos_df <- resultados_completos_df %>%
  left_join(modules, by = "id")

write.csv(resultados_completos_df, paste0(date_export, "_", condition_prefix, "_DEGs_long_format.csv"))

###########----------------------- 10. Annotation Integration ---------------

# Load annotation files (KOG, GO, Pfam functional annotations)
annotation <- read.delim("Creinhardtii_281_v5.6.annotation_info.csv", header = TRUE, sep = ",")

annotation <- annotation %>%
  select(id, KOG, GO, Pfam)

# Join wide format results with annotations
DESeq2_results_annotated <- left_join(DESeq2_results_wide, annotation, by = "id", multiple = "all")
DESeq2_results_final <- left_join(DESeq2_results_annotated, symbol_ndup2, by = "id", unmatched = "drop")

# Keep Phytozome ID if symbol is NA
DESeq2_results_final$symbol[is.na(DESeq2_results_final$symbol)] <- DESeq2_results_final$id[is.na(DESeq2_results_final$symbol)]

# Prepare final dataset for visualization
DESeq2_results_final2 <- DESeq2_results_final[,-1]
rownames(DESeq2_results_final2) <- DESeq2_results_final$symbol
DESeq2_results_final2 <- DESeq2_results_final2[,-44]
DESeq2_results_final2 <- DESeq2_results_final2[,-43]

###########----------------------- 11. Time-course Volcano Plots ---------------

# Volcano plot for 72h time point highlighting GPD genes
EnhancedVolcano(DESeq2_results_final2,
                lab = rownames(DESeq2_results_final2),
                x = 'Salt_72h_log2FoldChange',
                y = 'Salt_72h_pvalue',
                selectLab = c('GPD1','GPD2','GPD3','GPD4', 'Cre09.g387763'),
                xlab = bquote(~Log[2]~ 'fold change'),
                boxedLabels = TRUE,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

ggsave(filename = paste0(date_export, "_", condition_prefix, "_GPD_72h_DEGvolcano.pdf"), 
       width = 6, height = 6)


###########----------------------- 12. GPD Genes Heatmap ---------------

# Filter data for GPD genes (glycerol-3-phosphate dehydrogenase family)
GPD_genes_exp <- c("Cre01.g053000", "Cre01.g053150", "Cre12.g511150", "Cre10.g421700", "Cre09.g387763")
gpd_exp_filtered <- resultados_completos_df[resultados_completos_df$id %in% GPD_genes_exp, ]

# Set time factor levels
niveles_tiempo <- c("2h", "4h", "8h", "12h", "24h", "48h", "72h")
gpd_exp_filtered$time <- factor(gpd_exp_filtered$time, levels = niveles_tiempo)

# Create heatmap matrix
matriz_heatmap <- reshape2::dcast(
  gpd_exp_filtered, 
  id ~ time, 
  value.var = "log2FoldChange"
)

# Convert to matrix
rownames(matriz_heatmap) <- matriz_heatmap$id
matriz_heatmap <- as.matrix(matriz_heatmap[, -1])

# Generate heatmap and save using pdf()
pdf(file = paste0(date_export, "_", condition_prefix, "_gpd_exp_heatmap.pdf"), 
    width = 6, height = 4)

pheatmap(
  matriz_heatmap, 
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(min(matriz_heatmap, na.rm = TRUE), 
               max(matriz_heatmap, na.rm = TRUE), 
               length.out = 51)
)

dev.off()

###########----------------------- 13. UpSet Plot Analysis ---------------

# Create binary matrix for upregulated genes
upregulated_genes <- resultados_completos_df[resultados_completos_df$expression == "up", ]
upregulated_matrix <- as.data.frame.matrix(table(upregulated_genes$id, upregulated_genes$time))
colnames(upregulated_matrix) <- paste0("Time_", colnames(upregulated_matrix))

pdf(file = paste0(date_export, "_", condition_prefix, "_DEGs_up_intersection.pdf"), 
    width = 7, height = 6, onefile = TRUE)

upset(upregulated_matrix,
      nsets = 7,
      mainbar.y.label = "Intersections (Upregulated)",
      sets.x.label = "Genes per time (Upregulated)",
      order.by = "freq")

dev.off()

# Create binary matrix for downregulated genes
downregulated_genes <- resultados_completos_df[resultados_completos_df$expression == "down", ]
downregulated_matrix <- as.data.frame.matrix(table(downregulated_genes$id, downregulated_genes$time))
colnames(downregulated_matrix) <- paste0("Time_", colnames(downregulated_matrix))

pdf(file = paste0(date_export, "_", condition_prefix, "_DEGs_down_intersection.pdf"), 
    width = 7, height = 6, onefile = TRUE)

upset(downregulated_matrix,
      nsets = 7,
      mainbar.y.label = "Intersections (Downregulated)",
      sets.x.label = "Genes per time (Downregulated)",
      order.by = "freq")

dev.off()

###########----------------------- 14. Data Export ---------------

write_delim(res_condition_df,
            file = paste0(date_export, "_", condition_prefix, "_DEG_multiFact.csv"),
            delim = ",")

write.csv(DESeq2_results_final, paste0(date_export, "_", condition_prefix, "_DESEq2_results_annotated.csv"))

###########----------------------- 16. Session Information Export ---------------

session_info_output <- sessionInfo()

sink(paste0(date_export, "_", condition_prefix, "_DESEq2_SessionInfo.txt"))
print(session_info_output)
sink()

saveRDS(session_info_output, 
        file = paste0(date_export, "_", condition_prefix, "_DESEq2_SessionInfo.rds"))

loaded_packages <- session_info_output$otherPkgs
base_packages <- session_info_output$basePkgs

package_versions_df <- data.frame(
  Package = c(names(loaded_packages), base_packages),
  Version = c(sapply(loaded_packages, function(x) as.character(x$Version)), 
              rep("base", length(base_packages))),
  stringsAsFactors = FALSE
)

write.csv(package_versions_df, 
          file = paste0(date_export, "_", condition_prefix, "_DESEq2_PackageVersions.csv"),
          row.names = FALSE)

