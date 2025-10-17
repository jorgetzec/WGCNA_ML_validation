###########--------------------------------------------------------------------###########
#           Weighted Gene Co-expression Network Analysis (WGCNA)
#          Salt Stress Response in Chlamydomonas reinhardtii

#  Script for: Integrating co-expression network analysis and machine learning to reveal 
#  the regulatory landscape of GPD genes in Chlamydomonas reinhardtii under salinity stress
#
#  This script performs WGCNA to identify co-expressed gene modules under
#  salt stress conditions (200 mM NaCl) across a time-course experiment.
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
  "ggplot2", "ggrepel", "patchwork", "gridExtra", "pheatmap"
)

bioc_pkgs <- c("DESeq2", "WGCNA"
)

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
library(DESeq2)
library(WGCNA)

allowWGCNAThreads()

###########----------------------- 1. Data Loading ---------------

counts_full_df <- read.delim("Cre_rawCounts.csv", header = TRUE, sep = ",", row.names = 1)
head(counts_full_df, 3); dim(counts_full_df)

countData_raw <- counts_full_df[, 5:28]
head(countData_raw, 3); dim(countData_raw)

###########----------------------- 2. Sample Metadata Preparation ---------------

sample_ids <- colnames(countData_raw)

colData <- data.frame(
  id = sample_ids,
  condition = ifelse(grepl("^ck", sample_ids), "control", "treated"),
  replicates = gsub(".*?(\\d+)$", "\\1", sample_ids),
  time = ifelse(grepl("^ck", sample_ids), "0",
                gsub("salt(\\d+)h_.*", "\\1", sample_ids)),
  stringsAsFactors = FALSE
)

colData$treatment <- paste0(colData$time, "_h")
colData$time <- as.numeric(colData$time)
colData$replicates <- as.integer(colData$replicates)
colData$replicates <- as.factor(colData$replicates)
colData$condition <- as.factor(colData$condition)
colData$treatment <- as.factor(colData$treatment)

rownames(colData) <- colData$id

head(colData); dim(colData)

if(!all(rownames(colData) == colnames(countData_raw))) {
  stop("Mismatch between colData rownames and countData_raw colnames.")
}

###########----------------------- 3. Data Filtering and Variance Stabilization ---------------

# Remove low-expressed genes (mean count > 10) to reduce noise and improve
# statistical power in network construction
countData_filtered <- countData_raw %>%
  filter(rowMeans(.) > 10)
dim(countData_filtered)

# Variance Stabilizing Transformation (VST) for WGCNA
# Design ~1 and blind=TRUE ensure transformation is unsupervised, 
# appropriate for co-expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData_filtered,
                              colData = colData,
                              design = ~ 1)

vsd <- vst(dds, blind = TRUE)
countData_vst <- assay(vsd)
head(countData_vst, 3); dim(countData_vst)

###########----------------------- 4. Sample Outlier Detection and Removal ---------------

# Hierarchical clustering to identify outlier samples
htree_before_removal <- hclust(dist(t(countData_vst)), method = "average")
plot(htree_before_removal, main = "Sample Clustering (VST data, Before Outlier Removal)", xlab = "", sub = "")
abline(h=60, col="red")

# Remove samples identified as outliers based on clustering analysis
# These samples showed abnormal expression profiles compared to biological replicates
samples_to_remove <- c("salt2h_1", "salt72h_1")

countData_vst_filtered <- countData_vst[, !colnames(countData_vst) %in% samples_to_remove]
colData_filtered <- colData[!rownames(colData) %in% samples_to_remove, ]

dim(countData_vst_filtered)
dim(colData_filtered)

if(!all(colnames(countData_vst_filtered) == rownames(colData_filtered))) {
  stop("Mismatch between filtered VST data colnames and filtered colData rownames.")
}

datExpr <- t(countData_vst_filtered)
head(datExpr[,1:5]); dim(datExpr)

htree_after_removal <- hclust(dist(datExpr), method = "average")
plot(htree_after_removal, main = "Sample Clustering (VST data, After Outlier Removal)", xlab = "", sub = "")

###########----------------------- 5. Data Distribution Visualization ---------------

message("Dimensions of raw count data: ", paste(dim(countData_raw), collapse = " x "))
message("Dimensions of VST transformed data (before sample removal): ", paste(dim(countData_vst), collapse = " x "))
message("Dimensions of VST transformed data (after sample removal, for WGCNA): ", paste(dim(countData_vst_filtered), collapse = " x "))
message("Dimensions of transposed expression data for WGCNA (datExpr): ", paste(dim(datExpr), collapse = " x "))

plot_density_raw <- reshape2::melt(countData_raw, variable.name = "treatment", value.name = "reads") %>%
  ggplot(aes(x = reads, group = treatment)) +
  geom_density(aes(fill = treatment, color = treatment), alpha = 0.01) +
  xlab("Raw Reads") + ylab("Density") + theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"))
plot_density_raw

# Density: VST transformed, filtered counts (used for WGCNA)
plot_density_vst_filtered <- reshape2::melt(as.data.frame(countData_vst_filtered), variable.name = "treatment", value.name = "ExprVST") %>%
  ggplot(aes(x = ExprVST, group = treatment)) +
  geom_density(aes(fill = treatment, color = treatment), alpha = 0.01) +
  xlab("VST Transformed Counts (Filtered)") + ylab("Density") + theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"))

plot_density_vst_filtered

plot_samples_distro <- reshape2::melt(as.data.frame(countData_vst_filtered), value.name = "reads") %>%
  ggplot(aes(x = variable, y = reads)) +
  ggdist::stat_slabinterval(side = "both") +
  xlab("Samples per treatment time (h)") +
  ylab("Distribution - Normalized Expression (VST)") + theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"))
plot_samples_distro

plot_figure_panel <- ((plot_density_raw | plot_density_vst_filtered) / plot_samples_distro) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

plot_figure_panel
ggsave(paste0(date_export, "_", condition_prefix, "_density_comparison_3panels.pdf"), plot_figure_panel, width = 9, height = 9)

###########----------------------- 6. Sample Correlation Heatmap ---------------

cor_matrix <- cor(countData_vst_filtered, method = "pearson")

annotation_df <- data.frame(time = colData$treatment)
rownames(annotation_df) <- rownames(colData)

colors_time <- c("0_h"= "#E41A1C","2_h"="#377EB8", "4_h"= "#4DAF4A", "8_h"= "#984EA3","12_h"="#0000EE", "24_h"="#FF7F00", "48_h"="#00EE00", "72_h"="black")

ann_colors <- list(time = colors_time)

hm_corSample <- pheatmap(cor_matrix,
                         annotation_col = annotation_df,
                         annotation_colors = ann_colors,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         display_numbers = FALSE,
                         number_format = "%.2f",
                         fontsize_number = 8,
                         #main = "Correlation between samples (Pearson)",
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         border_color = NA)


###########----------------------- 7. Soft-Thresholding Power Selection ---------------

# Test range of soft-thresholding powers to achieve scale-free topology
# Goal: R² > 0.8 while maintaining reasonable mean connectivity
power_options <- c(c(1:10), seq(from = 12, to = 30, by = 2))

sft <- pickSoftThreshold(datExpr,
                         powerVector = power_options,
                         networkType = "signed hybrid",  # Distinguishes positive/negative correlations
                         verbose = 5)

sft_data <- sft$fitIndices

plot_sft_rsq <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text_repel(size = 4, segment.color = "black", segment.size = 0.5) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = expression('Scale free topology model fit, signed ' * R^2)) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color="black"), axis.title = element_text(color="black"))

plot_sft_meank <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text_repel(size = 4, segment.color = "black", segment.size = 0.5) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color="black"), axis.title = element_text(color="black"))

plot_sft_combined <- grid.arrange(plot_sft_rsq, plot_sft_meank, nrow = 1)

soft_power_val <- sft_data %>%
  filter(SFT.R.sq >= 0.8) %>%
  pull(Power) %>%
  head(1)

if (is.na(soft_power_val)) {
  warning("No power value met the R^2 threshold. Choosing the power that maximizes R^2.")
  soft_power_val <- sft_data$Power[which.max(sft_data$SFT.R.sq)]
}
message(paste("Selected soft power:", soft_power_val)) 

###########----------------------- 8. Network Construction and Module Detection ---------------

# Resolve namespace conflict between stats::cor and WGCNA::cor
temp_cor_holder <- cor
cor <- WGCNA::cor

max_block_size <- ncol(datExpr) + 100

# Construct co-expression network and identify gene modules
# Key parameters:
#   - mergeCutHeight: 0.25 merges modules with >75% eigengene similarity
#   - minModuleSize: 30 genes minimum per module
#   - deepSplit: 4 (sensitive detection for smaller, distinct modules)
bwnet <- blockwiseModules(datExpr,
                          maxBlockSize = max_block_size,
                          networkType = "signed hybrid",
                          TOMType = "signed",
                          power = soft_power_val,
                          mergeCutHeight = 0.25,
                          minModuleSize = 30,
                          deepSplit = 4,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          saveTOMs = TRUE,
                          saveTOMFileBase = paste0(condition_prefix, "_TOM"),
                          minKMEtoStay = 0,
                          verbose = 3)

cor <- temp_cor_holder

module_colors_assigned <- bwnet$colors
table(module_colors_assigned)

module_gene_counts <- as.data.frame(table(module_colors_assigned)) %>%
  setNames(c("Module", "Genes")) %>%
  filter(Module != "grey") %>%
  arrange(desc(Genes))

module_gene_counts$ColorFill <- ifelse(
  as.character(module_gene_counts$Module) %in% standardColors(100),
  as.character(module_gene_counts$Module),
  "grey"
)

plot_module_counts_bar <- ggplot(module_gene_counts, aes(x = reorder(Module, -Genes), y = Genes, fill = ColorFill)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.5) +
  geom_text(aes(label = Genes), vjust = -0.5, size = 3.5) +
  scale_fill_identity(guide = "none") + # Use actual color names for fill
  labs(x = "Module", y = "Gene Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"))

plot_module_counts_bar
ggsave(filename = paste0(date_export, "_", condition_prefix, "_genes_per_module.pdf"), plot = plot_module_counts_bar, width = 8, height = 5) 

plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05,
                    main = NULL)

plotDendroAndColors(bwnet$dendrograms[[1]],
                    bwnet$colors,
                    "Merged Modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05,
                    main = NULL)


###########----------------------- 9. Gene-Module Assignment for Candidate Genes ---------------

gene_module_key <- tibble::enframe(bwnet$colors, name = "Genes", value = "module")

genes_GPD <- tibble(
  SYMBOL = c("GPD2", "GPD3", "GPD5", "GPD4", "GPD1"),
  Genes = c("Cre01.g053000", "Cre01.g053150", "Cre09.g387763", "Cre10.g421700", "Cre12.g511150")
)

genes_GPD_modules <- merge(genes_GPD, gene_module_key, by = "Genes", all.x = TRUE)

print("GPD genes and their modules:")
print(genes_GPD_modules)

###########----------------------- 10. Module Expression Profile Visualization ---------------

module_to_plot <- "magenta"
gene_to_highlight <- "Cre01.g053000"
gene_highlight_color <- "magenta"

genes_in_module_list <- gene_module_key %>%
  filter(module == module_to_plot) %>%
  pull(Genes)

module_expression_data <- countData_vst_filtered[rownames(countData_vst_filtered) %in% genes_in_module_list, , drop = FALSE]
mean_profile_for_module <- colMeans(module_expression_data)
module_expression_data_with_mean <- rbind(module_expression_data,
                                          "MeanProfileInternalID" = mean_profile_for_module)
module_expr_long_df <- as.data.frame(module_expression_data_with_mean) %>%
  rownames_to_column("Gene") %>%
  melt(id.vars = "Gene", variable.name = "id", value.name = "Expression")

module_expr_long_df <- module_expr_long_df %>%
  mutate(
    PlotGroup = dplyr::case_when(
      Gene == "MeanProfileInternalID" ~ "Mean Profile",
      Gene == gene_to_highlight       ~ "Highlighted Gene",
      TRUE                            ~ "Other Genes"
    ),
    PlotGroup = factor(PlotGroup, levels = c("Other Genes", "Highlighted Gene", "Mean Profile"))
  )

sample_order_for_plot <- colnames(countData_vst_filtered)
module_expr_long_df$id <- factor(module_expr_long_df$id, levels = sample_order_for_plot)

color_palette_module_plot <- c(
  "Other Genes" = "grey80",
  "Highlighted Gene" = gene_highlight_color,
  "Mean Profile" = "black"
)

linewidth_palette_module_plot <- c(
  "Other Genes" = 0.3,
  "Highlighted Gene" = 1.0,
  "Mean Profile" = 1.0
)

linetype_palette_module_plot <- c(
  "Other Genes" = "solid",
  "Highlighted Gene" = "solid",
  "Mean Profile" = "dashed"
)

legend_labels_module_plot <- c(
  "Other Genes" = paste0("Other Genes (", module_to_plot, ")"),
  "Highlighted Gene" = gene_to_highlight,
  "Mean Profile" = paste0("Mean Profile (", module_to_plot, ")")
)

plot_module_profile_samples_viz <- ggplot(module_expr_long_df,
                                          aes(x = id, y = Expression, group = Gene)) +
  
  geom_line(data = . %>% filter(PlotGroup == "Other Genes"),
            aes(color = PlotGroup, linewidth = PlotGroup, , linetype = PlotGroup)) +
  
  geom_line(data = . %>% filter(PlotGroup == "Highlighted Gene"),
            aes(color = PlotGroup, linewidth = PlotGroup, linetype = PlotGroup)) +
  
  geom_line(data = . %>% filter(PlotGroup == "Mean Profile"),
            aes(color = PlotGroup, linewidth = PlotGroup, linetype = PlotGroup)) +
  
  scale_color_manual(values = color_palette_module_plot,
                     labels = legend_labels_module_plot,
                     name = "Line Type",
                     breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")) +
  scale_linewidth_manual(values = linewidth_palette_module_plot,
                         labels = legend_labels_module_plot,
                         name = "Line Type",
                         breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")) +
  scale_linetype_manual(values = linetype_palette_module_plot,
                        labels = legend_labels_module_plot,
                        name = "Line Type",
                        breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")) +
  guides(linewidth = "none") +
  labs(x = "Time treatment under 200 mM NaCl (h)", y = "Normalized expression (VST)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(plot_module_profile_samples_viz)

ggsave(filename = paste0(date_export, "_", condition_prefix, "_moduleProfile_Samples_", module_to_plot, ".pdf"),
       plot = plot_module_profile_samples_viz, width = 6, height = 5)


###########----------------------- 11. Module Membership Analysis ---------------

# Module Membership (MM) quantifies how strongly a gene correlates with 
# the module eigengene. High MM values indicate hub genes or "drivers" 
# that may be functionally important within the module.

module_of_interest <- "magenta"

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

module_eigengenes <- bwnet$MEs

modNames <- gsub("^ME", "", names(module_eigengenes))

# Calculate correlation between each gene and module eigengenes
geneModuleMembership <- cor(datExpr, module_eigengenes, use = 'p')
colnames(geneModuleMembership) <- paste0("MM.", modNames)
rownames(geneModuleMembership) <- colnames(datExpr)

MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
colnames(MMPvalue) <- paste0("p.MM.", modNames)

colMM   <- paste0("MM.", module_of_interest)
colPval <- paste0("p.MM.", module_of_interest)

df_plot <- data.frame(Gene   = rownames(geneModuleMembership),
                      MM     = geneModuleMembership[, colMM],
                      Pvalue = MMPvalue[, colPval],
                      stringsAsFactors = FALSE) %>% 
  mutate(negLogP = -log10(Pvalue))

highlightGenes <- c("Cre01.g053000", "Cre01.g053150", "Cre10.g421700", "Cre09.g387763")
p1 <- ggplot(df_plot, aes(x = MM, y = negLogP)) +
  geom_point(alpha = 0.5, shape = 16, color = "grey50") +
  geom_point(data = subset(df_plot, Gene %in% highlightGenes),
             aes(x = MM, y = negLogP), color = "red", size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth=1) +
  labs(
    title = paste("All Genes vs.", module_of_interest, "Module"),
    x     = paste("Module Membership (MM.", module_of_interest, ")", sep=""),
    y     = expression(-log[10](italic(p) * "-value"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

gene_module_assignments <- data.frame(
  Gene = names(bwnet$colors),
  ActualModule = bwnet$colors,
  stringsAsFactors = FALSE
)

df_mod_specific <- df_plot %>%
  inner_join(gene_module_assignments, by = "Gene") %>%
  filter(ActualModule == module_of_interest)
p2 <- ggplot(df_mod_specific, aes(x = MM, y = negLogP)) +
  geom_point(color = module_of_interest, alpha = 0.7, shape = 16) +
  geom_point(data = subset(df_mod_specific, Gene %in% highlightGenes),
             aes(x = MM, y = negLogP), color = "red", size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth=1) +
  labs(
    title = paste("Genes Assigned to", module_of_interest, "Module"),
    x     = paste("Module Membership (MM.", module_of_interest, ")", sep=""),
    y     = expression(-log[10](italic(p) * "-value"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

combined_plot <- p1 + p2 + plot_layout(ncol = 2)
print(combined_plot)

ggsave(
  filename = paste0(date_export, "_", condition_prefix, "_MM_vs_pMM_GPDs_", module_of_interest, ".pdf"),
  width = 8,
  height = 5
)

###########----------------------- 12. Hub Gene Identification ---------------

# Hub genes are highly connected within modules and often represent key
# regulatory elements. This custom function identifies the top 10 hub genes
# per module based on intramodular connectivity.

moduleColors <- module_colors_assigned
chooseTopHubInEachModule_modf <- function (datExpr, colorh, omitColors = "grey", power = 2, type = "signed hybr", top_n = 10, ...) {
  isIndex = FALSE
  modules = names(table(colorh))
  if (!(is.na(omitColors)[1])) 
    modules = modules[!is.element(modules, omitColors)]
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = seq_len(dim(datExpr)[2])
    isIndex = TRUE
  }
  hubs = list()
  for (m in modules) {
    adj = adjacency(datExpr[, colorh == m], power = power, 
                    type = type, ...)
    top_hubs_indices = order(rowSums(adj), decreasing = TRUE)[seq_len(min(top_n, length(rowSums(adj))))]
    top_hubs = colnames(adj)[top_hubs_indices]
    hubs[[m]] = top_hubs
  }
  if (isIndex) {
    hubs = as.numeric(hubs)
    names(hubs) = modules
  }
  return(hubs)
}


as.data.frame(chooseTopHubInEachModule_modf(
  datExpr, 
  moduleColors, 
  omitColors = "grey", 
  power = soft_power_val, 
  type = "signed hybrid"))

###########----------------------- 13. Intramodular Connectivity Analysis ---------------

# Calculate two connectivity measures:
#   - kWithin (kIM): connectivity within the assigned module
#   - kTotal: connectivity across the entire network
# These metrics complement Module Membership for hub gene identification

connectivity <- intramodularConnectivity.fromExpr(datExpr, moduleColors,
                                                  corOptions = "use = 'p'",
                                                  weights = NULL,
                                                  distFnc = "dist", 
                                                  distOptions = "method = 'euclidean'",
                                                  networkType = "signed hybrid",
                                                  power = soft_power_val,
                                                  scaleByMax = FALSE,
                                                  ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                                  getWholeNetworkConnectivity = TRUE)


rownames(connectivity) <- colnames(datExpr)
connectivity$module_color <- moduleColors
connectivity$gene <- row.names(connectivity)
connect <- connectivity[order(connectivity$module_color,-connectivity$kWithin),]

allClusters_kME <- signedKME(datExpr,module_eigengenes)
head(allClusters_kME)

allClusters_kME <- rownames_to_column(allClusters_kME, var = "gene")
allClusters_kME_kWithin <- left_join(connect,allClusters_kME, by="gene") 
allClusters_kME_kWithin %>%
  filter(module_color=="magenta") %>%
  arrange(desc(kWithin))%>%
  select(gene, kTotal, kWithin,module_color, kMEmagenta)%>%
  head(10)

allClusters_kME_kWithin %>%
  filter(module_color=="black") %>%
  arrange(desc(kWithin))%>%
  select(gene, kTotal, kWithin,module_color,kMEblack)%>%
  head(10)

allClusters_kME_kWithin %>%
  filter(gene=="Cre01.g053000") %>%
  select(kWithin,module_color,gene, kMEmagenta)

allClusters_kME_kWithin %>%
  filter(gene=="Cre01.g053150") %>%
  select(kWithin,module_color,gene, kMEblack)


###########----------------------- 14. Connectivity Visualization (kIM vs kME) ---------------

module_to_analyze <- "magenta"

genes_to_highlight_and_label <- list(
  "Cre01.g053000" = "GPD2",
  "Cre09.g410050" = "Cre09.g410050",
  "Cre03.g174400" = "CDO1"
)

stopifnot(exists("allClusters_kME_kWithin"))

kME_column_name <- paste0("kME", module_to_analyze)

stopifnot(kME_column_name %in% colnames(allClusters_kME_kWithin))

connectivity_data_for_module <- allClusters_kME_kWithin %>%
  filter(module_color == module_to_analyze) %>%
  select(gene, kWithin, kME_value = all_of(kME_column_name))

stopifnot(nrow(connectivity_data_for_module) > 0)

highlight_data <- connectivity_data_for_module %>%
  filter(gene %in% names(genes_to_highlight_and_label)) %>%
  mutate(label_to_show = genes_to_highlight_and_label[gene])

message(paste("Top 10 hub genes by kWithin in module:", module_to_analyze))
print(
  connectivity_data_for_module %>%
    arrange(desc(kWithin)) %>%
    select(gene, kWithin, kME_value) %>%
    head(10)
)

highlight_color <- module_to_analyze

plot_kIM_vs_kME <- ggplot(connectivity_data_for_module, aes(x = kWithin, y = kME_value)) +
  geom_point(color = "grey80", size = 2, alpha = 0.7) +
  geom_point(data = highlight_data,
             aes(x = kWithin, y = kME_value),
             color = highlight_color, size = 4) +
  ggrepel::geom_text_repel(data = highlight_data,
                           aes(x = kWithin, y = kME_value, label = label_to_show),
                           size = 3.5,
                           box.padding = 0.5,
                           point.padding = 0.5,
                           segment.color = 'grey50',
                           max.overlaps = Inf) +
  labs(
    y = paste("kME", tools::toTitleCase(module_to_analyze)), 
    x = "kIM (Intramodular Connectivity)") +  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black")
  )

print(plot_kIM_vs_kME)


###########----------------------- 15. Network Export for Visualization ---------------

# Export network files for visualization in Cytoscape and/or NetworkX (Python)
# Topological Overlap Matrix (TOM) quantifies connection strength between genes

annot_file_path <- 'Cre_proteome_uniprot.csv'
cre_annotations_full <- readr::read_csv(annot_file_path, show_col_types = FALSE)
annotation_map_df <- cre_annotations_full %>%
  dplyr::select(Gene_id, Gene_names) %>%
  dplyr::distinct(Gene_id, .keep_all = TRUE)

load(bwnet$TOMFiles[1])
TOM_matrix <- if (inherits(TOM, "dist")) as.matrix(TOM) else TOM
rm(TOM)

rownames(TOM_matrix) <- colnames(datExpr)
colnames(TOM_matrix) <- colnames(datExpr)

current_gene_ids_in_datExpr <- colnames(datExpr)
TOM_aligned_matrix <- TOM_matrix[current_gene_ids_in_datExpr, current_gene_ids_in_datExpr]
module_colors_for_network <- module_colors_assigned[current_gene_ids_in_datExpr]

alt_node_names_for_network <- annotation_map_df$Gene_names[match(current_gene_ids_in_datExpr, annotation_map_df$Gene_id)]
alt_node_names_for_network[is.na(alt_node_names_for_network)] <- current_gene_ids_in_datExpr[is.na(alt_node_names_for_network)]

message("Exporting full network to Cytoscape files...")
# Threshold 0.05 retains strong connections while keeping file size manageable
full_network_threshold <- 0.05

exportNetworkToCytoscape(
  TOM_aligned_matrix,
  edgeFile = paste0(date_export, "_", condition_prefix, "_CytoscapeInput-edges-FULL.txt"),
  nodeFile = paste0(date_export, "_", condition_prefix, "_CytoscapeInput-nodes-FULL.txt"),
  weighted = TRUE,
  threshold = full_network_threshold,
  nodeNames = current_gene_ids_in_datExpr,
  altNodeNames = alt_node_names_for_network,
  nodeAttr = module_colors_for_network
)
message(paste("Full network exported. Edge threshold:", full_network_threshold))

# Export module-specific network for detailed analysis
# Threshold 0.17 selected to capture biologically relevant connections for GPD genes in specific module
# while maintaining network interpretability
module_color_to_export <- "magenta"
module_edge_threshold  <- 0.17

message(paste("\nExporting network for module:", module_color_to_export, "..."))

genes_in_target_module <- current_gene_ids_in_datExpr[module_colors_for_network == module_color_to_export]

module_tom_subset <- TOM_aligned_matrix[genes_in_target_module, genes_in_target_module]

module_alt_node_names_subset <- alt_node_names_for_network[match(genes_in_target_module, current_gene_ids_in_datExpr)]
module_colors_subset <- module_colors_for_network[genes_in_target_module]

cytoscape_output_module <- exportNetworkToCytoscape(
  module_tom_subset,
  edgeFile = paste0(date_export, "_", condition_prefix, "_CytoscapeInput-edges-MODULE_", module_color_to_export, ".txt"),
  nodeFile = paste0(date_export, "_", condition_prefix, "_CytoscapeInput-nodes-MODULE_", module_color_to_export, ".txt"),
  weighted = TRUE,
  threshold = module_edge_threshold,
  nodeNames = genes_in_target_module,
  altNodeNames = module_alt_node_names_subset,
  nodeAttr = module_colors_subset
)
message(paste("Module", module_color_to_export, "network exported. Edge threshold:", module_edge_threshold))

message("\nNetwork export process completed.")

###########----------------------- 17. Data Export for Downstream Analysis ---------------

# Export normalized expression data and module assignments for:
#   - Machine learning validation
#   - Functional enrichment analysis
#   - Integration with other datasets

write.csv(countData_vst_filtered, 
          file = paste0(condition_prefix, "_VST_corrected_fulldata.csv"),
          row.names = TRUE)

gene_modules_df <- data.frame(
  module = module_colors_assigned,
  row.names = names(module_colors_assigned)
)
write.csv(gene_modules_df, 
          file = paste0(condition_prefix, "_gene_modules.csv"),
          row.names = TRUE)

message("Data exported for machine learning analysis.")

###########----------------------- 18. Save R Session ---------------
save.image(file = paste0(date_export, "_", condition_prefix, "_WGCNA_session.RData"))

###########----------------------- 19. Session Information Export ---------------

# Export computational environment details for reproducibility
# Critical for scientific publications and sharing with collaborators
# Includes: R version, OS, loaded packages with versions

session_info_output <- sessionInfo()

sink(paste0(date_export, "_", condition_prefix, "_SessionInfo.txt"))
print(session_info_output)
sink()

saveRDS(session_info_output, 
        file = paste0(date_export, "_", condition_prefix, "_SessionInfo.rds"))

loaded_packages <- session_info_output$otherPkgs
base_packages <- session_info_output$basePkgs

package_versions_df <- data.frame(
  Package = c(names(loaded_packages), base_packages),
  Version = c(sapply(loaded_packages, function(x) as.character(x$Version)), 
              rep("base", length(base_packages))),
  stringsAsFactors = FALSE
)

write.csv(package_versions_df, 
          file = paste0(date_export, "_", condition_prefix, "_PackageVersions.csv"),
          row.names = FALSE)

message("Session information exported for reproducibility.")
message("WGCNA analysis for salt stress completed.")

