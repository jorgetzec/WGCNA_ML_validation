###########--------------------------------------------------------------------###########
#           Gene Ontology Enrichment Analysis (topGO)
#          Salt Stress Response in Chlamydomonas reinhardtii

#  Script for: Integrating co-expression network analysis and machine learning to reveal 
#  the regulatory landscape of GPD genes in Chlamydomonas reinhardtii under salinity stress
#
#  This script performs Gene Ontology (GO) enrichment analysis for WGCNA modules
#  using topGO to identify biological processes, molecular functions, and cellular
#  components associated with salt stress response.
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
  "stringr", "scales", "tidyverse"
)

bioc_pkgs <- c("topGO", "GO.db", "biomaRt", "Rgraphviz")

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
library(scales)
library(tidyverse)
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)

###########----------------------- 1. Data Loading ---------------

# Load expression data for background gene set
counts_full_df <- read.delim("Cre_rawCounts.csv", header = TRUE, sep = ",", row.names = 1)
bg_genes <- rownames(counts_full_df)

# Load WGCNA module assignments
modules <- read.csv("salt_gene_modules.csv", header = TRUE, sep = ",")
colnames(modules)[1] <- "id"

unique_modules <- unique(modules$module)
message("Found ", length(unique_modules), " modules: ", paste(unique_modules, collapse = ", "))

###########----------------------- 2. GO Annotation Setup ---------------

# Connect to Phytozome database for GO annotations
message("Setting up GO annotation database...")
db <- useMart(biomart = "phytozome_mart", 
              dataset = "phytozome", 
              host = "https://phytozome-next.jgi.doe.gov")

message("Retrieving GO annotations for ", length(bg_genes), " genes...")
go_ids <- getBM(attributes = c('gene_name1', 'go_id', 'gene_description', 'kegg_id', 'kegg_desc'), 
                filters = 'gene_name_filter', 
                values = bg_genes, 
                mart = db)

# Build gene-to-GO mapping for topGO
gene_2_GO <- as.list(setNames(as.character(go_ids$go_id), as.character(go_ids$gene_name1)))
gene_2_GO <- lapply(gene_2_GO, function(x) unlist(strsplit(x, split = "[,]")))
gene_2_GO <- gene_2_GO[!sapply(gene_2_GO, function(x) any(is.na(x) | x == ""))]

message("GO annotations retrieved for ", length(gene_2_GO), " genes")

###########----------------------- 3. GO Enrichment Function ---------------

perform_go_enrichment <- function(module_name, ontology = "BP", top_nodes = 50) {
  message("Processing module: ", module_name, " (", ontology, ")")
  
  module_genes <- modules$id[modules$module == module_name]
  
  keep <- module_genes %in% names(gene_2_GO)
  candidate_list <- module_genes[keep]
  
  if (length(candidate_list) < 5) {
    message("Warning: Module ", module_name, " has only ", length(candidate_list), " genes with GO annotations")
    return(NULL)
  }
  
  geneList <- factor(as.integer(bg_genes %in% candidate_list))
  names(geneList) <- bg_genes
  
  # Create topGO data object with weight01 algorithm
  GOdata <- new('topGOdata', 
                ontology = ontology,
                allGenes = geneList, 
                nodeSize = 10, 
                annot = annFUN.gene2GO, 
                gene2GO = gene_2_GO)
  
  weight_fisher_result <- runTest(GOdata, 
                                 algorithm = 'weight01', 
                                 statistic = 'fisher')
  
  all_res <- GenTable(GOdata, 
                      weightFisher = weight_fisher_result, 
                      orderBy = 'weightFisher', 
                      topNodes = top_nodes)
  
  # Multiple testing correction
  p.adj <- round(p.adjust(all_res$weightFisher, method = "BH"), digits = 4)
  
  all_res_final <- cbind(all_res, p.adj)
  all_res_final <- all_res_final[order(all_res_final$p.adj), ]
  
  all_res_final$module <- module_name
  all_res_final$ontology <- ontology
  
  return(all_res_final)
}

###########----------------------- 4. GO Enrichment Analysis ---------------

# Analyze BP and MF ontologies for all modules
ontologies <- c("BP", "MF")
all_go_results <- list()

for (module in unique_modules) {
  for (ontology in ontologies) {
    result <- perform_go_enrichment(module, ontology)
    if (!is.null(result)) {
      all_go_results[[paste(module, ontology, sep = "_")]] <- result
    }
  }
}

combined_results <- bind_rows(all_go_results, .id = "module_ontology")

###########----------------------- 5. Filter and Export Results ---------------

# Filter significant results (p.adj < 0.05)
significant_results <- combined_results %>%
  filter(p.adj < 0.05) %>%
  arrange(module, ontology, p.adj)

top5_results <- significant_results %>%
  group_by(module, ontology) %>%
  slice_head(n = 5) %>%
  ungroup()

write.csv(significant_results, 
          paste0(date_export, "_", condition_prefix, "_GO_enrichment_significant.csv"),
          row.names = FALSE)

write.csv(top5_results, 
          paste0(date_export, "_", condition_prefix, "_GO_enrichment_top5.csv"),
          row.names = FALSE)

message("Exported ", nrow(significant_results), " significant GO terms")
message("Exported top 5 results for ", length(unique(top5_results$module_ontology)), " module-ontology combinations")

###########----------------------- 6. Visualization Functions ---------------

create_go_lollipop <- function(data, module_name, ontology) {
  if (nrow(data) == 0) {
    message("No data to plot for ", module_name, " ", ontology)
    return(NULL)
  }
  
  plot_data <- data %>%
    arrange(p.adj) %>%
    slice_head(n = 5) %>%
    mutate(module = factor(module, levels = c("magenta", "black")))
  
  p <- ggplot(plot_data, aes(x = Significant, y = reorder(Term, Significant))) +
    geom_segment(aes(x = 0, xend = Significant, yend = Term), color = "gray50", linewidth = 1) +
    geom_point(aes(size = Significant, fill = p.adj), shape = 21, stroke = 0.7) +
    scale_fill_gradient2(low = "red", high = "white", midpoint = 0.5) +
    scale_size_continuous(range = c(2, 8)) +
    labs(x = "Count", y = "GO Term", size = "Count", fill = "p.adj") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      legend.title = element_text(size = 12, color = "black"),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      strip.text = element_text(size = 16),
      strip.text.y.right = element_text(angle = 0, color = "black", face = "bold"),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
  
  return(p)
}


###########----------------------- 7. Generate Lollipop Plots ---------------

# Focus on GOD genes module: magenta and black
target_modules <- c("magenta", "black")
target_ontologies <- c("MF", "BP")

lollipop_plots <- list()

for (module in target_modules) {
  for (ontology in target_ontologies) {
    module_data <- significant_results %>%
      filter(module == !!module, ontology == !!ontology)
    
    if (nrow(module_data) > 0) {
      p <- create_go_lollipop(module_data, module, ontology)
      if (!is.null(p)) {
        lollipop_plots[[paste(module, ontology, sep = "_")]] <- p
      }
    }
  }
}

for (plot_name in names(lollipop_plots)) {
  ggsave(
    filename = paste0(date_export, "_", condition_prefix, "_GO_lollipop_", plot_name, ".pdf"),
    plot = lollipop_plots[[plot_name]],
    width = 6, height = 2
  )
}

message("Generated ", length(lollipop_plots), " lollipop plots for magenta and black modules (MF and BP)")


###########----------------------- 8. Summary Statistics ---------------

# Generate summary statistics per module-ontology combination
summary_stats <- significant_results %>%
  group_by(module, ontology) %>%
  summarise(
    n_significant = n(),
    min_padj = min(p.adj, na.rm = TRUE),
    max_padj = max(p.adj, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(module, ontology)

write.csv(summary_stats, 
          paste0(date_export, "_", condition_prefix, "_GO_summary_stats.csv"),
          row.names = FALSE)

message("\n=== GO Enrichment Summary ===")
message("Total significant GO terms: ", nrow(significant_results))
message("Modules analyzed: ", length(unique(significant_results$module)))
message("Ontologies analyzed: ", paste(unique(significant_results$ontology), collapse = ", "))

###########----------------------- 9. Session Information Export ---------------

# Export computational environment for reproducibility
session_info_output <- sessionInfo()

sink(paste0(date_export, "_", condition_prefix, "_GO_SessionInfo.txt"))
print(session_info_output)
sink()

saveRDS(session_info_output, 
        file = paste0(date_export, "_", condition_prefix, "_GO_SessionInfo.rds"))

loaded_packages <- session_info_output$otherPkgs
base_packages <- session_info_output$basePkgs

package_versions_df <- data.frame(
  Package = c(names(loaded_packages), base_packages),
  Version = c(sapply(loaded_packages, function(x) as.character(x$Version)), 
              rep("base", length(base_packages))),
  stringsAsFactors = FALSE
)

write.csv(package_versions_df, 
          file = paste0(date_export, "_", condition_prefix, "_GO_PackageVersions.csv"),
          row.names = FALSE)

message("Gene Ontology enrichment analysis completed.")
message("Session information exported for reproducibility.")

