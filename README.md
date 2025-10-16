# WGCNA and Machine Learning Validation Analysis of Salt Stress Response in Chlamydomonas reinhardtii

## Overview

This repository contains the computational workflow for analyzing salt stress response in *Chlamydomonas reinhardtii* using Weighted Gene Co-expression Network Analysis (WGCNA) and machine learning validation. The analysis integrates co-expression network analysis and machine learning to uncover GPD (glycerol-3-phosphate dehydrogenase) gene regulation under salinity stress conditions.

## Research Context

**Study:** Integrating co-expression network analysis and machine learning to uncover GPD genes regulation in *Chlamydomonas reinhardtii* under salinity stress

**Authors:** Tzec-Interián et al.

**Date:** October 2025

**Journal:** PeerJ

## Repository Contents

### Analysis Scripts

| Script | Purpose | Description |
|--------|---------|-------------|
| `Cre_Salt_WGCNA.R` | Main WGCNA Analysis | Identifies co-expression modules from RNA-seq data under salt stress (200 mM NaCl) |
| `Cre_Salt_DESeq2.R` | Differential Expression | Performs differential expression analysis across time-course salt stress conditions |
| `Cre_Salt_GeneOntology.R` | GO Enrichment | Performs Gene Ontology enrichment analysis for WGCNA modules using topGO |
| `Cre_salt_MachineLearning_RandomForest.R` | ML Validation | Validates WGCNA module assignments using Random Forest classification |
| `Cre_Salt_ModulePreservation.R` | Module Preservation | Evaluates module stability against null models through permutation analysis |

### Data Files

| File | Description |
|------|-------------|
| `Cre_rawCounts.csv` | Raw RNA-seq count data for all samples |
| `Cre_proteome_uniprot.csv` | Gene annotations and protein information |

### Session Information

| File | Description |
|------|-------------|
| `20251014_salt_SessionInfo.txt` | R session info from WGCNA analysis |
| `20251014_salt_ML_RandomForest_SessionInfo.txt` | R session info from ML validation |
| `20251014_salt_ModulePreservation_SessionInfo.txt` | R session info from preservation analysis |


## Analysis Steps

### Step 1: WGCNA Analysis (`Cre_Salt_WGCNA.R`)

**Purpose:** Identifies co-expression modules from RNA-seq data under salt stress conditions.

**Key Features:**
- Data quality control and outlier removal
- Variance Stabilizing Transformation (VST)
- Soft threshold selection for scale-free topology
- Dynamic module detection with mergeCutHeight = 0.25
- Module characterization and hub gene identification
- GPD gene connectivity analysis (kIM, kME, hub gene identification)
- Network export for Cytoscape visualization

**Outputs:**
- VST-transformed expression data
- Gene-to-module assignments (`salt_gene_modules.csv`)
- Module eigengenes and connectivity metrics
- Network files for visualization
- Complete R session for reproducibility

### Step 2: Differential Expression Analysis (`Cre_Salt_DESeq2.R`)

**Purpose:** Identifies differentially expressed genes across time-course salt stress conditions.

**Key Features:**
- Multi-factor analysis (control vs treatment)
- Time-course contrast analysis (2h, 4h, 8h, 12h, 24h, 48h, 72h)
- Integration with WGCNA module assignments
- Enhanced volcano plots highlighting GPD genes
- GPD gene expression heatmaps
- Venn diagrams and UpSet plots for gene intersections
- Comprehensive annotation integration

**Dependencies:** Requires WGCNA output (`salt_gene_modules.csv`)

**Outputs:**
- Wide format results (matrix style)
- Long format results (analysis style)
- Annotated differential expression results
- Publication-ready visualizations
- Session information for reproducibility

### Step 3: Gene Ontology Enrichment Analysis (`Cre_Salt_GeneOntology.R`)

**Purpose:** Performs Gene Ontology enrichment analysis for WGCNA modules to identify biological functions.

**Key Features:**
- GO enrichment analysis using topGO for all three ontologies (BP, MF, CC)
- Automatic processing of all WGCNA modules
- Multiple testing correction (Benjamini-Hochberg)
- Top 5 enriched terms per module per ontology
- Publication-ready bar plots for enriched terms
- Integration with biomaRt for GO annotations

**Dependencies:** Requires WGCNA output (`salt_gene_modules.csv`) and raw count data (`Cre_rawCounts.csv`)

**Outputs:**
- Significant GO terms (p.adj < 0.05)
- Top 5 results per module-ontology combination
- Bar plots for enriched terms
- Summary statistics
- Session information for reproducibility

### Step 4: Machine Learning Validation (`Cre_salt_MachineLearning_RandomForest.R`)

**Purpose:** Validates WGCNA module assignments using supervised machine learning.

**Key Features:**
- Independent train/test split (80/20)
- Class balancing with upsampling
- Random Forest with 5-fold cross-validation
- Comprehensive performance metrics (ROC, AUC, confusion matrix)
- UMAP dimensionality reduction visualization
- Misclassification analysis for GPD-relevant modules

**Performance Metrics:**
- Overall accuracy and Kappa statistic
- Per-module AUC (one-vs-rest approach)
- Sensitivity and specificity
- Detailed misclassification reports

### Step 5: Module Preservation Analysis (`Cre_Salt_ModulePreservation.R`)

**Purpose:** Evaluates module stability against null models through permutation analysis.

**Key Features:**
- Null model generation by permuting time points within genes
- Module preservation analysis with 100 permutations
- Z-summary and medianRank statistics
- Module eigengene correlation analysis
- Statistical significance assessment

**Interpretation Guidelines:**
- Z-summary > 10: Strongly preserved
- Z-summary 2-10: Moderately preserved  
- Z-summary < 2: Not preserved

## Installation and Setup

### System Requirements

- **R Version:** 4.0.0 or higher (tested with R 4.4.1)
- **Operating System:** Windows, macOS, or Linux
- **Memory:** Minimum 8GB RAM recommended
- **Storage:** ~2GB for complete analysis

### Required R Packages

**Core Analysis:**
```r
# CRAN packages
install.packages(c("dplyr", "tidyr", "tibble", "readr", "reshape2",
                   "ggplot2", "ggrepel", "patchwork", "gridExtra", "pheatmap",
                   "randomForest", "pROC", "caret", "umap"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "WGCNA"))
```

### Data Preparation

1. Ensure `Cre_rawCounts.csv` and `Cre_proteome_uniprot.csv` are in the working directory
2. Verify data format matches expected structure (genes as rows, samples as columns)

## Running the Analysis

### Complete Workflow

```r
# Step 1: Run WGCNA analysis
source("Cre_Salt_WGCNA.R")

# Step 2: Run differential expression analysis
source("Cre_Salt_DESeq2.R")

# Step 3: Run Gene Ontology enrichment analysis
source("Cre_Salt_GeneOntology.R")

# Step 4: Run Random Forest validation  
source("Cre_salt_MachineLearning_RandomForest.R")

# Step 5: Run module preservation analysis
source("Cre_Salt_ModulePreservation.R")
```

### Individual Scripts

Each script can be run independently if required input files are available:

```r
# WGCNA analysis (requires raw count data)
source("Cre_Salt_WGCNA.R")

# Differential expression analysis (requires WGCNA outputs)
source("Cre_Salt_DESeq2.R")

# Gene Ontology enrichment (requires WGCNA outputs and raw count data)
source("Cre_Salt_GeneOntology.R")

# ML validation (requires WGCNA outputs)
source("Cre_salt_MachineLearning_RandomForest.R")

# Module preservation (requires WGCNA outputs)
source("Cre_Salt_ModulePreservation.R")
```

## Reproducibility

### Random Seeds

All random processes use fixed seeds for reproducibility:
- **WGCNA:** `randomSeed = 1234`
- **ML Training:** `set.seed(123)`
- **Data Splitting:** `set.seed(123)`
- **UMAP:** `set.seed(42)`
- **Module Preservation:** `randomSeed = 1234`

### Session Information

Complete session information is exported for each analysis step, including:
- R version and platform details
- Loaded packages with versions
- System environment details

## Citation

When using this workflow, please cite:

1. **WGCNA:** Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559
2. **Random Forest:** Breiman L (2001) Random Forests. Machine Learning 45:5-32
3. **DESeq2:** Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550


## Contact and Support
For questions about this analysis or repository, please refer to the original publication or contact the corresponding author.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Version History
- **v1.0** (October 2025): Initial release
  - Complete WGCNA analysis pipeline
  - Random Forest validation with comprehensive metrics
  - Module preservation analysis with null models
  - Publication-ready visualizations and documentation

---
**Note:** This repository contains the computational analysis supporting the peer-reviewed publication. All scripts are designed for reproducibility and transparency in scientific research.
