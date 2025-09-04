# 00_libraries.R
# Master libraries file for the RNA-seq pipeline
# Source this file in every script: source("00_libraries.R")

# Core data science
library(tidyverse)       # ggplot2, dplyr, readr, tidyr, etc.

# Differential expression & RNA-seq
library(edgeR)           # DGEList objects, normalization
library(limma)           # linear modeling, avearrays
library(tximport)        # import Kallisto results
library(rhdf5)           # read Kallisto bootstrap HDF5 files

# Annotation & databases
library(ensembldb)       # Ensembl database access
library(EnsDb.Hsapiens.v86) # Human reference database (change to your organism)
library(biomaRt)         # alternative annotation resource
library(msigdbr)         # MSigDB gene sets

# Functional enrichment & GSEA
library(GSEABase)        # GSEA framework
library(Biobase)         # required by Bioconductor packages
library(GSVA)            # Gene Set Variation Analysis
library(gprofiler2)      # g:Profiler enrichment
library(clusterProfiler) # functional enrichment toolkit
library(enrichplot)      # GSEA visualization

# Visualization
library(ggplot2)
library(gplots)          # heatmap.2
library(RColorBrewer)    # palettes
library(gameofthrones)   # fun palettes
library(d3heatmap)       # interactive heatmaps
library(plotly)          # interactive plots
library(cowplot)         # combine ggplots

# Tables & reporting
library(DT)              # interactive HTML tables
library(gt)              # pretty static tables

# Utilities
library(matrixStats)     # row/column statistics
library(datapasta)       # paste/import data snippets
library(beepr)           # fun sound alerts
library(htmlwidgets)