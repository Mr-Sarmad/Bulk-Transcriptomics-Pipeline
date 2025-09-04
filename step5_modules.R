# ---- Load libraries ----
library(tidyverse)
library(gplots)
library(RColorBrewer)

# ---- Load inputs from Step 4 ----
step4 <- readRDS("results/step4_diffGenes.rds")
diffGenes.df <- step4$diffGenes

# Convert tibble back to matrix for clustering/heatmap
diffGenes <- diffGenes.df %>%
  column_to_rownames("geneID") %>%
  as.matrix()

# ---- Heatmap settings ----
myheatcolors <- rev(brewer.pal(name = "RdBu", n = 11))

# ---- Clustering ----
clustRows <- hclust(as.dist(1 - cor(t(diffGenes), method = "pearson")), method = "complete")
clustColumns <- hclust(as.dist(1 - cor(diffGenes, method = "spearman")), method = "complete")

# ---- Module assignment ----
module.assign <- cutree(clustRows, k = 2)  # Adjust k if needed
module.color <- rainbow(length(unique(module.assign)), start = 0.1, end = 0.9)
module.color <- module.color[as.vector(module.assign)]

# ---- Heatmap: all DEGs ----
png("results/step5_diffGenes_heatmap.png", width = 1000, height = 800)
heatmap.2(diffGenes, 
          Rowv = as.dendrogram(clustRows), 
          Colv = as.dendrogram(clustColumns),
          RowSideColors = module.color,
          col = myheatcolors, scale = 'row', labRow = NA,
          density.info = "none", trace = "none",
          cexRow = 1, cexCol = 1, margins = c(8, 20))
dev.off()

# ---- Module 2: Upregulated cluster ----
modulePick <- 2
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_up <- hclust(as.dist(1 - cor(t(myModule_up), method = "pearson")), method = "complete")

png("results/step5_module2_up.png", width = 1000, height = 800)
heatmap.2(myModule_up, 
          Rowv = as.dendrogram(hrsub_up), 
          Colv = NA, labRow = NA,
          col = myheatcolors, scale = "row",
          density.info = "none", trace = "none",
          RowSideColors = module.color[module.assign %in% modulePick],
          margins = c(8, 20))
dev.off()

# ---- Module 1: Downregulated cluster ----
modulePick <- 1
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_down <- hclust(as.dist(1 - cor(t(myModule_down), method = "pearson")), method = "complete")

png("results/step5_module1_down.png", width = 1000, height = 800)
heatmap.2(myModule_down, 
          Rowv = as.dendrogram(hrsub_down), 
          Colv = NA, labRow = NA,
          col = myheatcolors, scale = "row",
          density.info = "none", trace = "none",
          RowSideColors = module.color[module.assign %in% modulePick],
          margins = c(8, 20))
dev.off()

# ---- Save outputs ----
saveRDS(list(
  diffGenes = diffGenes,
  module.assign = module.assign,
  module.colors = module.color,
  myModule_up = myModule_up,
  myModule_down = myModule_down
), file = "results/step5_modules.rds")

print("Step 5 complete!")
