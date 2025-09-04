# ---- Load libraries ----
library(tidyverse)
library(limma)
library(edgeR)
library(DT)
library(gt)
library(plotly)

# ---- Load inputs from Step 2 ----
myDGEList.filtered.norm <- readRDS("results/step2_DGEList_norm.rds")
targets   <- read_tsv("studydesign.txt")
sampleLabels <- targets$sample
group <- factor(targets$group)

# ---- Design matrix and voom transformation ----
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)

# ---- Linear modeling with limma ----
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = disease - healthy, levels = design)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

# ---- Top hits table ----
myTopHits <- topTable(ebFit, adjust = "BH", coef = 1, number = Inf, sort.by = "logFC")
myTopHits.df <- as_tibble(myTopHits, rownames = "geneID")

datatable.topHits <- datatable(myTopHits.df,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          caption = 'Table 1: Differentially expressed genes (all results)',
          options = list(keys = TRUE, searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

# ---- Volcano plot ----
vplot <- ggplot(myTopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "grey", linewidth = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", colour = c("#2C467A", "#BE684D"), linewidth = 1) +
  labs(title = "Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw()

volcano.interactive <- ggplotly(vplot)

# ---- Filter DEGs (adj.P.Val < 0.01 & |logFC| > 2) ----
results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.01, lfc = 2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] != 0, ]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

datatable.diffGenes <- datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 2: Significant DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = 2:ncol(diffGenes.df), digits = 2)

# ---- Save outputs ----
saveRDS(list(
  v.DEGList = v.DEGList.filtered.norm,
  fit = fit,
  ebFit = ebFit,
  myTopHits = myTopHits.df,
  diffGenes = diffGenes.df
), file = "results/step4_diffGenes.rds")

ggsave("results/step4_volcano.png", vplot, width = 7, height = 6)
htmlwidgets::saveWidget(volcano.interactive, "results/step4_volcano_interactive.html")
DT::saveWidget(datatable.topHits, "results/step4_topHits_table.html")
DT::saveWidget(datatable.diffGenes, "results/step4_diffGenes_table.html")

print("Step 4 complete!")
