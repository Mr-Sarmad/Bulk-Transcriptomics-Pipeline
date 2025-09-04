# ---- Load libraries ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

# ---- Load input from Step 1 ----
Txi_gene <- readRDS("results/step1_txi.rds")
targets   <- read_tsv("studydesign.txt")
sampleLabels <- targets$sample

# ---- Your original code ----
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = -1,
                                  names_to = "samples",
                                  values_to = "expression")

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point",
               shape = 95, size = 10, color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm_vals <- cpm(myDGEList)
keepers <- rowSums(cpm_vals > 1) >= 5
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = -1,
                                           names_to = "samples",
                                           values_to = "expression")

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point",
               shape = 95, size = 10, color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = -1,
                                                names_to = "samples",
                                                values_to = "expression")

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point",
               shape = 95, size = 10, color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

final_plot <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

# ---- Save outputs for pipeline ----
saveRDS(myDGEList.filtered.norm, file = "results/step2_DGEList_norm.rds")
ggsave("results/step2_QCplots.png", final_plot, width = 12, height = 6)

print("Step 2 complete!")
