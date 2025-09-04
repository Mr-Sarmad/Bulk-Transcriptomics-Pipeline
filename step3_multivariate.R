# ---- Load libraries ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

# ---- Load inputs from Step 2 ----
myDGEList.filtered.norm <- readRDS("results/step2_DGEList_norm.rds")
targets   <- read_tsv("studydesign.txt")
sampleLabels <- targets$sample
group <- factor(targets$group)

# ---- Log2 CPM for PCA ----
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# ---- PCA ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=FALSE, retx=TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) + 
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

pca.interactive <- ggplotly(pca.plot)

# ---- Fold change table ----
mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = (HS01 + HS02 + HS03 + HS04 + HS05)/5, 
                    disease.AVG = (CL08 + CL10 + CL11 + CL12 + CL13)/5,
                    LogFC = (disease.AVG - healthy.AVG)) %>%
  mutate_if(is.numeric, round, 2)

datatable.obj <- datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

# ---- Save outputs ----
saveRDS(list(
  log2.cpm = log2.cpm.filtered.norm,
  log2.df  = log2.cpm.filtered.norm.df,
  pca.res  = pca.res,
  pca.df   = pca.res.df,
  foldchange.table = mydata.df
), file = "results/step3_multivariate.rds")

ggsave("results/step3_PCA.png", pca.plot, width=7, height=6)
htmlwidgets::saveWidget(pca.interactive, "results/step3_PCA_interactive.html")
DT::saveWidget(datatable.obj, "results/step3_datatable.html")

print("Step 3 complete!")
