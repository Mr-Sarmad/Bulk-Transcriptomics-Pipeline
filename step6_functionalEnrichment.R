#!/usr/bin/env Rscript
# Step 6: Functional enrichment (GSVA) - pipeline ready

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(limma)
  library(gplots)
  library(RColorBrewer)
  library(GSVA)
  # We'll try multiple ways to read GMTs; prefer clusterProfiler::read.gmt if available
})

option_list <- list(
  make_option(c("-g","--gmt"), type="character", default="refs/c2.cp.v2023.2.Hs.symbols.gmt",
              help="GMT file path (gene sets) [default %default]"),
  make_option(c("-i","--input"), type="character", default="results/step4_diffGenes.rds",
              help="Input RDS from step4 [default %default]"),
  make_option(c("-p","--prefix"), type="character", default="results/step6",
              help="Output prefix (files: <prefix>_gsva.rds, <prefix>_topPaths.csv, <prefix>_heatmap.png) [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- Load objects from Step 4 ----
step4 <- readRDS(opt$input)
# Expect v.DEGList (voom object) and design/ebFit in step4 object
if(!("v.DEGList" %in% names(step4))) {
  stop("step4 RDS does not contain 'v.DEGList'. Please ensure step4 saved v.DEGList (voom object).")
}
vobj <- step4$v.DEGList
# If design exists in step4, use it; otherwise read studydesign and recreate
design <- NULL
if("fit" %in% names(step4) && !is.null(step4$fit$design)) {
  design <- step4$fit$design
} else {
  # fallback: try to recreate design from studydesign.txt
  if(file.exists("studydesign.txt")) {
    targets <- read_tsv("studydesign.txt")
    group <- factor(targets$group)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
  } else {
    stop("Cannot find design. Provide a step4 RDS with fit$design or a studydesign.txt file.")
  }
}

# ---- Read GMT (try multiple readers) ----
gmt_path <- opt$gmt
if(!file.exists(gmt_path)) stop("GMT file not found at: ", gmt_path)

# prefer clusterProfiler::read.gmt if available (returns data.frame)
gs <- NULL
if(requireNamespace("clusterProfiler", quietly=TRUE)) {
  gs.df <- clusterProfiler::read.gmt(gmt_path)
  # convert to list-of-character vectors: names -> gene vector
  gs <- split(gs.df$gene, gs.df$gs_name)
} else if(requireNamespace("GSEABase", quietly=TRUE)) {
  # try GSEABase::getGmt or readGmt; wrap in tryCatch
  gmt_obj <- tryCatch({
    if(exists("getGmt", where = asNamespace("GSEABase"))) {
      GSEABase::getGmt(gmt_path)
    } else {
      GSEABase::readGmt(gmt_path)
    }
  }, error = function(e) NULL)
  if(!is.null(gmt_obj)) {
    # convert to list
    gs <- lapply(GSEABase::geneIds(gmt_obj), as.character)
    names(gs) <- names(GSEABase::geneIds(gmt_obj))
  } else {
    stop("Could not read GMT using GSEABase.")
  }
} else {
  stop("Install 'clusterProfiler' or 'GSEABase' to read GMT files, or provide GMT in a different format.")
}

# ---- Run GSVA ----
exprMat <- vobj$E        # voom expression matrix (genes x samples)
# Ensure gene identifiers in exprMat match gene ids in gs (likely gene symbols)
gsva.res <- gsva(exprMat, gs, min.sz = 5, max.sz = 500, verbose = TRUE, parallel.sz=1)

# ---- Linear modeling on GSVA scores (limma) ----
fit.gsva <- lmFit(gsva.res, design)
# If contrast exists in step4 (e.g., contrasts), try to reuse it; otherwise assume disease - healthy
contrast.matrix <- NULL
if(!is.null(step4$fit) && !is.null(step4$fit$coefficients)) {
  # Attempt to recreate a simple contrast if design has "disease" and "healthy"
  if("disease" %in% colnames(design) && "healthy" %in% colnames(design)) {
    contrast.matrix <- makeContrasts(infection = disease - healthy, levels = design)
  }
}
if(is.null(contrast.matrix)) {
  # fallback: if only two columns in design, subtract second - first
  if(ncol(design) == 2) {
    contrast.matrix <- makeContrasts(infection = design[,2] - design[,1], levels = design)
  } else {
    # use fit.gsva directly without contrasts
    contrast.matrix <- NULL
  }
}

if(!is.null(contrast.matrix)) {
  fit.gsva2 <- contrasts.fit(fit.gsva, contrast.matrix)
  ebFit.gsva <- eBayes(fit.gsva2)
} else {
  ebFit.gsva <- eBayes(fit.gsva)
}

# ---- Results: top pathways and decideTests ----
topPaths <- topTable(ebFit.gsva, adjust="BH", coef = 1, number = Inf, sort.by = "logFC")
resSets <- decideTests(ebFit.gsva, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
summary_res <- summary(resSets)
diffSets <- gsva.res[resSets[,1] != 0, , drop = FALSE]

# ---- Clustering + heatmap of differential gene sets ----
hr <- hclust(as.dist(1 - cor(t(diffSets), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(diffSets, method = "spearman")), method = "complete")
mycl <- cutree(hr, k = 2)
mycol <- rainbow(length(unique(mycl)), start = 0.1, end = 0.9)
mycol <- mycol[as.vector(mycl)]
myheatcol <- colorRampPalette(colors = c("yellow", "white", "blue"))(100)

png(paste0(opt$prefix, "_GSVA_diffSets_heatmap.png"), width = 1000, height = 800)
heatmap.2(diffSets,
          Rowv = as.dendrogram(hr),
          Colv = NA,
          col = myheatcol, scale = "row",
          density.info = "none", trace = "none",
          cexRow = 0.9, cexCol = 1, margins = c(10, 14))
dev.off()

# ---- Save outputs ----
saveRDS(list(gsva = gsva.res,
             topPaths = topPaths,
             decideTests = resSets,
             diffSets = diffSets,
             clustering = list(hr = hr, hc = hc, cut = mycl)),
        file = paste0(opt$prefix, "_gsva.rds"))

write_csv(as.data.frame(topPaths) %>% rownames_to_column("pathway"), paste0(opt$prefix, "_topPaths.csv"))
write_csv(as.data.frame(summary_res) %>% rownames_to_column("category"), paste0(opt$prefix, "_decideTests_summary.csv"))

message("Step 6 complete. Outputs written with prefix: ", opt$prefix)
