# Load libraries
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# ---- INPUT ----
targets <- read_tsv("studydesign.txt") # read in your study design
path <- file.path(targets$sample, "abundance.tsv") # kallisto outputs

# ---- Transcript-to-gene mapping ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx) %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select(target_id, gene_name)

# ---- Import counts ----
Txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,  # summarise to gene-level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# ---- SAVE output for Snakemake ----
saveRDS(Txi_gene, file = "results/step1_txi.rds")
