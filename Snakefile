# RNA-seq Snakemake pipeline (full dependency tracking)

configfile: "config.yaml"

# Final outputs
rule all:
    input:
        "results/step6_functionalEnrichment.rds",
        "results/step6_functionalEnrichment_heatmap.png"


# Step 1: Import transcript-to-gene mapping and tximport
rule step1_tximport:
    input:
        design = "studydesign.txt",
        abundance = expand("{sample}/abundance.tsv", sample=config["samples"])
    output:
        "results/step1_txi.rds"
    script:
        "step1_tximport.R"


# Step 2: Data wrangling
rule step2_dataWrangling:
    input:
        txi = "results/step1_txi.rds",
        design = "studydesign.txt"  # also needed here
    output:
        "results/step2_DGEList_norm.rds"
    script:
        "step2_dataWrangling.R"


# Step 3: Multivariate (PCA, fold change table)
rule step3_multivariate:
    input:
        dge = "results/step2_DGEList_norm.rds",
        design = "studydesign.txt"  # also required here
    output:
        rds     = "results/step3_multivariate.rds",
        pca_png = "results/step3_PCA.png",
        pca_html= "results/step3_PCA_interactive.html",
        dt_html = "results/step3_datatable.html"
    script:
        "step3_multivariate.R"


# Step 4: Differential expression
rule step4_diffGenes:
    input:
        dge = "results/step2_DGEList_norm.rds",
        design = "studydesign.txt"  # also required
    output:
        rds     = "results/step4_diffGenes.rds",
        volcano = "results/step4_volcano.png",
        degs    = "results/step4_diffGenes_table.html"
    script:
        "step4_diffGenes.R"


# Step 5: Modules (clustering + heatmaps)
rule step5_modules:
    input:
        degs = "results/step4_diffGenes.rds",
        design = "studydesign.txt"  # also required
    output:
        rds          = "results/step5_modules.rds",
        heatmap_all  = "results/step5_heatmap_all.png",
        heatmap_up   = "results/step5_heatmap_up.png",
        heatmap_down = "results/step5_heatmap_down.png"
    script:
        "step5_modules.R"


# Step 6: Functional enrichment
rule step6_functionalEnrichment:
    input:
        modules = "results/step5_modules.rds",
        gmt = "refs/c2.cp.v2023.2.Hs.symbols.gmt"  # GMT file required
    output:
        rds     = "results/step6_functionalEnrichment.rds",
        heatmap = "results/step6_functionalEnrichment_heatmap.png"
    script:
        "step6_functionalEnrichment.R"
