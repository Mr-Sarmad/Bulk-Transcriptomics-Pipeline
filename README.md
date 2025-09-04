
```markdown
# 🧬 RNA-seq Analysis Pipeline

A **Snakemake-based workflow** for RNA-seq data analysis.  
This pipeline takes you from **raw FASTQ reads → mapping → quantification → differential expression → functional enrichment**, all in one automated process.  
⚡ **Why this pipeline?**  
It is designed to be **computationally friendly**, meaning it can run efficiently on **minimal or low-spec systems** without requiring large HPC clusters. Perfect for lightweight setups or personal machines.
🚀 **Key Features**
- Supports **single-end** and **paired-end** reads.  
- Works with **Kallisto** for alignment and **tximport** for quantification.  
- Performs **DEG analysis** with DESeq2/limma.  
- Includes **PCA, clustering, and enrichment analysis**.  
- Lightweight – designed to run on **minimal or low-spec machines**.  

---

## ⚙️ Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/rnaseq-pipeline.git
cd rnaseq-pipeline
````

### 1️⃣ Create Conda environment

```bash
conda env create -f environment.yaml
conda activate rnaseq
```

### 2️⃣ Install Python dependencies

```bash
pip install -r requirements.txt
```

---

## 📝 Configuration

The pipeline is controlled through `config.yaml`:

```yaml
samples:
  - HS01
  - HS02
  - HS03
  - HS04
  - HS05
  - CL08
  - CL10
  - CL11
  - CL12
  - CL13

reference:
  fasta: refs/genome.fa
  gtf: refs/annotation.gtf

design: studydesign.txt
```

The experimental design (`studydesign.txt`) should look like:

```text
sample	condition
HS01	Treated
HS02	Treated
HS03	Treated
HS04	Control
HS05	Control
CL08	Control
CL10	Control
CL11	Control
CL12	Control
CL13	Control
```

## 🚀 Running the Pipeline

### 1. **Build genome index**

```bash
bash mapping.sh index refs/genome.fa refs/genome_index
```

### 2. **Run the workflow**

```bash
snakemake --cores 8
```

### 3. **Visualize the DAG**

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

---

## 📊 Outputs

* `results/step1_txi.rds` → gene-level counts
* `results/step3_PCA.png` → PCA plot
* `results/step4_diffGenes_table.html` → DEG table
* `results/step5_heatmap_all.png` → clustering heatmap
* `results/step6_functionalEnrichment_heatmap.png` → enriched pathways

---

## 📜 License

MIT License – free to use and modify with attribution.

---

## 👨‍💻 Author

Developed by **Sarmad Jutt**
📧 Email: [sarmadjutt136@gmail.com](mailto:sarmadjutt136@gmail.com)
📍 Faisalabad, Pakistan

```
