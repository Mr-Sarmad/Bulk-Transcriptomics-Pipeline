#!/usr/bin/env bash

# ========================================================
# RNA-seq Mapping with Kallisto
# Supports both single-end and paired-end data
# ========================================================

# Usage:
#   bash mapping.sh index transcriptome.fa
#   bash mapping.sh single index.idx sample.fastq.gz output_dir
#   bash mapping.sh paired index.idx sample_1.fastq.gz sample_2.fastq.gz output_dir
#
# Example:
#   bash mapping.sh index Homo_sapiens.GRCh38.cdna.all.fa
#   bash mapping.sh paired transcripts.idx reads_1.fq.gz reads_2.fq.gz sample1
#   bash mapping.sh single transcripts.idx reads.fq.gz sample2
# ========================================================

set -euo pipefail

MODE=$1

if [[ "$MODE" == "index" ]]; then
    FASTA=$2
    echo "🔹 Building kallisto index from $FASTA..."
    kallisto index -i transcripts.idx "$FASTA"
    echo "✅ Index saved as transcripts.idx"

elif [[ "$MODE" == "single" ]]; then
    INDEX=$2
    READS=$3
    OUTDIR=$4

    mkdir -p "$OUTDIR"
    echo "🔹 Running kallisto quant (single-end)..."
    kallisto quant -i "$INDEX" -o "$OUTDIR" --single -l 200 -s 20 "$READS"
    echo "✅ Quantification complete: $OUTDIR/abundance.tsv"

elif [[ "$MODE" == "paired" ]]; then
    INDEX=$2
    READS1=$3
    READS2=$4
    OUTDIR=$5

    mkdir -p "$OUTDIR"
    echo "🔹 Running kallisto quant (paired-end)..."
    kallisto quant -i "$INDEX" -o "$OUTDIR" "$READS1" "$READS2"
    echo "✅ Quantification complete: $OUTDIR/abundance.tsv"

else
    echo "❌ Unknown mode: $MODE"
    echo "Usage:"
    echo "  bash mapping.sh index transcriptome.fa"
    echo "  bash mapping.sh single index.idx sample.fastq.gz output_dir"
    echo "  bash mapping.sh paired index.idx sample_1.fastq.gz sample_2.fastq.gz output_dir"
    exit 1
fi
