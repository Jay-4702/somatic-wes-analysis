#!/usr/bin/env bash
set -euo pipefail

SWES="${SWES:-$HOME/swes}"
REF_DIR="$SWES/ref"
CHRS="chr6_and_chr17"
FASTA="$REF_DIR/ref_genome.fa"

mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== [1.3] Reference Genome Setup ==="

# 1. Download and Extract Reference FASTA
if [[ -f "$FASTA" ]]; then
    echo "[SKIP] Reference FASTA already exists: $FASTA"
else
    echo "[DOWNLOAD] Fetching reference genome archive..."
    wget -c "http://genomedata.org/pmbio-workshop/references/genome/$CHRS/ref_genome.tar"
    echo "[EXTRACT] Unpacking reference genome..."
    tar -xvf ref_genome.tar
    rm -f ref_genome.tar
    gunzip -f ref_genome.fa.gz 2>/dev/null || true
    echo "[SUCCESS] FASTA extracted."
fi

# 2. Check and build Samtools FASTA index (.fai)
if [[ -f "$FASTA.fai" ]]; then
    echo "[SKIP] Samtools index (.fai) already exists."
else
    echo "[INDEX] Building Samtools FASTA index..."
    samtools faidx "$FASTA"
fi

# 3. Check and build GATK Sequence Dictionary (.dict)
DICT="${FASTA%.fa}.dict"
if [[ -f "$DICT" ]]; then
    echo "[SKIP] GATK Sequence Dictionary (.dict) already exists."
else
    echo "[INDEX] Building GATK Sequence Dictionary..."
    gatk CreateSequenceDictionary -R "$FASTA" -O "$DICT"
fi

# 4. Check and build BWA Index (.bwt)
if [[ -f "$FASTA.bwt" ]]; then
    echo "[SKIP] BWA index files already exist."
else
    echo "[INDEX] Building BWA reference indexes..."
    bwa index "$FASTA"
fi

echo "=== [1.3] Reference Setup Completed Successfully ==="
