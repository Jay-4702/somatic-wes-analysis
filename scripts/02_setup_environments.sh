#!/usr/bin/env bash
set -euo pipefail

echo "=== [Phase 2] Conda Environment Isolation Setup ==="

# 1. Ensure mamba is available in the base environment
if ! command -v mamba &> /dev/null; then
    echo "[INSTALL] Installing mamba into base environment..."
    conda install -n base -c conda-forge mamba -y
fi

# Use flexible priority so bioconda and conda-forge libraries can interoperate smoothly
conda config --set channel_priority flexible

# Helper function to create environment idempotently
create_env() {
    local env_name=$1
    shift
    local packages=("$@")

    if conda info --envs | awk '{print $1}' | grep -Fxq "$env_name"; then
        echo "[SKIP] Environment '$env_name' already exists."
    else
        echo "[CREATE] Building environment: '$env_name'..."
        mamba create -y -n "$env_name" -c conda-forge -c bioconda "${packages[@]}"
    fi
}

# 2.1 Core alignment / BAM / VCF / QC utilities (Unpinned fastp & multiqc for glibc compatibility)
create_env "swes-core" \
    "bwa=0.7.18" "samtools=1.21" "bcftools=1.21" "htslib=1.21" \
    "bedtools=2.31.1" "fastp" "mosdepth=0.3.8" "multiqc"

# 2.2 GATK4 (with OpenJDK 17)
create_env "swes-gatk" \
    "gatk4=4.6.1.0" "openjdk=17"

# 2.3 Strelka2 + Manta (Python 2.7 Island)
create_env "swes-strelka" \
    "strelka=2.9.10" "manta=1.6.0" "python=2.7"

# 2.4 Python 3 Data Science & Variant Interpretation Stack
create_env "swes-py" \
    "python=3.11" "numpy" "scipy" "pandas" "pysam" "pyfaidx" \
    "requests" "matplotlib" "seaborn" "scikit-learn"

echo "[CLEANUP] Purging package caches to reclaim disk space..."
mamba clean --all --yes || true
conda clean --all --yes || true
pip cache purge 2>/dev/null || true

echo "=== [Phase 2] All Environments Built Successfully ==="
