#!/usr/bin/env bash
set -euo pipefail

SWES="${SWES:-$HOME/swes}"
THREADS="${THREADS:-4}"
REF="$SWES/ref/ref_genome.fa"

mkdir -p "$SWES/calls/mutect2" "$SWES/calls/strelka2" "$SWES/logs"

echo "=== [Phase 4] Somatic Variant Calling ==="

NORM_BAM="$SWES/align/Exome_Norm.final.bam"
TUMOR_BAM="$SWES/align/Exome_Tumor.final.bam"

# ---------------------------------------------------------
# 4.1 GATK Mutect2 Calling
# ---------------------------------------------------------
echo "--- Step 4.1: Running GATK4 Mutect2 ---"
M2_VCF="$SWES/calls/mutect2/mutect2.unfiltered.vcf.gz"
M2_FILTERED="$SWES/calls/mutect2/mutect2.filtered.vcf.gz"

if [[ -f "$M2_FILTERED" && -f "${M2_FILTERED}.tbi" ]]; then
    echo "[SKIP] Mutect2 filtered VCF already exists. Reusing existing results."
else
    if [[ ! -f "$M2_VCF" ]]; then
        echo "[MUTECT2] Calling raw somatic variants on chr6 & chr17..."
        # Explicitly passing contigs instead of an external BED file
        conda run -n swes-gatk gatk --java-options "-Xmx4g" Mutect2 \
            -R "$REF" \
            -I "$TUMOR_BAM" \
            -I "$NORM_BAM" \
            -normal "HCC1395BL_DNA" \
            -L chr6 -L chr17 \
            --germline-resource "$SWES/ref/gnomad.chr6_17.vcf.gz" \
            --panel-of-normals "$SWES/ref/pon.chr6_17.vcf.gz" \
            -O "$M2_VCF"
    fi

    echo "[MUTECT2] Filtering somatic variant calls..."
    conda run -n swes-gatk gatk --java-options "-Xmx4g" FilterMutectCalls \
        -R "$REF" \
        -V "$M2_VCF" \
        -O "$M2_FILTERED"
    echo "[SUCCESS] Mutect2 calling and filtering complete."
fi

