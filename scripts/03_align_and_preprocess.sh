#!/usr/bin/env bash
set -euo pipefail

SWES="${SWES:-$HOME/swes}"
THREADS="${THREADS:-4}"
REF="$SWES/ref/ref_genome.fa"

# 1. Convert targets file to a standard BED format if needed
TARGETS_BED="$SWES/ref/targets.chr6_17.bed"
if [[ ! -f "$TARGETS_BED" ]]; then
    echo "[SETUP] Creating clean BED targets file for GATK..."
    grep -v "^@" "$SWES/ref/targets.chr6_17.interval_list" \
        | awk 'BEGIN{OFS="\t"} NF>=3 {print $1, $2-1, $3}' > "$TARGETS_BED"
fi

mkdir -p "$SWES/align" "$SWES/logs" "$SWES/tmp"

echo "=== [Phase 3] Alignment & Streamed Preprocessing ==="

# ---------------------------------------------------------
# 3.1 Read Quality Control (fastp)
# ---------------------------------------------------------
echo "--- Step 3.1: Running Fastp Read QC ---"
for S in Norm Tumor; do
    JSON_OUT="$SWES/logs/fastp_${S}.json"
    if [[ -f "$JSON_OUT" ]]; then
        echo "[SKIP] Fastp QC report already exists for Exome_${S}. Reusing existing results."
    else
        echo "[QC] Processing Exome_${S} reads via fastp..."
        conda run -n swes-core fastp \
            -i "$SWES/fastq/Exome_${S}/Exome_${S}_R1.fastq.gz" \
            -I "$SWES/fastq/Exome_${S}/Exome_${S}_R2.fastq.gz" \
            -o /dev/null -O /dev/zero \
            --json "$JSON_OUT" --html "$SWES/logs/fastp_${S}.html" \
            --thread "$THREADS" --disable_adapter_trimming
        echo "[SUCCESS] Saved Fastp report for Exome_${S}."
    fi
done

# ---------------------------------------------------------
# 3.2 Zero-Intermediate Alignment Pipe
# ---------------------------------------------------------
align_sample() {
    local TAG=$1 SM=$2 PU=$3 LB=$4
    local OUT_BAM="$SWES/align/${TAG}.dedup.bam"
    local FINAL_BAM="$SWES/align/${TAG}.final.bam"
    local RG="@RG\tID:${TAG}\tPL:ILLUMINA\tPU:${PU}\tLB:${LB}\tSM:${SM}"

    # Check if we already have the final BQSR bam OR the intermediate dedup bam
    if [[ -f "$FINAL_BAM" && -f "${FINAL_BAM}.bai" ]]; then
        echo "[SKIP] Final processed BAM already exists for ${TAG}. Reusing existing results."
        return 0
    elif [[ -f "$OUT_BAM" && -f "${OUT_BAM}.bai" ]]; then
        echo "[SKIP] Deduplicated BAM already exists for ${TAG} ($OUT_BAM). Reusing existing results."
        return 0
    fi

    echo "[ALIGN] Streaming BWA-MEM -> Fixmate -> Sort -> MarkDup for ${TAG}..."
    conda run -n swes-core bash -c "
        bwa mem -t $THREADS -Y -K 100000000 -R '$RG' \
            '$REF' \
            '$SWES/fastq/Exome_${TAG#Exome_}/${TAG}_R1.fastq.gz' \
            '$SWES/fastq/Exome_${TAG#Exome_}/${TAG}_R2.fastq.gz' \
        | samtools fixmate -u -m -@ 2 - - \
        | samtools sort -u -@ 2 -m 768M -T '$SWES/tmp/${TAG}.sort' - \
        | samtools markdup -@ 2 --reference '$REF' \
            -f '$SWES/logs/${TAG}.markdup.txt' - '$OUT_BAM'
        samtools index -@ $THREADS '$OUT_BAM'
    "
    echo "  -> Saved: $OUT_BAM ($(du -h "$OUT_BAM" | cut -f1))"
}

echo "--- Step 3.2: Executing Streamed Alignments ---"
align_sample Exome_Norm  HCC1395BL_DNA C1TD1ACXX-CGATGT.7 exome_norm_lib1
align_sample Exome_Tumor HCC1395_DNA   C1TD1ACXX-ATCACG.7 exome_tumor_lib1

# ---------------------------------------------------------
# 3.3 Base Quality Score Recalibration (BQSR)
# ---------------------------------------------------------
run_bqsr() {
    local TAG=$1
    local IN_BAM="$SWES/align/${TAG}.dedup.bam"
    local RECAL_TABLE="$SWES/align/${TAG}.bqsr.table"
    local FINAL_BAM="$SWES/align/${TAG}.final.bam"

    if [[ -f "$FINAL_BAM" && -f "${FINAL_BAM}.bai" ]]; then
        echo "[SKIP] Final BQSR BAM already exists for ${TAG} ($FINAL_BAM). Reusing existing results."
        return 0
    fi

    if [[ -f "$RECAL_TABLE" ]]; then
        echo "[SKIP] BQSR recalibration table already exists for ${TAG} ($RECAL_TABLE). Reusing table."
    else
        echo "[BQSR] Step 1/2: Calculating error recalibration table for ${TAG}..."
        conda run -n swes-gatk gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" BaseRecalibrator \
            -R "$REF" \
            -I "$IN_BAM" \
            -L "$TARGETS_BED" --interval-padding 100 \
            --known-sites "$SWES/ref/gnomad.chr6_17.vcf.gz" \
            --known-sites "$SWES/ref/mills.chr6_17.vcf.gz" \
            --known-sites "$SWES/ref/known_indels.chr6_17.vcf.gz" \
            -O "$RECAL_TABLE"
    fi

    echo "[BQSR] Step 2/2: Applying BQSR and quality quantization for ${TAG}..."
    conda run -n swes-gatk gatk --java-options "-Xmx4g" ApplyBQSR \
        -R "$REF" \
        -I "$IN_BAM" \
        --bqsr-recal-file "$RECAL_TABLE" \
        -L "$TARGETS_BED" --interval-padding 100 \
        --preserve-qscores-less-than 6 \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        -O "$FINAL_BAM"

    echo "[CLEANUP] Removing intermediate deduplicated BAM for ${TAG} to save storage..."
    rm -f "$IN_BAM" "${IN_BAM}.bai"
    echo "  -> Final BAM ready: $FINAL_BAM ($(du -h "$FINAL_BAM" | cut -f1))"
}

echo "--- Step 3.3: Running BQSR ---"
run_bqsr Exome_Norm
run_bqsr Exome_Tumor

# ---------------------------------------------------------
# 3.4 Alignment QC & Coverage Verification
# ---------------------------------------------------------
echo "--- Step 3.4: Generating Alignment QC Metrics ---"
for TAG in Exome_Norm Exome_Tumor; do
    FLAGSTAT_OUT="$SWES/logs/${TAG}.flagstat"
    MOSDEPTH_OUT="$SWES/logs/${TAG}.mosdepth.summary.txt"
    
    if [[ -f "$FLAGSTAT_OUT" && -f "$MOSDEPTH_OUT" ]]; then
        echo "[SKIP] Alignment QC metrics already exist for ${TAG}. Reusing existing results."
    else
        echo "[QC] Running samtools flagstat & mosdepth for ${TAG}..."
        conda run -n swes-core samtools flagstat -@ "$THREADS" "$SWES/align/${TAG}.final.bam" > "$FLAGSTAT_OUT"
        conda run -n swes-core mosdepth --by "$TARGETS_BED" --no-per-base --fast-mode \
            --thresholds 1,10,20,50,100 "$SWES/logs/${TAG}" "$SWES/align/${TAG}.final.bam"
    fi
done

echo "[MULTIQC] Compiling aggregated QC report..."
conda run -n swes-core multiqc -o "$SWES/logs" -n qc_report "$SWES/logs/" --force 2>/dev/null || true

echo "=== [Phase 3] Alignment & Preprocessing Completed Successfully ==="
