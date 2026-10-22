#!/usr/bin/env bash
set -euo pipefail

SWES="${SWES:-$HOME/swes}"
REF_DIR="$SWES/ref"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== [1.4] Stream-Slicing GATK Resource VCFs ==="

# Publicly accessible GATK Best Practices GCS mirror
GB="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
# Open Broad reference bucket mirror
BR="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

slice_vcf() {
    local url=$1
    local out=$2
    
    if [[ -f "$out" && -f "$out.tbi" ]]; then
        echo "[SKIP] $out already sliced and indexed."
        return 0
    fi

    echo "[SLICING] Streaming chr6 and chr17 from $(basename "$url")..."
    tabix -h "$url" chr6 > "${out%.gz}.tmp"
    tabix    "$url" chr17 >> "${out%.gz}.tmp"
    
    echo "[COMPRESS] Bgzipping and indexing $out..."
    bgzip -c "${out%.gz}.tmp" > "$out"
    tabix -f -p vcf "$out"
    rm -f "${out%.gz}.tmp"
    
    echo "  -> Saved: $out ($(du -h "$out" | cut -f1))"
}

# 1. gnomAD (already completed - will safely SKIP)
slice_vcf "$GB/af-only-gnomad.hg38.vcf.gz"                      "gnomad.chr6_17.vcf.gz"
# 2. Panel of Normals (already completed - will safely SKIP)
slice_vcf "$GB/1000g_pon.hg38.vcf.gz"                           "pon.chr6_17.vcf.gz"
# 3. ExAC Common (already completed - will safely SKIP)
slice_vcf "$GB/small_exac_common_3.hg38.vcf.gz"                 "exac_common.chr6_17.vcf.gz"
# 4. Mills Gold Standard Indels (Updated Public URL)
slice_vcf "$BR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "mills.chr6_17.vcf.gz"
# 5. Known Indels (Updated Public URL)
slice_vcf "$BR/Homo_sapiens_assembly38.known_indels.vcf.gz"     "known_indels.chr6_17.vcf.gz"

echo "=== [1.4] All Resources Successfully Sliced ==="
