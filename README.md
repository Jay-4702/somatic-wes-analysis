# Storage-Optimized Somatic WES Pipeline (SWES)

An end-to-end, reproducible pipeline for Somatic Whole Exome Sequencing (WES) variant calling and functional annotation, engineered specifically for resource-constrained local environments (WSL2/Windows). 

This project processes paired Tumor/Normal exome data (HCC1395 breast cancer cell line) to identify somatic mutations, annotate their clinical impact, and perform pathway enrichment analysis.

## 🧬 Pipeline Architecture & Engineering Highlights

Running clinical genomics pipelines on standard hardware usually results in disk space exhaustion and memory crashes. This pipeline circumvents those limitations using advanced data engineering patterns:

1. **Zero-Intermediate Alignment Streaming**
   Instead of writing massive intermediate SAM and BAM files to disk, this pipeline streams `bwa mem` directly into `samtools fixmate`, `sort`, and `markdup` via Unix pipes. This reduces the disk I/O footprint by **>70%**.
   
2. **API Fault Tolerance & Resiliency**
   Functional annotation is handled via a custom Python client querying the Ensembl VEP REST API. To combat frequent 500/503 server overloads on public endpoints, the client implements an exponential backoff/retry loop with a graceful offline `.csv` fallback, ensuring the pipeline never crashes mid-run.

3. **Strict Environment Isolation**
   Dependency conflicts are resolved using isolated Conda environments (`swes-core`, `swes-gatk`, `swes-py`). 

4. **Idempotency & State Checkpointing**
   Every bash script evaluates the output state before execution. If a step is interrupted due to power loss or hardware limits, the pipeline safely resumes exactly where it left off without overwriting existing progress or re-running expensive compute steps.

---

## ⚙️ Workflow Breakdown

The pipeline is modularized into discrete Bash execution steps:

- **`01_setup_and_slice.sh`**: Data ingestion, reference genome indexing, and slicing panel of normals (PoN) to target chromosomes (chr6 & chr17) to fit local memory constraints.
- **`02_setup_environments.sh`**: Resolves and builds Conda environments.
- **`03_align_and_preprocess.sh`**: Streamed BWA-MEM alignment, duplicate marking, and GATK Base Quality Score Recalibration (BQSR).
- **`04_somatic_variant_calling.sh`**: Tumor-normal matched calling using **GATK4 Mutect2** with local de novo assembly for high-confidence SNVs and Indels. *(Note: Strelka2 execution was bypassed to accommodate WSL2 C++ binary compatibility limits).*
- **`05_annotate_variants.sh`**: Native Linux `zcat | awk` extraction of `PASS` variants, followed by fault-tolerant functional annotation via the Ensembl VEP REST API.

---

## 📊 Exploratory Data Analysis & Biological Interpretation

The pipeline concludes with a Python-based Jupyter Notebook (`notebooks/06_Exploratory_Analysis.ipynb`) that extracts biological meaning from the variant calls using `pandas`, `matplotlib`, and `seaborn`.

**Key Visualizations Include:**
* **Variant Consequences**: Distribution of functional impacts (missense, frameshift, synonymous, etc.).
* **Clinical Impact Severity**: Profiling the percentage of HIGH, MODERATE, and LOW impact mutations.
* **Cancer Driver Identification**: Bar charts of the top 10 most frequently mutated genes in the HCC1395 cell line.
* **KEGG Pathway Enrichment**: Integration with the **Enrichr REST API** to map mutated genes to biological pathways (e.g., p53 signaling, PI3K-Akt) using standardized `-log10(P-value)` transformations.

---

## 🛠️ Tech Stack

* **Bioinformatics Tools**: GATK4, BWA, Samtools, Picard
* **Data Engineering**: Bash scripting, Unix pipes, Conda
* **Data Science**: Python 3.11, Pandas, Matplotlib, Seaborn, Jupyter
* **Cloud / APIs**: Ensembl VEP REST API, Enrichr API

---
*Developed by Jaydeep Desai*
