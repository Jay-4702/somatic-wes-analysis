# Storage-Optimized Somatic WES Pipeline (SWES)

An end-to-end, reproducible pipeline for Somatic Whole Exome Sequencing (WES) variant calling, engineered specifically for resource-constrained environments.

## Architecture & Engineering Highlights

1. **Zero-Intermediate Alignment Streaming**
   Instead of writing large intermediate SAM and BAM files to disk, this pipeline streams `bwa mem` directly into `samtools fixmate`, `sort`, and `markdup` via Unix pipes. This reduces disk IO footprint by >70%, making it possible to run clinical datasets on local hardware.
   
2. **Idempotency & Checkpointing**
   Every phase of the pipeline evaluates output state before execution. If a step fails due to external factors (e.g., power loss or API timeouts), the pipeline safely resumes exactly where it left off without overwriting existing progress.

3. **Strict Environment Isolation**
   Dependency conflicts (like Python 2.7 requirements for older callers vs Python 3.11 for modern annotation) are resolved using four strict Conda environments (`swes-core`, `swes-gatk`, `swes-strelka`, `swes-py`).

4. **API Fault Tolerance**
   Functional annotation is handled via a custom Python client querying the Ensembl VEP REST API. The client implements an exponential backoff/retry loop to gracefully handle 500/503 server overloads without crashing the pipeline.

## Usage
The pipeline is broken into modular bash scripts:
* `01_3_setup_reference.sh` / `01_4_slice_resources.sh`: Data ingestion and genome indexing.
* `02_setup_environments.sh`: Resolves and builds Conda environments.
* `03_align_and_preprocess.sh`: Streamed alignment and GATK BQSR.
* `04_somatic_variant_calling.sh`: GATK Mutect2 tumor-normal matched calling.
* `05_annotate_variants.sh`: Fault-tolerant VEP REST API annotation.
