# Storage-Optimized Somatic Whole Exome Sequencing Pipeline (SWES)

A modular, reproducible bioinformatics pipeline for matched tumor-normal whole-exome sequencing (WES), somatic variant calling, functional annotation, and exploratory pathway analysis, designed for resource-constrained local environments using Ubuntu on Windows Subsystem for Linux 2 (WSL2).

**Author:** Jaydeep Desai  
**Project type:** Cancer genomics · Bioinformatics · Pipeline engineering  
**Primary variant caller:** GATK4 Mutect2  
**Annotation:** Ensembl Variant Effect Predictor (VEP) REST API  
**Downstream analysis:** Python, Jupyter Notebook, and Enrichr  

**Scope and interpretation:** The documented Mutect2 analysis was restricted to chromosomes 6 and 17. The reported output contained 272 variants classified as PASS. These results must not be interpreted as a complete exome-wide variant catalogue. The reference genome assembly, upstream alignment implementation, and several quality-control details remain to be verified.

## 1. Project Overview

Somatic mutations contribute to cancer development, progression, and variation in cellular behavior. Whole-exome sequencing provides a way to investigate genetic alterations within protein-coding regions of the genome. Comparing tumor DNA with matched normal DNA helps distinguish candidate somatic alterations from variants also present in the normal sample.

The Storage-Optimized Somatic WES Pipeline (SWES) implements a computational workflow for processing matched tumor-normal sequencing data from the HCC1395 breast cancer cell-line system.

The pipeline combines sequence-processing tools, somatic variant calling, functional annotation, and exploratory biological analysis. Its engineering design addresses practical constraints encountered when running genomics software on a local Windows computer through WSL2, including storage consumption, software dependency conflicts, interrupted execution, and unreliable external API services.

The project has two complementary objectives:
1. Biological analysis: Identify candidate somatic single-nucleotide variants (SNVs) and insertions/deletions (indels), annotate their potential functional consequences, and explore their distribution across genes and biological pathways.
2. Pipeline engineering: Organize the analysis into modular scripts, isolate software dependencies, reduce unnecessary intermediate files where supported by the implementation, and improve recovery from failed or interrupted steps.

The current implementation should be regarded as a research and exploratory analysis workflow. It is not a clinically validated diagnostic pipeline.

## 2. Scientific Background and Motivation

### 2.1 Why matched tumor-normal sequencing?
Tumor samples contain a mixture of genetic alterations, including potentially acquired somatic mutations and inherited germline variants. Analyzing tumor sequencing data alongside a matched normal sample provides additional evidence for distinguishing these categories.

The pipeline uses GATK4 Mutect2 for tumor-normal variant calling. Mutect2 evaluates sequencing evidence to identify candidate somatic variants, including SNVs and indels.

However, matched-normal calling does not eliminate all false positives. Sequencing artifacts, contamination, low coverage, mapping errors, sample mismatches, and other technical issues can still affect the results. Appropriate filtering, quality control, and validation remain essential.

### 2.2 Why optimize for local computing?
Whole-exome sequencing workflows can create large alignment files and require substantial memory, storage, and computational time. On resource-constrained systems, unnecessary intermediate files and poorly managed dependencies can prevent an otherwise valid workflow from completing.

SWES addresses these practical issues through a modular Bash-based architecture, isolated Conda environments, resource-specific reference preparation, and an annotation client designed to tolerate transient API failures.

### 2.3 Project objectives
The main objectives are to:
- Process matched tumor-normal sequencing data.
- Prepare reference resources for somatic variant calling.
- Align sequencing reads and process alignment files.
- Identify candidate somatic variants using Mutect2.
- Annotate variant consequences through Ensembl VEP.
- Summarize variant consequences and predicted impact categories.
- Explore frequently mutated genes and pathway enrichment.
- Improve operational reliability on a local WSL2 environment.
- Maintain a reproducible record of tools, inputs, outputs, and limitations.

## 3. Pipeline Architecture

The documented workflow is organized into modular execution stages.

```mermaid
flowchart TD
    A[Matched tumor-normal FASTQ data] --> B[Reference and resource preparation]
    B --> C[Conda environment setup]
    C --> D[Read alignment and preprocessing]
    D --> E[Mutect2 somatic variant calling]
    E --> F[VCF filtering and PASS extraction]
    F --> G[Ensembl VEP REST API annotation]
    G --> H[Annotated variant data]
    H --> I[Jupyter Notebook exploratory analysis]
    I --> J[Variant consequence summaries]
    I --> K[Gene-level exploration]
    I --> L[Enrichr pathway enrichment]
