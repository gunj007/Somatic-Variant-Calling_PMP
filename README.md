# [Phased Methylation Pattern (PMP)](https://github.com/gunj007/Somatic-Variant-Calling_PMP/blob/main/scripts/pmp.ipynb) 
# Somatic-Variant-Calling

### --- Tools to Download ---

- **Quality Check**: FastQC & MultiQC 
- **Alignment**: BWA
- **Variant Calling**: Mutect2
- **Variant Annotation**: Funcotator 
- **Somatic Mutation Filtering**: GATK

### STEPS:

## 1. Sample Collection and Preparation
- Obtain tumor and matched normal samples.
- Perform quality control (QC) on raw sequencing data (FastQC).

## 2. Read Alignment
- Align reads to a reference genome using tools like **BWA**, **HISAT2**, or **STAR**.
- Post-alignment QC with **samtools** or **Picard** (mark duplicates, evaluate mapping quality).

## 3. Base Quality Score Recalibration (BQSR)
- Recalibrate base quality scores using **GATK's BaseRecalibrator**.

## 4. Somatic Variant Calling
- Call somatic variants using tools:
  - **Mutect2** (GATK)

- Filter variants based on:
  - Tumor vs. normal ratio (Variant Allele Fraction, VAF).
  - Variant quality and read depth.

## 5. Post-Variant Calling Processing
- Annotate variants using tool  **Funcotator*
- Apply filtering to remove germline variants and low-confidence calls.

## 6. Validation and Interpretation
- Validate somatic mutations with **Sanger sequencing** or **ddPCR**.
- Visualize mutations using **IGV** or **UCSC Genome Browser**.

### WORKFLOW:

> Script 1: QC (Quality Control)
Script 2: BWA Alignment and Preprocessing
Script 3: Variant Calling and Annotatio


