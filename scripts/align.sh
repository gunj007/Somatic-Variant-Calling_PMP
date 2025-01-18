#!/bin/bash

# Load required modules or activate conda environment
conda activate samtools

# Define variables
TRIMMED_DATA_DIR=~/bio/projectVC/trim
SAMTOOLS_DIR=~/bio/projectVC/sambam
HG38_DIR=~/bio/projectVC/hg38
THREADS=4

# Create directories if they don't exist
mkdir -p "$SAMTOOLS_DIR"

# Define input and output files
TUMOR_R1="$TRIMMED_DATA_DIR/tumor_R1.fastq.gz"
TUMOR_R2="$TRIMMED_DATA_DIR/tumor_R2.fastq.gz"
NORMAL_R1="$TRIMMED_DATA_DIR/normal_R1.fastq.gz"
NORMAL_R2="$TRIMMED_DATA_DIR/normal_R2.fastq.gz"
TUMOR_SAM="$SAMTOOLS_DIR/tumor.sam"
NORMAL_SAM="$SAMTOOLS_DIR/normal.sam"

# Perform alignment with BWA
bwa mem -t "$THREADS" -R '@RG\tID:tumor\tSM:tumor\tPL:Illumina' "$HG38_DIR/hg38.fa" "$TUMOR_R1" "$TUMOR_R2" > "$TUMOR_SAM"
bwa mem -t "$THREADS" -R '@RG\tID:normal\tSM:normal\tPL:Illumina' "$HG38_DIR/hg38.fa" "$NORMAL_R1" "$NORMAL_R2" > "$NORMAL_SAM"

# Convert SAM to BAM, sort, and mark duplicates
samtools view -Sb "$TUMOR_SAM" | samtools sort -o "$SAMTOOLS_DIR/tumor_sorted.bam"
samtools view -Sb "$NORMAL_SAM" | samtools sort -o "$SAMTOOLS_DIR/normal_sorted.bam"

# Generate read mapping statistics
samtools flagstat "$SAMTOOLS_DIR/normal_sorted.bam" > "$SAMTOOLS_DIR/readmappednormal.txt"
samtools flagstat "$SAMTOOLS_DIR/tumor_sorted.bam" > "$SAMTOOLS_DIR/readmappedtumor.txt"

# Mark duplicates
gatk MarkDuplicates -I "$SAMTOOLS_DIR/tumor_sorted.bam" -O "$SAMTOOLS_DIR/tumor_dedup.bam" -M "$SAMTOOLS_DIR/tumor.metrics.txt"
gatk MarkDuplicates -I "$SAMTOOLS_DIR/normal_sorted.bam" -O "$SAMTOOLS_DIR/normal_dedup.bam" -M "$SAMTOOLS_DIR/normal.metrics.txt"

echo "Alignment, preprocessing, and read mapping statistics completed!"
