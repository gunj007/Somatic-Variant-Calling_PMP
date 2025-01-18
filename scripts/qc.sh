#!/bin/bash

# Load required modules or activate conda environment
conda activate samtools

# Define variables
RAW_DATA_DIR=~/bio/projectVC/raw_data
TRIMMED_DATA_DIR=~/bio/projectVC/trim
THREADS=4

# Create directories if they don't exist
mkdir -p "$TRIMMED_DATA_DIR"

# Perform QC and adapter trimming using fastp
for file in "$RAW_DATA_DIR"/*R1_001.fastq.gz; do
  BASE_NAME=$(basename "$file" R1_001.fastq.gz)
  R1="$RAW_DATA_DIR/${BASE_NAME}R1_001.fastq.gz"
  R2="$RAW_DATA_DIR/${BASE_NAME}R2_001.fastq.gz"
  R1_OUT="$TRIMMED_DATA_DIR/${BASE_NAME}R1.fastq.gz"
  R2_OUT="$TRIMMED_DATA_DIR/${BASE_NAME}R2.fastq.gz"

  echo "Processing $BASE_NAME"
  fastp -i "$R1" -I "$R2" -o "$R1_OUT" -O "$R2_OUT" --detect_adapter_for_pe -c -D --thread "$THREADS"
done

echo "QC and adapter trimming completed!"
