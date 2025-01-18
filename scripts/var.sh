#!/bin/bash

# Load required modules or activate conda environment
conda activate samtools

# Define variables
SAMTOOLS_DIR=~/bio/projectVC/sambam
HG38_DIR=~/bio/projectVC/hg38
MUTECT_DIR=~/bio/projectVC/mutect
FUNCOTATOR_DATASOURCES=~/bio/projectVC/funcotator_dataSources.v1.7.20200521s
THREADS=4
JAVA_OPTS="-Xmx50G"

# Create directories if they don't exist
mkdir -p "$MUTECT_DIR"

# Define input and output files
TUMOR_BAM="$SAMTOOLS_DIR/tumor_dedup.bam"
NORMAL_BAM="$SAMTOOLS_DIR/normal_dedup.bam"
PON_VCF="$MUTECT_DIR/norm_pon.vcf.gz"
GNOMAD_VCF="$MUTECT_DIR/af-only-gnomad.hg38.vcf.gz"
EXOME_INTERVALS="$MUTECT_DIR/exome_calling_regions.v1.1.interval_list"

# Perform Mutect2 variant calling
gatk Mutect2 -R "$HG38_DIR/hg38.fa" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    --germline-resource "$GNOMAD_VCF" \
    --panel-of-normals "$PON_VCF" \
    -O "$MUTECT_DIR/somatic_variants_mutect2.vcf.gz" \
    --f1r2-tar-gz "$MUTECT_DIR/somatic.tar.gz"

# Filter Mutect2 calls
gatk FilterMutectCalls -V "$MUTECT_DIR/somatic_variants_mutect2.vcf.gz" \
    -R "$HG38_DIR/hg38.fa" \
    --contamination-table "$MUTECT_DIR/contamination.table" \
    --ob-priors "$MUTECT_DIR/read-orientation-model.tar.gz" \
    -O "$MUTECT_DIR/somatic_variants_filtered_mutect2.vcf"

# Annotate variants with Funcotator
gatk Funcotator -V "$MUTECT_DIR/somatic_variants_filtered_mutect2.vcf" \
    -R "$HG38_DIR/hg38.fa" \
    --ref-version hg38 \
    --data-sources-path "$FUNCOTATOR_DATASOURCES" \
    --output "$MUTECT_DIR/somatic_variants_functotated.vcf" \
    --output-file-format VCF

echo "Variant calling and annotation completed!"
