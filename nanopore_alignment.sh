#!/bin/bash

# Define variables
REFERENCE_GENOME="ref_gen/reference.fasta"
OUTPUT_PREFIX="consensus_output_batch_1"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_PREFIX}"

# Get a list of unique sample names
SAMPLE_NAMES=$(ls *.fastq | sed 's/.fastq//')

for SAMPLE_NAME in ${SAMPLE_NAMES}; do

    # Adapter trimming with Porechop
    porechop -i "${SAMPLE_NAME}.fastq" -o "${OUTPUT_PREFIX}/${SAMPLE_NAME}_trimmed.fastq"
    
    # Alignment with Minimap2
    minimap2 -ax map-ont $REFERENCE_GENOME "${OUTPUT_PREFIX}/${SAMPLE_NAME}_trimmed.fastq" | samtools sort -o "${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam" -
    
    # Index the BAM file
    samtools index "${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam"
    
    # Generate VCF file
    bcftools mpileup -Ou -f ${REFERENCE_GENOME} "${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam" | bcftools call -mv -Oz -o "${OUTPUT_PREFIX}/${SAMPLE_NAME}_variants.vcf.gz"

    # Index the VCF file
    bcftools index "${OUTPUT_PREFIX}/${SAMPLE_NAME}_variants.vcf.gz"

    # Alignment statistics
    samtools flagstat "${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam" >> "${OUTPUT_PREFIX}/${SAMPLE_NAME}_alignment_stats.txt"
    
    # Pileup to check variants
    samtools mpileup -aa -A -d 0 -Q 20 "${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam" | ivar consensus -p "${OUTPUT_PREFIX}/${SAMPLE_NAME}_consensus.fasta" -t 0.01 -q 8
    
done
