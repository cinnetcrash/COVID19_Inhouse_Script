#!/bin/bash

# Define variables
REFERENCE_GENOME="reference_genome/reference.fasta"
OUTPUT_PREFIX="Consensus_output"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_PREFIX}

# Get a list of unique sample names
SAMPLE_NAMES=$(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//')

bowtie2-build ${REFERENCE_GENOME} reference_genome/reference

# Quality trimming
# Read trimming with Cutadapt (no adapter trimming)
for SAMPLE_NAME in ${SAMPLE_NAMES}; do
  cutadapt -q 20,20 -o ${OUTPUT_PREFIX}/${SAMPLE_NAME}_R1_trimmed.fastq.gz -p ${OUTPUT_PREFIX}/${SAMPLE_NAME}_R2_trimmed.fastq.gz ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz

  # Human reads deletion
  bowtie2 -p 14 -x GRCh38_noalt_as/GRCh38_noalt_as -1 ${OUTPUT_PREFIX}/${SAMPLE_NAME}_R1_trimmed.fastq.gz -2 ${OUTPUT_PREFIX}/${SAMPLE_NAME}_R2_trimmed.fastq.gz --un-conc-gz ${OUTPUT_PREFIX}/${SAMPLE_NAME}_host_removed

  # Alignment with Bowtie2
  bowtie2 -x reference_genome/reference -1 ${OUTPUT_PREFIX}/${SAMPLE_NAME}_host_removed.1.gz -2 ${OUTPUT_PREFIX}/${SAMPLE_NAME}_host_removed.2.gz -p 8 --very-sensitive-local | samtools sort -o ${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam -

  # Index the BAM file
  samtools index ${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam

  # Variant calling
  freebayes -f ${REFERENCE_GENOME} ${OUTPUT_PREFIX}/${SAMPLE_NAME}_sorted.bam > ${OUTPUT_PREFIX}/${SAMPLE_NAME}_raw.vcf

  # Filter variants using bcftools
  bcftools filter -i 'QUAL>20' ${OUTPUT_PREFIX}/${SAMPLE_NAME}_raw.vcf | bgzip > ${OUTPUT_PREFIX}/${SAMPLE_NAME}_filtered.vcf.gz
  bcftools index ${OUTPUT_PREFIX}/${SAMPLE_NAME}_filtered.vcf.gz

  # Consensus calling
  bcftools consensus -f ${REFERENCE_GENOME} ${OUTPUT_PREFIX}/${SAMPLE_NAME}_filtered.vcf.gz > ${OUTPUT_PREFIX}/${SAMPLE_NAME}_consensus.fasta
done

