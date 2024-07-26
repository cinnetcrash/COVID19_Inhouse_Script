#!/bin/bash

# Set required paths
illumina_reads_dir="read directory"  # Update this to your Illumina reads directory
trimmed_output_dir="trimmed out"  # Output directory for trimmed files
dehumanized_output_dir="dehumanized "  # Output directory for dehumanized files
ref_fasta="Ref seq path /MN908947.3.no_poly_A.fa"  # Path to the reference genome FASTA file

# Create output directories if they don't exist
mkdir -p "$trimmed_output_dir"
mkdir -p "$dehumanized_output_dir"

# Define cutadapt and readItAndKeep commands
cutadapt_cmd="cutadapt"  # Ensure cutadapt is correctly installed and accessible
readItAndKeep_cmd="readItAndKeep"  # Ensure readItAndKeep is correctly installed and accessible

# Processing Illumina reads
for file in "$illumina_reads_dir"/*_1.fastq.gz; do
  # Assuming paired-end data
  file_r2="${file%_1.fastq.gz}_2.fastq.gz"  # Construct the name of the mate pair file
  sample_name=$(basename "${file%_1.fastq.gz}")
  
  # Output files for trimmed reads
  trimmed_r1="$trimmed_output_dir/${sample_name}_1_trimmed.fastq.gz"
  trimmed_r2="$trimmed_output_dir/${sample_name}_2_trimmed.fastq.gz"
  
  # Trim adapters and 20 bases from both ends for R1 and R2 using cutadapt
  echo "Trimming adapters and 30 bases from both ends for $sample_name"
  $cutadapt_cmd -u 20 -u -20 -o "$trimmed_r1" -p "$trimmed_r2" "$file" "$file_r2"

  # Dehumanize the trimmed reads
  outprefix="$dehumanized_output_dir/${sample_name}"
  echo "Dehumanizing reads for sample: $sample_name"
  $readItAndKeep_cmd --tech illumina \
                     --ref_fasta "$ref_fasta" \
                     --reads1 "$trimmed_r1" \
                     --reads2 "$trimmed_r2" \
                     --outprefix "$outprefix" \
                     --min_map_length 50 \
                     --min_map_length_pc 50.0

done

echo "Illumina processing completed."

