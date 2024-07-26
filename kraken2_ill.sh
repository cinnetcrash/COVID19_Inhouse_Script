#!/bin/bash

# Define directories
ILLUMINA_DIR="/media/cinnet/PortableSSD1/Paper1icinNanoporeVerileri/karsilastirma_plate2_192_fastq_and_fast5_batch2/192_CleanPlex_MiSeq/All_Files_Ill_Nano/Illumina_samples"  # Update with the correct path to your Illumina reads directory
KRAKEN_DB="/home/cinnet/../../media/cinnet/PortableSSD1/kraken2_db/minikraken2_v2_8GB_201904_UPDATE/"  # Update with the correct path to your Kraken2 database
KRAKEN_OUTPUT_DIR="/media/cinnet/PortableSSD1/Paper1icinNanoporeVerileri/karsilastirma_plate2_192_fastq_and_fast5_batch2/192_CleanPlex_MiSeq/All_Files_Ill_Nano"  # Update with the desired output directory for Kraken2 results

# Create output directory if it doesn't exist
mkdir -p "$KRAKEN_OUTPUT_DIR"

# Function to process paired-end Illumina reads with Kraken2
process_illumina_reads() {
  local READ1=$1
  local READ2=$2
  local SAMPLE=$3

  echo "Running Kraken2 for $SAMPLE"

  # Run Kraken2
  kraken2 --db "$KRAKEN_DB" --threads 12 --report "$KRAKEN_OUTPUT_DIR/$SAMPLE.report" --paired "$READ1" "$READ2" --use-names --memory-mapping > "$KRAKEN_OUTPUT_DIR/$SAMPLE.output"

  echo "Kraken2 completed for $SAMPLE"
}

# Find all paired-end reads in the Illumina directory
for READ1 in "$ILLUMINA_DIR"/*_1.fastq.gz; do
  READ2="${READ1/_R1.fastq.gz/_2.fastq.gz}"
  SAMPLE=$(basename "$READ1" _1.fastq.gz)
  
  if [[ -f "$READ2" ]]; then
    process_illumina_reads "$READ1" "$READ2" "$SAMPLE"
  else
    echo "Paired read for $READ1 not found. Skipping."
  fi
done

echo "All samples processed."

