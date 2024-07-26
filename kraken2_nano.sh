#!/bin/bash

# Define directories
NANOPORE_DIR="/media/cinnet/VERBATIM SD/Paper1/Revisions/03_NanoFilted/01_Nanofilt_Data_NoGuppy"  # Update with the correct path to your Nanopore reads directory
KRAKEN_DB="/home/cinnet/../../media/cinnet/PortableSSD1/kraken2_db/minikraken2_v2_8GB_201904_UPDATE/"  # Update with the correct path to your Kraken2 database
KRAKEN_OUTPUT_DIR="/media/cinnet/VERBATIM SD/Paper1/Revisions/03_NanoFilted/01_Nanofilt_Data_NoGuppy/Nanopore_Kraken2_reports_Noguppy"  # Update with the desired output directory for Kraken2 results

# Create output directory if it doesn't exist
mkdir -p "$KRAKEN_OUTPUT_DIR"

# Function to process Nanopore reads with Kraken2
process_nanopore_reads() {
  local READ=$1
  local SAMPLE=$2

  echo "Running Kraken2 for $SAMPLE"

  # Run Kraken2
  kraken2 --db "$KRAKEN_DB" --threads 12 --report "$KRAKEN_OUTPUT_DIR/$SAMPLE.report" --use-names --memory-mapping "$READ" > "$KRAKEN_OUTPUT_DIR/$SAMPLE.output"

  echo "Kraken2 completed for $SAMPLE"
}

# Find all single-end reads in the Nanopore directory
for READ in "$NANOPORE_DIR"/*.fastq; do
  SAMPLE=$(basename "$READ" .fastq)
  process_nanopore_reads "$READ" "$SAMPLE"
done

echo "All samples processed."

