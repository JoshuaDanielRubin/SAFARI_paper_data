#!/bin/bash

BAM_DIR="alignments/with_md"
OUTPUT_DIR="alignments/profs"
mkdir -p "$OUTPUT_DIR"

# Check if there are BAM files in the directory
shopt -s nullglob
bam_files=("$BAM_DIR"/*.bam)
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $BAM_DIR directory."
    exit 1
fi

for bam_file in "${bam_files[@]}"; do
    echo "$bam_file"
    base_name=$(basename "$bam_file" .bam)
    output_file="$OUTPUT_DIR/$base_name.prof"
    
    if [ -f "$output_file" ]; then
        echo "Output file for $base_name already exists, skipping..."
        continue
    fi
    
    if [[ "$base_name" == *"single"* ]]; then
        ./bam2prof/src/bam2prof -q -minl 20 -minq 10 -both "$bam_file" > "$output_file"
    else
        ./bam2prof/src/bam2prof -q -minl 20 -minq 10 -both "$bam_file" > "$output_file"
    fi
done

