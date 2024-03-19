#!/bin/bash

# Base directories
BAM_DIR_BASE="alignments/with_md/chag"
OUTPUT_DIR_BASE="alignments/profs/chag"
mkdir -p "$OUTPUT_DIR_BASE"

# Check and iterate over subdirectories excluding 'under_ten'
for sub_dir in "$BAM_DIR_BASE"/*/; do
    if [[ "$sub_dir" == *"under_ten"* ]]; then
        continue
    fi

    sub_dir_name=$(basename "$sub_dir")
    BAM_DIR="$BAM_DIR_BASE/$sub_dir_name/"
    OUTPUT_DIR="$OUTPUT_DIR_BASE/$sub_dir_name"
    mkdir -p "$OUTPUT_DIR"

    # Check if there are BAM files in the directory
    shopt -s nullglob
    bam_files=("$BAM_DIR"/*safari*.bam)
    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "No BAM files found in $BAM_DIR directory."
        continue
    fi

    for bam_file in "${bam_files[@]}"; do
        base_name=$(basename "$bam_file" .bam)
        output_file="$OUTPUT_DIR/$base_name.prof"
        echo $output_file

        if [ -f "$output_file" ]; then
            echo "Output file for $base_name already exists, skipping..."
            continue
        fi

        if [[ "$base_name" == *"single"* ]]; then
            ./bam2prof/src/bam2prof -q -single "$bam_file" > "$output_file"
        else
            ./bam2prof/src/bam2prof -q -double "$bam_file" > "$output_file"
        fi
    done
done

