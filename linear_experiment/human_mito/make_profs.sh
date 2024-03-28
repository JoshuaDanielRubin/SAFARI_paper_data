#!/bin/bash

# Base directories
BAM_DIR_BASE="alignments/with_md/vin"
OUTPUT_DIR_BASE="alignments/profs/vin"
mkdir -p "$OUTPUT_DIR_BASE"

# Function to process BAM files
process_bam_file() {
    local bam_file="$1"
    local BAM_DIR="$2"
    local OUTPUT_DIR="$3"

    base_name=$(basename "$bam_file" .bam)
    output_file="$OUTPUT_DIR/$base_name.prof"
    echo $output_file

    if [ -f "$output_file" ]; then
        echo "Output file for $base_name already exists, skipping..."
        return
    fi

    if [[ "$base_name" == *"single"* ]]; then
        ./bam2prof/src/bam2prof -minl 20 -q -single -minq 0 -length 5 "$bam_file" > "$output_file"
    else
        ./bam2prof/src/bam2prof -minl 20 -q -both -minq 0 -length 5 "$bam_file" > "$output_file"
    fi
}

export -f process_bam_file
export BAM_DIR_BASE OUTPUT_DIR_BASE

# Check and iterate over subdirectories excluding 'under_ten', then parallelize
find "$BAM_DIR_BASE" -mindepth 1 -maxdepth 1 -type d ! -name '*under_ten*' | while read sub_dir; do
    sub_dir_name=$(basename "$sub_dir")
    BAM_DIR="$BAM_DIR_BASE/$sub_dir_name/"
    OUTPUT_DIR="$OUTPUT_DIR_BASE/$sub_dir_name"
    mkdir -p "$OUTPUT_DIR"

    # Check if there are BAM files in the directory
    shopt -s nullglob
    bam_files=("$BAM_DIR"/*{safari,giraffe}*.bam)
    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "No BAM files found in $BAM_DIR directory."
        continue
    fi

    # Export BAM_DIR and OUTPUT_DIR for access in the exported function
    export BAM_DIR OUTPUT_DIR

    # Parallel execution
    printf "%s\n" "${bam_files[@]}" | parallel -j 40 --bar process_bam_file {} "$BAM_DIR" "$OUTPUT_DIR"
done

