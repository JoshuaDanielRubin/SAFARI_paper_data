#!/bin/bash

# Define directories and reference FASTA files
source_dir="alignments/chag"
target_dir="alignments/with_md/chag"
reference_fasta1="H2a2a1.fa"
reference_fasta2="rCRS.fa"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found. Please install samtools."
    exit 1
fi

# Check for parallel
if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel could not be found. Please install GNU Parallel."
    exit 1
fi

# Index reference FASTA files (if not already indexed)
samtools faidx "$reference_fasta1"
samtools faidx "$reference_fasta2"

# Function to process BAM files
process_bam() {
    bam_file="$1"
    subfolder=$(dirname "$bam_file" | sed "s|^$source_dir/||") # Extract subfolder from the path
    new_subfolder="$target_dir/$subfolder"
    mkdir -p "$new_subfolder" # Create corresponding subfolder in target_dir

    filename=$(basename "$bam_file" .bam)
    md_bam_file="$new_subfolder/${filename}_sorted_md.bam"

    if [ -f "$md_bam_file" ]; then
        echo "$md_bam_file already exists, skipping."
        return 0
    fi

    # Select reference based on filename keywords
    reference_fasta="$reference_fasta2" # Default reference
    if [[ "$filename" == *"safari"* ]] || [[ "$filename" == *"giraffe"* ]]; then
        reference_fasta="$reference_fasta1"
    fi

    # Sort BAM, add MD tag, and index
    samtools sort "$bam_file" | \
    samtools calmd -b - "$reference_fasta" > "$md_bam_file"

    if [ $? -eq 0 ]; then
        echo "Sorted and added MD to $md_bam_file"
    else
        echo "Failed processing $bam_file with reference $reference_fasta."
        [ -f "$md_bam_file" ] && rm "$md_bam_file"
    fi
}

export -f process_bam
export source_dir target_dir reference_fasta1 reference_fasta2

# Process BAM files in parallel, including those in specific subfolders
find "$source_dir" \( -path "*/linear_results" -o -path "*/k*_*" \) -type f -name '*.bam' | parallel -j 40 process_bam {}

echo "Processing complete."

