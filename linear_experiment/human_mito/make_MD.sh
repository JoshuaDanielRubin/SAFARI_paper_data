#!/bin/bash

# Define directories and reference FASTA files
source_dir="alignments"
target_dir="alignments/with_md"
reference_fasta1="H2a2a1.fa"
reference_fasta2="rCRS.fa"

# Check for samtools and parallel
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found. Please install samtools."
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel could not be found. Please install GNU Parallel."
    exit 1
fi

# Create the target directory if it doesn't exist
mkdir -p "$target_dir"

# Function to process BAM files
process_bam() {
    bam_file="$1"
    filename=$(basename "$bam_file")
    new_file="$target_dir/$filename"

    if [ -f "$new_file" ]; then
        echo "$new_file already exists, skipping."
        return 0
    fi

    # Select reference based on filename keywords
    if [[ "$filename" == *"safari"* ]] || [[ "$filename" == *"giraffe"* ]]; then
        reference_fasta="$reference_fasta1"
    else
        reference_fasta="$reference_fasta2"
    fi

    # Temporary file for filtered BAM
    filtered_bam=$(mktemp)

    # Filter and process BAM file
    samtools view -h "$bam_file" | \
    awk -F'\t' 'BEGIN {OFS = FS}
    {
        if ($0 ~ /^@/) {
            print
        } else {
            cigar = $6
            if (cigar !~ /^[0-9]+[ID]/ && cigar !~ /[ID][0-9]+$/) print
        }
    }' | \
    samtools view -Sb - > "$filtered_bam"

    # Attempt to recalibrate MD tag
    if samtools calmd -b "$filtered_bam" "$reference_fasta" > "$new_file"; then
        echo "Processed $bam_file --> $new_file using $reference_fasta"
        rm "$filtered_bam"
        if ! samtools quickcheck "$new_file"; then
            echo "Warning: $new_file might be corrupted or incomplete."
        fi
    else
        echo "Failed processing $bam_file with reference file $reference_fasta."
        rm "$filtered_bam"
    fi
}

export -f process_bam
export source_dir target_dir reference_fasta1 reference_fasta2

# Process BAM files in parallel
find "$source_dir" -name '*.bam' | parallel -j 20 process_bam {}

echo "Processing complete."

