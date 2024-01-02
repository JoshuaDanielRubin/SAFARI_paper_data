#!/bin/bash

# Define the directories and reference FASTA files
source_dir="alignments"
target_dir="alignments/with_md"
reference_fasta1="simulations/gen_0.fa"
reference_fasta2="rCRS.fa"

# Ensure samtools is available
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found. Please install samtools."
    exit
fi

# Create the target directory if it doesn't exist
mkdir -p "$target_dir"

# Iterate through all BAM files in the source directory
for bam_file in "$source_dir"/*.bam; do
    # Extract the filename without the directory
    filename=$(basename "$bam_file")

    # Construct the new file path in the target directory
    new_file="$target_dir/$filename"

    # Check if the processed file already exists
    if [ -f "$new_file" ]; then
        echo "$new_file already exists, skipping."
        continue
    fi

    # Determine the correct reference file by inspecting the header of the BAM file
    if samtools view -H "$bam_file" | grep -q "SN:generation_0"; then
        reference_fasta="$reference_fasta1"
    else
        reference_fasta="$reference_fasta2"
    fi

    # Create a temporary filtered BAM file
    filtered_bam=$(mktemp)

    # Filter out reads where CIGAR indicates an indel at the start or end
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

    # Recalculate the MD tags using the filtered BAM and save the new BAM file
    if ! samtools calmd -b "$filtered_bam" "$reference_fasta" > "$new_file"
    then
        # Output an error message if the command fails
        echo "Error processing $bam_file with reference file."
        rm "$filtered_bam"  # Remove the temporary file
        continue  # Skip to the next iteration
    fi

    rm "$filtered_bam"  # Remove the temporary file after use

    # Output progress
    echo "Processed $bam_file --> $new_file"
done

echo "Processing complete."

